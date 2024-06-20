
# this rule gets a raw h5 file and runs seacells to aggreate cells to metacells
# it returns an h5 file with aggregated counts
rule seacells:
    input:
        counts="results/counts_filtered/{sample}.genes_cells_filtered.h5",
        celltypes="results/celltyping/{sample}.cts_final.txt",
    output:
        counts="results/metacells/{sample}.genes_cells_filtered_seacells.h5",
        barcodes="results/metacells/{sample}.genes_cells_filtered_seacells_assignment.tsv", 
        #dont care about soft mappings for now
    params:
        prefix="{sample}.genes_cells_filtered",
        outdir="results/metacells/",
        custom_script="workflow/scripts/run_seacells.py",
        various_params=config["tools"]["metacells"]["seacells"]["params"]
    container:
        config["tools"]["metacells"]["seacells"]["container"]
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/metacells/{sample}.genes_cells_filtered_seacells.log",
    benchmark:
        "logs/benchmark/metacells/{sample}.genes_cells_filtered_seacells.benchmark"
    shell:
        """
        python {params.custom_script} \
            -i {input.counts} \
            -o {params.outdir} \
            -p {params.prefix} \
            -c {input.celltypes} \
            {params.various_params} \
            &> {log}
        """

# for other rules inbetween, wildcard {sample} becomes {sample}.genes_cells_filtered_seacells    
 
use rule sctransform_preprocessing as sctransform_preprocessing_filtered_seacells with:
    input:
        hdf5_file="results/metacells/{sample}.genes_cells_filtered_seacells.h5",
    output:
        outfile="results/counts_corrected/{sample}_seacells.corrected.RDS",
        highly_variable="results/counts_corrected/{sample}_seacells.corrected.variable_genes.h5",
    params:
        sample="{sample}_seacells",
        number_genes=config["tools"]["sctransform_preprocessing"]["number_genes_metacells"],
        min_var=config["tools"]["sctransform_preprocessing"]["min_var_metacells"],
        n_nn=config["tools"]["sctransform_preprocessing"]["n_nn_metacells"],
        outDir="results/counts_corrected/",
        custom_script=workflow.source_path("../scripts/sctransform_preprocessing.R"),


rule evaluate_seacells:
    input:
        celltypes="results/celltyping/{sample}.cts_final.txt",
        metacelltypes="results/celltyping/{sample}_seacells.cts_final.txt",
        assignment="results/metacells/{sample}.genes_cells_filtered_seacells_assignment.tsv",
    output:
        table="results/metacells/{sample}_metacell_celltype_counts.tsv",
        plot="results/metacells/{sample}_seacells_celltype_dens.png",
        report="workflow/report/rules/metacell/{sample}_celltyping_summary.rst"
    params:
        prefix="{sample}",
        outdir="results/metacells/",
        custom_script="workflow/scripts/cmp_seacells_celltypes.py",
        
    container:
        config["tools"]["metacells"]["seacells"]["container"]
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["low"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    log:
        "logs/metacells/{sample}_metacell_evaluation.log",
    benchmark:
        "logs/benchmark/metacells/{sample}_metacell_evaluation.benchmark",
    shell:
        """
        python {params.custom_script} \
            -i {input.celltypes} \
            -m {input.metacelltypes} \
            -a {input.assignment} \
            -o {params.outdir} \
            -p {params.prefix} \
            -r {output.report} \
            &> {log}
        """
