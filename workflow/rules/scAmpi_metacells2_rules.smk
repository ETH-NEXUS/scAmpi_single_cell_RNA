# this rule gets a raw h5 file and runs metacells2 to aggreate cells to metacells
# it returns an h5 file with aggregated counts
rule metacells2:
    input:
        counts="results/counts_filtered/{sample}.genes_cells_filtered.h5",
        celltypes="results/celltyping/{sample}.cts_final.txt",
    output:
        counts="results/metacells/{sample}.genes_cells_filtered_metacells2.h5",
        barcodes="results/metacells/{sample}.genes_cells_filtered_metacells2_assignment.tsv",
    params:
        prefix="{sample}.genes_cells_filtered",
        outdir="results/metacells/",
        custom_script="workflow/scripts/metacell_run_metacells2.py",
        various_params=config["tools"]["metacells"]["metacells2"]["params"],
    container:
        config["tools"]["metacells"]["metacells2"]["container"]
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/metacells/{sample}.genes_cells_filtered_metacells2.log",
    benchmark:
        "logs/benchmark/metacells/{sample}.genes_cells_filtered_metacells2.benchmark"
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


use rule sctransform_preprocessing as sctransform_preprocessing_filtered_metacells2 with:
    input:
        hdf5_file="results/metacells/{sample}.genes_cells_filtered_metacells2.h5",
    output:
        outfile="results/counts_corrected/{sample}_metacells2.corrected.RDS",
        highly_variable="results/counts_corrected/{sample}_metacells2.corrected.variable_genes.h5",
    params:
        sample="{sample}_metacells2",
        number_genes=config["tools"]["sctransform_preprocessing"][
            "number_genes_metacells"
        ],
        min_var=config["tools"]["sctransform_preprocessing"]["min_var_metacells"],
        n_nn=config["tools"]["sctransform_preprocessing"]["n_nn_metacells"],
        outDir="results/counts_corrected/",
        #custom_script=workflow.source_path("../scripts/sctransform_preprocessing.R"),
        custom_script="workflow/scripts/sctransform_preprocessing.R",
        smooth_pc="20",
        max_count="0",  # no 
        min_cells_per_gene="10",
    log:
        "logs/sctransform_preprocessing/{sample}_metacells2.log",


rule evaluate_metacells2:
    input:
        celltypes="results/celltyping/{sample}.cts_final.txt",
        metacelltypes="results/celltyping/{sample}_metacells2.cts_final.txt",
        assignment="results/metacells/{sample}.genes_cells_filtered_metacells2_assignment.tsv",
    output:
        table="results/metacells/{sample}_metacells2_celltype_counts.tsv",
        plot="results/metacells/{sample}_metacells2_celltype_dens.png",
        report="workflow/report/rules/metacells2/{sample}_celltyping_summary.rst",
    params:
        prefix="{sample}_metacells2",
        outdir="results/metacells/",
        custom_script="workflow/scripts/metacell_cmp_celltypes.py",
        ct_config=config["resources"]["celltype_config"],
    container:
        config["tools"]["metacells"]["seacells"]["container"]
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["low"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    log:
        "logs/metacells/{sample}_metacell_evaluation_metacells2.log",
    benchmark:
        "logs/benchmark/metacells/{sample}_metacell_evaluation_metacells2.benchmark"
    shell:
        """
        python {params.custom_script} \
            -i {input.celltypes} \
            -m {input.metacelltypes} \
            -a {input.assignment} \
            -c {params.ct_config} \
            -o {params.outdir} \
            -p {params.prefix} \
            -r {output.report} \
            &> {log}
        """
