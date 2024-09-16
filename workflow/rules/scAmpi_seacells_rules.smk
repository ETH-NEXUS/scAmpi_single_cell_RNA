# this rule gets a raw h5 file and runs seacells to aggreate cells to metacells
# it returns an h5 file with aggregated counts
rule seacells:
    input:
        counts="results/counts_filtered/{sample}.genes_cells_filtered.h5",
        celltypes="results/celltyping/{sample}.cts_final.txt",
    output:
        counts="results/metacells/{sample}.genes_cells_filtered_seacells.h5",
        barcodes="results/metacells/{sample}.genes_cells_filtered_seacells_assignment.tsv",
    params:
        prefix="{sample}.genes_cells_filtered",
        outdir="results/metacells/",
        custom_script="workflow/scripts/metacell_run_seacells.py",
        various_params=config["tools"]["metacells"]["seacells"]["params"],
    container:
        "docker://mlienhard/seacells"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["medium"],
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


use rule sctransform_preprocessing as sctransform_preprocessing_filtered_seacells with:
    input:
        hdf5_file="results/metacells/{sample}.genes_cells_filtered_seacells.h5",
    output:
        outfile="results/counts_corrected/{sample}_seacells.corrected.RDS",
        highly_variable="results/counts_corrected/{sample}_seacells.corrected.variable_genes.h5",
    params:
        sample="{sample}_seacells",
        number_genes=config["tools"]["sctransform_preprocessing"][
            "number_genes_metacells"
        ],
        min_var=config["tools"]["sctransform_preprocessing"]["min_var_metacells"],
        n_nn=config["tools"]["sctransform_preprocessing"]["n_nn_metacells"],
        outdir="results/counts_corrected/",
        custom_script=workflow.source_path("../scripts/sctransform_preprocessing.R"),
        smooth_pc="20",
        patch="--patch_vst workflow/scripts/vst_check.R",  # leave empty to not apply patch
    log:
        "logs/sctransform_preprocessing/{sample}_seacells.log",
    benchmark:
        "logs/benchmark/sctransform_preprocessing/{sample}.benchmark"


rule evaluate_seacells:
    input:
        celltypes="results/celltyping/{sample}.cts_final.txt",
        metacelltypes="results/celltyping/{sample}_seacells.cts_final.txt",
        assignment="results/metacells/{sample}.genes_cells_filtered_seacells_assignment.tsv",
    output:
        table="results/metacells/{sample}_seacells_celltype_counts.tsv",
        plot="results/metacells/{sample}_seacells_celltype_hist.png",
        report="workflow/report/rules/seacells/{sample}_celltyping_summary.rst",
    params:
        prefix="{sample}_seacells",
        outdir=lambda w, output: dirname(output.table),
        custom_script="workflow/scripts/metacell_cmp_celltypes.py",
        ct_config=config["resources"]["celltype_config"],
    container:
        "docker://mlienhard/seacells"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["low"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    log:
        "logs/metacells/{sample}_metacell_evaluation_seacells.log",
    benchmark:
        "logs/benchmark/metacells/{sample}_metacell_evaluation_seacells.benchmark"
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


rule metacell_stats:
    input:
        celltypes=expand(
            "results/metacells/{sample}_seacells_celltype_counts.tsv", 
            sample=sample_ids),
    output:
        purity="results/metacells/seacells_celltype_purity.png",
        subpurity="results/metacells/seacells_subcelltype_purity.png",
        ncells="results/metacells/seacells_cells_per_metacell.png",
    params:
        custom_script="workflow/scripts/metacell_stats.py",
        outprefix=lambda w, output: join(dirname(output.purity), "seacells_"),
    container:
        "docker://mlienhard/seacells"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["low"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["low"],
    log:
        "logs/metacells/metacell_stats.log",
    benchmark:
        "logs/benchmark/metacells/metacell_stats.benchmark",
    shell:
        """
        python {params.custom_script} \
            {params.outprefix} \
            {input} \
            &> {log}
        """


