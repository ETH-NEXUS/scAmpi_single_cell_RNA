# use the output from BD Rhapsody™ Sequence Analysis Pipeline as an input for scampi
# ToDo check if it makes sense to add the bd tool directly as part of the pipeline
# at the moment: link relevant output of bdpipeline into one folder & rename files to follow the schema
# {sample}_features.original.tsv.gz, {sample}_matrix.original.mtx, {sample}_barcodes.tsv
# in the rule adapt_gene_id the first column of the feature file that is containing gene symbols is replaced by ensembl ids 
# duplicated entries and na entries are then removed 


rule unzip_bdrhapsody_matrix:
    input:
        matrix_dir=config["inputOutput"]["input_bdr_matrix"],
    output:
        features_file="results/bdr_matrix/{sample}.features.original.tsv",
        matrix_file="results/bdr_matrix/{sample}.matrix.original.mtx",
        barcodes_file="results/bdr_matrix/{sample}.barcodes.tsv"
    params:
        br_out="results/bdr_matrix/",
        mySample="{sample}",
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["low"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    log:
        "logs/unzip_bdrhapsody_matrix/{sample}.log",
    benchmark:
        "logs/benchmark/unzip_bdrhapsody_matrix/{sample}.benchmark"
    shell:
        "gunzip -c {input.matrix_dir}/{params.mySample}_features.tsv.gz > {params.br_out}/{params.mySample}.features.original.tsv; "
        "gunzip -c {input.matrix_dir}/{params.mySample}_barcodes.tsv.gz > {params.br_out}/{params.mySample}.barcodes.tsv; "
        "gunzip -c {input.matrix_dir}/{params.mySample}_matrix.mtx.gz > {params.br_out}/{params.mySample}.matrix.original.mtx; "


rule adapt_gene_id:
    input:
        features_file="results/bdr_matrix/{sample}.features.original.tsv",
        matrix_file="results/bdr_matrix/{sample}.matrix.original.mtx"
    output:
        features_file="results/bdr_matrix/{sample}.features.tsv",
        matrix_file="results/bdr_matrix/{sample}.matrix.mtx",
        unmatched="results/bdr_matrix/{sample}.unmatched_genes.txt"
    params:
        singularity=config["tools"]["singularity"],
        call=config["tools"]["adapt_gene_id"]["call"]
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["medium"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/adapt_gene_id/{sample}.log",
    benchmark:
        "logs/benchmark/adapt_gene_id/{sample}.benchmark"
    shell:
        "{params.singularity} {params.call} Rscript workflow/scripts/adapt_gene_id.R -f {input.features_file} -m {input.matrix_file}"
        "&> {log} "
