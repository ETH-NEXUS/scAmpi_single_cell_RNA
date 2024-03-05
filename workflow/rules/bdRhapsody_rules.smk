# use the output from BD Rhapsody™ Sequence Analysis Pipeline as an input for scampi
# ToDo check if it makes sense to add the bd tool directly as part of the pipeline
# at the moment: link relevant output of bdpipeline into one folder & rename files to follow the schema
# {sample}_features.tsv.gz, {sample}_matrix.mtx.gz, {sample}_barcodes.tsv.gz


rule unzip_bdrhapsody_matrix:
    input:
        matrix_dir=config["inputOutput"]["input_bdr_matrix"],
    output:
        features_file="results/bdr_matrix/{sample}.features.tsv",
        matrix_file="results/bdr_matrix/{sample}.matrix.mtx",
        barcodes_file="results/bdr_matrix/{sample}.barcodes.tsv",
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
        "gunzip -c {input.matrix_dir}/{params.mySample}_features.tsv.gz > {params.br_out}/{params.mySample}.features.tsv; "
        "gunzip -c {input.matrix_dir}/{params.mySample}_barcodes.tsv.gz > {params.br_out}/{params.mySample}.barcodes.tsv; "
        "gunzip -c {input.matrix_dir}/{params.mySample}_matrix.mtx.gz > {params.br_out}/{params.mySample}.matrix.mtx; "
