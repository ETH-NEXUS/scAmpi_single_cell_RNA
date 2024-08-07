try:
    from snakemake import logger
except ImportError:
    # change in sm8
    from snakemake.logging import logger# Retrieve the fastqs directory name (ie. uses cellranger sample name) corresponding to a given sample

def get_fastq_dir(wildcards):
    "return fastq directory of one sample"
    sample = wildcards.sample
    fastqs_dir = config["inputOutput"]["input_fastqs"]
    sample_fastq_dir = fastqs_dir + sample
    return sample_fastq_dir


# make symlinks to fastq files that work with cellranger
rule fastq_symlinks:
    input:
        lambda w: f"{config['inputOutput']['input_fastqs']}/{get_input_fastq(w)}"
    output:
        # attempt to check whether the filename is cellranger conform 
        target="results/input_fastq/{link_filename,[a-zA-Z0-9_-]+\\.fastq\\.gz}"
    shell:
        "ln -s {input} {output.target}"

# cellranger call to process the raw samples
rule cellranger_count:
    input:
        #fastqs_dir=config["inputOutput"]["input_fastqs"],
        fastq_files=lambda w: [f'results/input_fastq/{fn}' for fn in get_symlink_names(w)],
        reference=config["resources"]["reference_transcriptome"],
    output:
        success="results/cellranger_run/{sample}_success_cellranger.txt",
    params:
        fastqs_dir="../input_fastq/", # needs to be relative to cr_out
        cr_out="results/cellranger_run/",
        local_cores=config["tools"]["cellranger_count"]["local_cores"],
        variousParams=config["tools"]["cellranger_count"]["variousParams"],
    resources:
        mem_mb=config["tools"]["cellranger_count"]["mem_mb"],
        runtime=config["tools"]["cellranger_count"]["runtime"],
    threads: config["tools"]["cellranger_count"]["local_cores"]
    log:
        "logs/cellranger_count/{sample}.log",
    benchmark:
        "logs/benchmark/cellranger_run/{sample}.benchmark"
    # NOTE: cellranger count function cannot specify the output directory, the output is the path you call it from.
    # Therefore, a subshell is used here.
    shell:
        "(cd {params.cr_out}; "
        "{config[tools][cellranger_count][call]} count "
        "--id={wildcards.sample} "
        "--sample={wildcards.sample} "
        "--transcriptome={input.reference} "
        "--localcores={params.local_cores} "
        "--fastqs={params.fastqs_dir} "
        "--nosecondary "
        "{params.variousParams} "
        " 2>&1 | tee ../../{log}); "
        "date > {output.success} "


# Run cellranger v8. Some new parameters are required (e.g. --create-bam)
rule cellranger_count_8:
    input:
        fastqs_dir=get_fastq_dir,
        reference=config["resources"]["reference_transcriptome"],
    output:
        # the cellranger output cannot be specified directly because this would trigger Snakemake
        # to create the directory before cellranger starts.
        # Cellranger needs to create the directory itself (otherwise gives pipestance error and aborts).
        success="results/cellranger_run/{sample}_success_cellranger.txt",
    params:
        cr_out="results/cellranger_run/{sample}/",
        local_cores=config["tools"]["cellranger_count"]["local_cores"],
        metrics_summary="results/cellranger_run/{sample}.metrics_summary.csv",
        web_summary="results/cellranger_run/{sample}.web_summary.html",
        create_bam=config["tools"]["cellranger_count"]["create_bam"],
        # NOTE: no dots are allowed in sample names!
        variousParams=config["tools"]["cellranger_count"]["variousParams"],
    resources:
        mem_mb=config["tools"]["cellranger_count"]["mem_mb"],
        runtime=config["tools"]["cellranger_count"]["runtime"],
    threads: config["tools"]["cellranger_count"]["local_cores"]
    log:
        "logs/cellranger_count/{sample}.log",
    benchmark:
        "logs/benchmark/cellranger_run/{sample}.benchmark"
    shell:
        "{config[tools][cellranger_count][call]} count "
        "--id={wildcards.sample} "
        "--transcriptome={input.reference} "
        "--localcores={params.local_cores} "
        "--fastqs={input.fastqs_dir} "
        "--nosecondary "
        "--create-bam={params.create_bam} "
        "--output-dir={params.cr_out} "
        "{params.variousParams} "
        " 2>&1 | tee {log} ; "
        "date > {output}"


# get sample ID to output files
rule gunzip_and_link_cellranger:
    input:
        success="results/cellranger_run/{sample}_success_cellranger.txt",
    output:
        features_file="results/cellranger_run/{sample}.features.tsv",
        matrix_file="results/cellranger_run/{sample}.matrix.mtx",
        barcodes_file="results/cellranger_run/{sample}.barcodes.tsv",
    params:
        cr_out="results/cellranger_run/{sample}/",
        metrics_summary="results/cellranger_run/{sample}.metrics_summary.csv",
        web_summary="results/cellranger_run/{sample}.web_summary.html",
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["low"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    benchmark:
        "logs/benchmark/gunzip_and_link_cellranger/{sample}.benchmark"
    log:
        "logs/rules/gunzip_and_link_cellranger/{sample}.log",
    # unzip and symlink raw cellranger features file
    shell:
        "gunzip --keep {params.cr_out}/outs/filtered_feature_bc_matrix/features.tsv.gz 2>> {log} ; "
        "gunzip --keep {params.cr_out}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz 2>> {log} ; "
        "gunzip --keep {params.cr_out}/outs/filtered_feature_bc_matrix/matrix.mtx.gz 2>> {log} ; "
        'ln -sr "{params.cr_out}/outs/filtered_feature_bc_matrix/features.tsv" "{output.features_file}" 2>> {log} ; '
        'ln -sr "{params.cr_out}/outs/filtered_feature_bc_matrix/matrix.mtx" "{output.matrix_file}" 2>> {log} ; '
        'ln -sr "{params.cr_out}/outs/filtered_feature_bc_matrix/barcodes.tsv" "{output.barcodes_file}" 2>> {log} ; '
        'ln -sr "{params.cr_out}/outs/web_summary.html" "{params.web_summary}" 2>> {log} ; '
        'ln -sr "{params.cr_out}/outs/metrics_summary.csv" "{params.metrics_summary}" 2>> {log} '

        


# create hdf5 from count matrix, genes, and cell barcodes file
rule create_hdf5:
    input:
        genes_file="results/cellranger_run/{sample}.features.tsv",
        matrix_file="results/cellranger_run/{sample}.matrix.mtx",
        barcodes_file="results/cellranger_run/{sample}.barcodes.tsv",
    output:
        outfile="results/counts_raw/{sample}.h5",
    params:
        custom_script=workflow.source_path("../scripts/create_hdf5.py"),
    conda:
        "../envs/create_hdf5.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/create_hdf5/{sample}.log",
    benchmark:
        "logs/benchmark/create_hdf5/{sample}.benchmark"
    shell:
        "python {params.custom_script} "
        "-g {input.genes_file} "
        "-m {input.matrix_file} "
        "-b {input.barcodes_file} "
        "-o {output.outfile} "
        " &> {log} "


# identify doublets with scDblFinder
rule identify_doublets:
    input:
        infile="results/counts_raw/{sample}.h5",
    output:
        outfile="results/identify_doublets/{sample}.doublet_barcodes.txt",
    params:
        sample="{sample}",
        outdir="results/identify_doublets/",
        custom_script=workflow.source_path("../scripts/identify_doublets.R"),
    conda:
        "../envs/identify_doublets.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/identify_doublets/{sample}.log",
    benchmark:
        "logs/benchmark/identify_doublets/{sample}.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--hdf5File {input.infile} "
        "--sample {params.sample} "
        "--outdir {params.outdir} "
        "&> {log} "


# filter out whole genes and cells.
# Filter out cells if the fraction of reads mapped to MT genes is too high or the NODG is too low.
# Filter out genes if they are not protein coding, if they are mitochondrial genes, or if they encode for ribosomal proteins.
rule filter_genes_and_cells:
    input:
        infile="results/counts_raw/{sample}.h5",
        doublets="results/identify_doublets/{sample}.doublet_barcodes.txt",
    output:
        outfile="results/counts_filtered/{sample}.genes_cells_filtered.h5",
    params:
        nmads_fractionMT=config["tools"]["filter_genes_and_cells"]["nmads_fractionMT"],
        nmads_NODG=config["tools"]["filter_genes_and_cells"]["nmads_NODG"],
        threshold_fractionMT=config["tools"]["filter_genes_and_cells"][
            "threshold_fractionMT"
        ],
        threshold_NODG=config["tools"]["filter_genes_and_cells"]["threshold_NODG"],
        remove_doublets=config["tools"]["filter_genes_and_cells"]["remove_doublets"],
        minNumberCells=config["tools"]["filter_genes_and_cells"]["minNumberCells"],
        protein_coding_only=config["tools"]["filter_genes_and_cells"][
            "protein_coding_only"
        ],
        outDir="results/counts_filtered/",
        genomeVersion=config["tools"]["filter_genes_and_cells"]["genomeVersion"],
        sample="{sample}",
        custom_script=workflow.source_path("../scripts/filter_genes_and_cells.R"),
    conda:
        "../envs/filter_genes_and_cells.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/filter_genes_and_cells/{sample}.log",
    benchmark:
        "logs/benchmark/filter_genes_and_cells/{sample}.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--hdf5File {input.infile} "
        "--nmads_NODG {params.nmads_NODG} "
        "--nmads_fractionMT {params.nmads_fractionMT} "
        "--threshold_NODG {params.threshold_NODG} "
        "--threshold_fractionMT {params.threshold_fractionMT} "
        "--genomeVersion {params.genomeVersion} "
        "--doublet_barcodes {input.doublets} "
        "--remove_doublets {params.remove_doublets} "
        "--minNumberCells {params.minNumberCells} "
        "--protein_coding_only {params.protein_coding_only} "
        "--sample {params.sample} "
        "--outDir {params.outDir} "
        "&> {log} "

def get_sctransform_param(wildcards, param):
    if wildcards.sample.split('_')[-1]=='seacells':
        return config["tools"]["sctransform_preprocessing"][param+'_metacells']
    else:
        return config["tools"]["sctransform_preprocessing"][param]

# perform normalisation, cell cycle correction and other preprocessing using sctransform
rule sctransform_preprocessing:
    input:
        hdf5_file="results/counts_filtered/{sample}.genes_cells_filtered.h5",
    output:
        outfile="results/counts_corrected/{sample, '^(?!.*_seacells$).+'}.corrected.RDS",
        highly_variable="results/counts_corrected/{sample, '^(?!.*_seacells$).+'}.corrected.variable_genes.h5",
    params:
        sample="{sample}",
        number_genes=config["tools"]["sctransform_preprocessing"]["number_genes"],
        min_var=config["tools"]["sctransform_preprocessing"]["min_var"],
        n_nn=config["tools"]["sctransform_preprocessing"]["n_nn"],
        outDir="results/counts_corrected/",
        #custom_script=workflow.source_path("../scripts/sctransform_preprocessing.R"),
        custom_script="workflow/scripts/sctransform_preprocessing.R",
        smooth_pc="100", # default value
        max_count="0", # no effect
        min_cells_per_gene="0" # do not filter
    conda:
        "../envs/sctransform_preprocessing.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["medium"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/sctransform_preprocessing/{sample}.log",
    benchmark:
        "logs/benchmark/sctransform_preprocessing/{sample}.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--inHDF5 {input.hdf5_file} "
        "--sample {params.sample} "
        "--number_genes {params.number_genes} "
        "--min_var {params.min_var} "
        "--n_nn {params.n_nn} "
        "--max_pc_smooth {params.smooth_pc} "
        "--max_count {params.max_count} "
        "--min_cells_per_gene {params.min_cells_per_gene} "
        "--outdir {params.outDir} "
        "&> {log} "


# perform clustering with phenograph
rule phenograph:
    input:
        infile="results/counts_corrected/{sample}.corrected.variable_genes.h5",
    output:
        outfile="results/clustering/{sample}.clusters_phenograph.csv",
        distance_matrix="results/clustering/{sample}.distance_matrix.tsv",
        modularity_score="results/clustering/{sample}.modularity_score.txt",
    params:
        n_neighbours=config["tools"]["clustering"]["phenograph"]["n_neighbours"],
        min_cluster_size=config["tools"]["clustering"]["phenograph"]["min_cluster_size"],
        log_normalize=config["tools"]["clustering"]["phenograph"]["log_normalize"],
        custom_script=workflow.source_path("../scripts/apply_phenograph.py"),
    conda:
        "../envs/phenograph.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/phenograph/{sample}.log",
    benchmark:
        "logs/benchmark/phenograph/{sample}.benchmark"
    shell:
        "python {params.custom_script} "
        "--input_file {input.infile} "
        "--output_file {output.outfile} "
        "--distance_matrix {output.distance_matrix} "
        "--modularity_score {output.modularity_score} "
        "--n_neighbours {params.n_neighbours} "
        "--min_size {params.min_cluster_size} "
        "--log_normalize {params.log_normalize} "
        "--n_threads {threads} "
        "&> {log} "


# prepare sce object in RDS file for cell type classification
rule prepare_celltyping:
    input:
        RDS_file="results/counts_corrected/{sample}.corrected.RDS",
        cluster="results/clustering/{sample}.clusters_phenograph.csv",
        distanceMatrix="results/clustering/{sample}.distance_matrix.tsv",
        modularity_score="results/clustering/{sample}.modularity_score.txt",
    output:
        outfile="results/prep_celltyping/{sample}.prep_celltyping.RDS",
    params:
        outputDirec="results/prep_celltyping/",
        sampleName="{sample}",
        custom_script=workflow.source_path("../scripts/prepare_celltyping.R"),
    conda:
        "../envs/prepare_celltyping.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/prepare_celltyping/{sample}.log",
    benchmark:
        "logs/benchmark/prepare_celltyping/{sample}.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--in_sce {input.RDS_file} "
        "--phenograph_cluster {input.cluster} "
        "--outputDirec {params.outputDirec} "
        "--sampleName {params.sampleName} "
        "--distanceMatrix {input.distanceMatrix} "
        "--modularity_score {input.modularity_score} "
        "&> {log} "


# perform cell type classification
rule celltyping:
    input:
        infile="results/prep_celltyping/{sample}.prep_celltyping.RDS",
    output:
        outfile="results/celltyping/{sample}.celltyping.phenograph_celltype_association.txt",
        out_sce="results/celltyping/{sample}.celltyping.RDS",
        out_barcodes="results/celltyping/{sample}.cts_final.txt",
    params:
        min_genes=config["tools"]["celltyping"]["min_genes"],
        celltype_lists=config["resources"]["celltype_lists"],
        celltype_config=config["resources"]["celltype_config"],
        outputDirec="results/celltyping/",
        sampleName="{sample}",
        custom_script=workflow.source_path("../scripts/celltyping.R"),
    conda:
        "../envs/celltyping.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["medium"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/celltyping/{sample}.log",
    benchmark:
        "logs/benchmark/celltyping/{sample}.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--SCE {input.infile} "
        "--celltype_lists {params.celltype_lists} "
        "--celltype_config {params.celltype_config} "
        "--sampleName {params.sampleName} "
        "--min_genes {params.min_genes} "
        "--outputDirec {params.outputDirec} "
        "&> {log} "

def get_atypical_params(wildcard, key):
    what='cells'  
    if '_seacells' in wildcard['sample'] or '_metacells2' in wildcard['sample']:
        what = 'metacells'
    return config["tools"]["remove_atypical_cells"][what][key]

# filter out atypical cells from sce object
rule remove_atypical_cells:
    input:
        infile="results/celltyping/{sample}.celltyping.RDS",
        cluster_table="results/celltyping/{sample}.celltyping.phenograph_celltype_association.txt",
    output:
        out_sce="results/atypical_removed/{sample}.atypical_removed.RDS",
        out_table="results/atypical_removed/{sample}.atypical_removed.phenograph_celltype_association.txt",
        out_h5file="results/atypical_removed/{sample}.atypical_removed.h5",
    params:
        celltype_config=config["resources"]["celltype_config"],
        outputDirec="results/atypical_removed/",
        sample_name="{sample}",
        threshold_filter=lambda w:get_atypical_params(w,"threshold_filter"),
        min_threshold=lambda w:get_atypical_params(w,"min_threshold"),
        threshold_type=lambda w:get_atypical_params(w,"threshold_type"),
        custom_script=workflow.source_path("../scripts/remove_atypical_cells.R"),
    conda:
        "../envs/remove_atypical_cells.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/remove_atypical_cells/{sample}.log",
    benchmark:
        "logs/benchmark/remove_atypical_cells/{sample}.benchmark"
    shell: 
        "Rscript {params.custom_script} "
        "--sce_in {input.infile} "
        "--cluster_table {input.cluster_table} "
        "--celltype_config {params.celltype_config} "
        "--threshold_filter {params.threshold_filter} "
        "--min_threshold {params.min_threshold} "
        "--threshold_type {params.threshold_type} "
        "--outDir {params.outputDirec} "
        "--sample_name {params.sample_name} "
        "&> {log} "


# perform gsva gene set analysis
rule gsva:
    input:
        infile="results/atypical_removed/{sample}.atypical_removed.RDS",
    output:
        outfile="results/gsva/{sample}.gsetscore_hm.png",
    params:
        outputDirec="results/gsva/",
        sampleName="{sample}",
        genesets=config["resources"]["genesets"],
        custom_script=workflow.source_path("../scripts/gsva.R"),
    conda:
        "../envs/gsva.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["medium"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/gsva/{sample}.log",
    benchmark:
        "logs/benchmark/gsva/{sample}.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--SCE {input.infile} "
        "--geneset {params.genesets} "
        "--outputDirec {params.outputDirec} "
        "--sampleName {params.sampleName} "
        "&> {log} "


# generate plots about sample composition and gene expression
rule plotting:
    input:
        infile="results/atypical_removed/{sample}.atypical_removed.RDS",
    output:
        outfile="results/plotting/{sample}.celltype_barplot.png",
    params:
        outputDirec="results/plotting/",
        sampleName="{sample}",
        genes_of_interest=config["resources"]["priority_genes"],
        colour_config=config["resources"]["colour_config"],
        use_alias=config["tools"]["plotting"]["use_alias"],
        custom_script=workflow.source_path("../scripts//plotting.R"),
    conda:
        "../envs/plotting.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["medium"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/plotting/{sample}.log",
    benchmark:
        "logs/benchmark/plotting/{sample}.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--sce_in {input.infile} "
        "--genelist {params.genes_of_interest} "
        "--outDir {params.outputDirec} "
        "--sampleName {params.sampleName} "
        "--colour_config {params.colour_config} "
        "--toggle_label {params.use_alias} "
        "&> {log} "


# give out gene expression values per cluster
rule gene_exp:
    input:
        sce_in="results/atypical_removed/{sample}.atypical_removed.RDS",
    output:
        out="results/gene_exp/{sample, '^(?!.*_seacells$).+'}.gene_expression_per_cluster.tsv",
    params:
        sampleName="{sample}",
        outpath="results/gene_exp/",
        threshold_sample=config["tools"]["gene_exp"]["threshold_sample"],
        type_sample=config["tools"]["gene_exp"]["type_sample"],
        priority_genes=config["resources"]["priority_genes"],
        custom_script=workflow.source_path("../scripts/gene_exp.R"),
    conda:
        "../envs/gene_exp.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/gene_exp/{sample}.log",
    benchmark:
        "logs/benchmark/gene_exp/{sample}.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--sce_in {input.sce_in} "
        "--priority_genes {params.priority_genes} "
        "--filtering_threshold_sample {params.threshold_sample} "
        "--filter_type_sample {params.type_sample} "
        "--outDir {params.outpath} "
        "--sample_name {params.sampleName} "
        "&> {log} "


# This rule generates general quality control plots to the raw hdf5 expression files
rule generate_qc_plots_raw:
    input:
        infile="results/counts_raw/{sample}.h5",
    output:
        out="results/qc_plots/raw/{sample}.raw.histogram_library_sizes.png",
    params:
        custom_script=workflow.source_path("../scripts/generate_QC_plots.R"),
        outdir="results/qc_plots/raw/",
        sample_status="raw",
    conda:
        "../envs/generate_qc_plots.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/generate_qc_plots/{sample}.raw.log",
    benchmark:
        "logs/benchmark/generate_qc_plots/{sample}.raw.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--hdf5File {input.infile} "
        "--sample_name {wildcards.sample} "
        "--sample_status {params.sample_status} "
        "--outdir {params.outdir} "
        "&> {log} "


# This rule generates general quality control plots to the raw hdf5 expression files
rule generate_qc_plots_filtered:
    input:
        infile="results/counts_filtered/{sample}.genes_cells_filtered.h5",
    output:
        out="results/qc_plots/filtered/{sample}.genes_cells_filtered.histogram_library_sizes.png",
    params:
        custom_script=workflow.source_path("../scripts/generate_QC_plots.R"),
        outdir="results/qc_plots/filtered/",
        sample_status="genes_cells_filtered",
    conda:
        "../envs/generate_qc_plots.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["medium"],
        runtime=config["computingResources"]["runtime"]["low"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/generate_qc_plots/{sample}.genes_cells_filtered.log",
    benchmark:
        "logs/benchmark/generate_qc_plots/{sample}.genes_cells_filtered.benchmark"
    shell:
        "Rscript {params.custom_script} "
        "--hdf5File {input.infile} "
        "--sample_name {wildcards.sample} "
        "--sample_status {params.sample_status} "
        "--outdir {params.outdir} "
        "&> {log} "


# perform the differential expression analysis using a Wilcoxon test
checkpoint diff_exp_analysis:
    input:
        sce_in="results/atypical_removed/{sample}.atypical_removed.RDS",
        cell_types="results/atypical_removed/{sample}.atypical_removed.phenograph_celltype_association.txt",
    output:
        output=directory("results/diff_exp_analysis/{sample, '^(?!.*_seacells$).+'}"),
    params:
        sampleName="{sample}",
        malignant=config["inputOutput"]["malignant_cell_type"],
        threshold_comparison=config["tools"]["diff_exp_analysis"][
            "threshold_comparison"
        ],
        fdr_cut=config["tools"]["diff_exp_analysis"]["fdr_cut"],
        fc_cut=config["tools"]["diff_exp_analysis"]["fc_cut"],
        mindiff2second=config["tools"]["diff_exp_analysis"]["mindiff2second"],
        minNumberNonMalignant=config["tools"]["diff_exp_analysis"][
            "minNumberNonMalignant"
        ],
        outpath="results/diff_exp_analysis/{sample}/",
        custom_script=workflow.source_path("../scripts/diff_exp_analysis.R"),
    conda:
        "../envs/diff_exp_analysis.yaml"
    resources:
        mem_mb=config["computingResources"]["mem_mb"]["high"],
        runtime=config["computingResources"]["runtime"]["high"],
    threads: config["computingResources"]["threads"]["medium"]
    log:
        "logs/diff_exp_analysis/{sample}.log",
    benchmark:
        "logs/benchmark/diff_exp_analysis/{sample}.benchmark"
    shell:
        "mkdir {params.outpath} ; "
        "Rscript {params.custom_script} "
        "--sample_data {input.sce_in} "
        "--sampleName {params.sampleName} "
        "--cluster_table {input.cell_types} "
        "--malignant_tag {params.malignant} "
        "--fdr_cut {params.fdr_cut} "
        "--fc_cut {params.fc_cut} "
        "--mindiff2second {params.mindiff2second} "
        "--threshold_comparison {params.threshold_comparison} "
        "--minNumberNonMalignant {params.minNumberNonMalignant} "
        "--outdir {params.outpath} "
        "&> {log} "
