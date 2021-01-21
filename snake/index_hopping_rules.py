# Index hopping related rules for sc GEX data
# Adapted from singlecell pipeline
# Lourdes Rosano, March 2020

if not 'CELLRANGER_NOVA_IN' in globals():
    CELLRANGER_NOVA_IN = INPUTDIR
if not 'CELLRANGER_NOVA_OUT' in globals():
    CELLRANGER_NOVA_OUT = OUTDIR + 'cellranger_run_index_hopping/'

# Retrieve the fastqs directory name (ie. uses cellranger sample name) corresponding to a given sample
def getSampleFastqDir(wildcards):
    crsample = getCellRangerName(wildcards)
    fastq_dir = CELLRANGER_NOVA_IN + str(crsample) + '/'
    return fastq_dir

# Only for index hopping removal framework, cellranger call to process the raw GEX samples
rule cellranger_count_novaSeq:
    input:
        # Subdirectories containing each sample's fastqs use the cellranger sample name
        fastqs_dir = getSampleFastqDir,
        reference = config['resources']['reference_transcriptome']
    output:
        molecule_file = CELLRANGER_NOVA_OUT + '{sample}.molecule_info.h5',
        gene_mapper_file = CELLRANGER_NOVA_OUT + '{sample}.features.raw.with_gene_names.tsv',
        web_file = CELLRANGER_NOVA_OUT + '{sample}.web_summary.html',
        metrics_file = CELLRANGER_NOVA_OUT + '{sample}.metrics_summary.csv'
    params:
        cr_out = CELLRANGER_NOVA_OUT,
        cellranger_sampleName = getCellRangerName,
        local_cores = config['tools']['cellranger_count_novaSeq']['local_cores'],
        lsfoutfile = CELLRANGER_NOVA_OUT + '{sample}.cellranger_count_novaSeq.lsfout.log',
        lsferrfile = CELLRANGER_NOVA_OUT + '{sample}.cellranger_count_novaSeq.lsferr.log',
        scratch = config['tools']['cellranger_count_novaSeq']['scratch'],
        mem = config['tools']['cellranger_count_novaSeq']['mem'],
        time = config['tools']['cellranger_count_novaSeq']['time'],
        variousParams = config['tools']['cellranger_count_novaSeq']['variousParams'],
    threads:
        config['tools']['cellranger_count_novaSeq']['threads']
    benchmark:
        CELLRANGER_NOVA_OUT + '{sample}.cellranger_count_novaSeq.benchmark'
    # NOTE: cellranger count function cannot specify the output directory, the output it the path you call it from.
    # Therefore, a subshell is used here.
    # Also, symlink output molecule file in preparation for the index hopping removal step, and unzip and symlink raw cr features file in preparation for rule 'add_gene_names'
    # Also, create symlinks for the web and metrics summary files
    shell:
        '(cd {params.cr_out}; {config[tools][cellranger_count_novaSeq][call]} count --id={params.cellranger_sampleName} --transcriptome={input.reference} --localcores={params.local_cores} --fastqs={input.fastqs_dir} --nosecondary {params.variousParams}) ; gunzip {params.cr_out}{params.cellranger_sampleName}/outs/raw_feature_bc_matrix/features.tsv.gz ; ln -s "{params.cr_out}{params.cellranger_sampleName}/outs/raw_feature_bc_matrix/features.tsv" "{output.gene_mapper_file}" ; ln -s "{params.cr_out}{params.cellranger_sampleName}/outs/molecule_info.h5" "{output.molecule_file}" ; ln -s "{params.cr_out}{params.cellranger_sampleName}/outs/web_summary.html" "{output.web_file}" ; ln -s "{params.cr_out}{params.cellranger_sampleName}/outs/metrics_summary.csv" "{output.metrics_file}"'


# Require all molecule_info.h5 files (ie. from all samples in the run) to be ready before 'remove_index_hopping'
def getInputFiles_indexHopping(wildcards):
    return expand(CELLRANGER_NOVA_OUT + '{sample}.molecule_info.h5', sample = getSampleNames())


if not 'INDEX_HOPPING_OUT' in globals():
    INDEX_HOPPING_OUT = OUTDIR + 'index_hopping_removed/'

# Apply method for index hopping removal to all samples in current run
rule remove_index_hopping:
    input:
        inFiles = getInputFiles_indexHopping
    output:
        rdata_object = INDEX_HOPPING_OUT + 'cleaned_unfiltered.RData',
        successFile = INDEX_HOPPING_OUT + 'complete_index_hopping.txt'
    params:
        ih_out = INDEX_HOPPING_OUT,
        lsfoutfile = INDEX_HOPPING_OUT + 'remove_index_hopping.lsfout.log',
        lsferrfile = INDEX_HOPPING_OUT + 'remove_index_hopping.lsferr.log',
        scratch = config['tools']['remove_index_hopping']['scratch'],
        mem = config['tools']['remove_index_hopping']['mem'],
        time = config['tools']['remove_index_hopping']['time'],
        fdr_thres = config['tools']['remove_index_hopping']['fdr_thres'],
        libsize_thres = config['tools']['remove_index_hopping']['libsize_thres']
    threads:
        config['tools']['remove_index_hopping']['threads']
    benchmark:
        INDEX_HOPPING_OUT + 'remove_index_hopping.benchmark'
    shell:
        '{config[tools][remove_index_hopping][call]} {output.rdata_object} {params.ih_out} {params.fdr_thres} {params.libsize_thres} {input.inFiles} ; date > {output.successFile}'


# Prepare output files from index hopping removal for 'add_gene_names' and 'link_files'
rule prepare_indexHopping_files:
    input:
        successFile = INDEX_HOPPING_OUT + 'complete_index_hopping.txt',
    output:
        features_file = INDEX_HOPPING_OUT + '{sample}/features.index_hopping_removed.tsv',
        matrix_file = INDEX_HOPPING_OUT + '{sample}/matrix.index_hopping_removed.mtx',
        barcodes_file = INDEX_HOPPING_OUT + '{sample}/barcodes.index_hopping_removed.tsv'
    params:
        ih_out = INDEX_HOPPING_OUT,
        lsfoutfile = INDEX_HOPPING_OUT + '{sample}/prepare_indexHopping_files.lsfout.log',
        lsferrfile = INDEX_HOPPING_OUT + '{sample}/prepare_indexHopping_files.lsferr.log',
        scratch = config['tools']['prepare_indexHopping_files']['scratch'],
        mem = config['tools']['prepare_indexHopping_files']['mem'],
        time = config['tools']['prepare_indexHopping_files']['time']
    threads:
        config['tools']['prepare_indexHopping_files']['threads']
    benchmark:
        INDEX_HOPPING_OUT + '{sample}/prepare_indexHopping_files.benchmark'
    # unzip resulting files (features, barcodes, matrix) and change their names for rule 'add_gene_names'
    shell:
        'gunzip {params.ih_out}{wildcards.sample}/features.tsv.gz && mv {params.ih_out}{wildcards.sample}/features.tsv {output.features_file} ; gunzip {params.ih_out}{wildcards.sample}/barcodes.tsv.gz && mv {params.ih_out}{wildcards.sample}/barcodes.tsv {output.barcodes_file} ; gunzip {params.ih_out}{wildcards.sample}/matrix.mtx.gz && mv {params.ih_out}{wildcards.sample}/matrix.mtx {output.matrix_file}'


# Map ensemblIDs (default output format of index hopping) back to gene symbols, as required by pipeline
rule add_gene_names:
    input:
        features_file = INDEX_HOPPING_OUT + '{sample}/features.index_hopping_removed.tsv',
        gene_mapper_file = CELLRANGER_NOVA_OUT + '{sample}.features.raw.with_gene_names.tsv'
    output:
        features_file_with_genes = INDEX_HOPPING_OUT + '{sample}/features.index_hopping_removed.with_gene_names.tsv'
    params:
        lsfoutfile = INDEX_HOPPING_OUT + '{sample}/add_gene_names.lsfout.log',
        lsferrfile = INDEX_HOPPING_OUT + '{sample}/add_gene_names.lsferr.log',
        scratch = config['tools']['add_gene_names']['scratch'],
        mem = config['tools']['add_gene_names']['mem'],
        time = config['tools']['add_gene_names']['time']
    threads:
        config['tools']['add_gene_names']['threads']
    benchmark:
        INDEX_HOPPING_OUT + '{sample}/add_gene_names.benchmark'
    shell:
        '{config[tools][add_gene_names][call]} {input.features_file} {input.gene_mapper_file} {output.features_file_with_genes}'

# Create symlinks from index hopping removal folder to preprocessing (novaSeq) cellranger folder
# Preparation for final rule 'create_symlinks'
rule link_files:
    input:
        # NOTE: input file 'features_file_no_indexHopping' requires extra rule 'add_gene_names' compared to matrix and barcodes
        features_file_no_indexHopping = INDEX_HOPPING_OUT + '{sample}/features.index_hopping_removed.with_gene_names.tsv',
        matrix_file_no_indexHopping = INDEX_HOPPING_OUT + '{sample}/matrix.index_hopping_removed.mtx',
        barcodes_file_no_indexHopping = INDEX_HOPPING_OUT + '{sample}/barcodes.index_hopping_removed.tsv'
    output:
        features_file = CELLRANGER_NOVA_OUT + '{sample}.features.tsv',
        matrix_file = CELLRANGER_NOVA_OUT +'{sample}.matrix.mtx',
        barcodes_file = CELLRANGER_NOVA_OUT + '{sample}.barcodes.tsv'
    params:
        lsfoutfile = CELLRANGER_NOVA_OUT + '{sample}.link_files.lsfout.log',
        lsferrfile = CELLRANGER_NOVA_OUT + '{sample}.link_files.lsferr.log',
        scratch = config['tools']['link_files']['scratch'],
        mem = config['tools']['link_files']['mem'],
        time = config['tools']['link_files']['time']
    threads:
        config['tools']['link_files']['threads']
    benchmark:
        CELLRANGER_NOVA_OUT + '{sample}.link_files.benchmark'
    shell:
        'ln -s "{input.features_file_no_indexHopping}" "{output.features_file}" ; ln -s "{input.matrix_file_no_indexHopping}" "{output.matrix_file}" ; ln -s "{input.barcodes_file_no_indexHopping}" "{output.barcodes_file}"'


# Create symlinks from preprocessing (novaSeq) folder to analysis folder, to be able to continue with downstream analyses
# This rule creates a symlink of a symlink
rule create_symlinks:
    input:
        features_file_tmp = CELLRANGER_NOVA_OUT + '{sample}.features.tsv',
        matrix_file_tmp = CELLRANGER_NOVA_OUT +'{sample}.matrix.mtx',
        barcodes_file_tmp = CELLRANGER_NOVA_OUT + '{sample}.barcodes.tsv',
        web_summary_tmp = CELLRANGER_NOVA_OUT + '{sample}.web_summary.html',
        metrics_summary_tmp = CELLRANGER_NOVA_OUT + '{sample}.metrics_summary.csv'
    output:
        # path to sample-specific analysis folder (determined by the given sample)
        # NOTE: here, a fixed structure is assumed for the analysis directory
        features_file = '{root}{crsample}/singlecell_rna/analysis/cellranger_run/{sample}.features.tsv',
        matrix_file = '{root}{crsample}/singlecell_rna/analysis/cellranger_run/{sample}.matrix.mtx',
        barcodes_file = '{root}{crsample}/singlecell_rna/analysis/cellranger_run/{sample}.barcodes.tsv',
        web_summary = '{root}{crsample}/singlecell_rna/analysis/cellranger_run/{sample}.web_summary.html',
        metrics_summary = '{root}{crsample}/singlecell_rna/analysis/cellranger_run/{sample}.metrics_summary.csv'
    params:
        sample_out = '{root}{crsample}/singlecell_rna/analysis/cellranger_run/',
        lsfoutfile = '{root}{crsample}/singlecell_rna/analysis/cellranger_run/{sample}.create_symlinks.lsfout.log',
        lsferrfile = '{root}{crsample}/singlecell_rna/analysis/cellranger_run/{sample}.create_symlinks.lsferr.log',
        scratch = config['tools']['create_symlinks']['scratch'],
        mem = config['tools']['create_symlinks']['mem'],
        time = config['tools']['create_symlinks']['time']
    threads:
        config['tools']['create_symlinks']['threads']
    benchmark:
        '{root}{crsample}/singlecell_rna/analysis/cellranger_run/{sample}.create_symlinks.benchmark'
    shell:
        'mkdir -p {params.sample_out} ; ln -s "{input.features_file_tmp}" "{output.features_file}" ; ln -s "{input.matrix_file_tmp}" "{output.matrix_file}" ; ln -s "{input.barcodes_file_tmp}" "{output.barcodes_file}" ; ln -s "{input.web_summary_tmp}" "{output.web_summary}" ; ln -s "{input.metrics_summary_tmp}" "{output.metrics_summary}"'
