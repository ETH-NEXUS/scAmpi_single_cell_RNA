# cellranger call to process the raw samples
rule cellranger_count: 
    input:
        fastqs_dir = config['inputOutput']['input_fastqs'],
        reference = config['resources']['reference_transcriptome']
    output:
        features_file = 'results/cellranger_run/{sample}.features.tsv',
        matrix_file = 'results/cellranger_run/{sample}.matrix.mtx',
        barcodes_file = 'results/cellranger_run/{sample}.barcodes.tsv'
    params:
        cr_out = 'results/cellranger_run/',
        local_cores = config['tools']['cellranger_count']['local_cores'],
        lsfoutfile = 'results/cellranger_run/{sample}.cellranger_count.lsfout.log',
        lsferrfile = 'results/cellranger_run/{sample}.cellranger_count.lsferr.log',
        scratch = config['tools']['cellranger_count']['scratch'],
        mem = config['tools']['cellranger_count']['mem'],
        time = config['tools']['cellranger_count']['time'],
        variousParams = config['tools']['cellranger_count']['variousParams'],
        metrics_summary = 'results/cellranger_run/{sample}.metrics_summary.csv',
        web_summary = 'results/cellranger_run/{sample}.web_summary.html',
	mySample = '{sample}' # needs to be the prefix of all fastq files that belong to this sample. NOTE: no dots are allowed in sample names!
    threads:
        config['tools']['cellranger_count']['threads']
    benchmark:
        'results/cellranger_run/{sample}.cellranger_count.benchmark'
    # NOTE: cellranger count function cannot specify the output directory, the output it the path you call it from.
    # Therefore, a subshell is used here.
    shell:
	    '(cd {params.cr_out}; {config[tools][cellranger_count][call]} count --id={params.mySample} --sample={params.mySample} --transcriptome={input.reference} --localcores={params.local_cores} --fastqs={input.fastqs_dir} --nosecondary {params.variousParams}); gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv.gz ; gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ; gunzip {params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz ; ln -s "{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/features.tsv" "{output.features_file}"; ln -s "{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/matrix.mtx" "{output.matrix_file}"; ln -s "{params.cr_out}{params.mySample}/outs/filtered_feature_bc_matrix/barcodes.tsv" "{output.barcodes_file}" ; ln -s "{params.cr_out}{params.mySample}/outs/web_summary.html" "{params.web_summary}" ; ln -s "{params.cr_out}{params.mySample}/outs/metrics_summary.csv" "{params.metrics_summary}"'

# create hdf5 from count matrix, genes, and cell barcodes file
rule create_hdf5:
    input:
        genes_file = 'results/cellranger_run/{sample}.features.tsv',
        matrix_file = 'results/cellranger_run/{sample}.matrix.mtx',
        barcodes_file = 'results/cellranger_run/{sample}.barcodes.tsv'
    output:
        outfile = 'results/rawCounts/{sample}.h5'
    params:
        lsfoutfile = 'results/rawCounts/{sample}.create_hd5.lsfout.log',
        lsferrfile = 'results/rawCounts/{sample}.create_hd5.lsferr.log',
        scratch = config['tools']['create_hd5']['scratch'],
        mem = config['tools']['create_hd5']['mem'],
        time = config['tools']['create_hd5']['time']
    threads:
        config['tools']['create_hd5']['threads']
    benchmark:
        'results/rawCounts/{sample}.create_hd5.benchmark'
    shell:
        '{config[tools][create_hd5][call]} -g {input.genes_file} -m {input.matrix_file} -b {input.barcodes_file} -o {output.outfile}'

# identify doublets with scDblFinder
rule identify_doublets:
    input:
        infile = 'results/rawCounts/{sample}.h5'
    output:
        outfile = 'results/filteredCounts/{sample}.doublet_barcodes.txt'
    params:
        sample = '{sample}',
        outdir = 'results/filteredCounts/',
        lsfoutfile = 'results/filteredCounts/{sample}.identify_doublets.lsfout.log',
        lsferrfile = 'results/filteredCounts/{sample}.identify_doublets.lsferr.log',
        scratch = config['tools']['identify_doublets']['scratch'],
        mem = config['tools']['identify_doublets']['mem'],
        time = config['tools']['identify_doublets']['time']
    threads:
        config['tools']['identify_doublets']['threads']
    benchmark:
        'results/filteredCounts/{sample}.identify_doublets.benchmark'
    shell:
        "{config[tools][identify_doublets][call]} " +
        "--hdf5File {input.infile} " +
        "--sample {params.sample} " +
        "--outdir {params.outdir}"


# filter out whole genes and cells.
# Filter out cells if the fraction of reads mapped to MT genes is too high or the NODG is too low.
# Filter out genes if they are not protein coding, if they are mitochondrial genes, or if they encode for ribosomal proteins.
rule filter_genes_and_cells:
    input:
        infile = 'results/rawCounts/{sample}.h5',
        doublets = 'results/filteredCounts/{sample}.doublet_barcodes.txt'
    output:
        outfile = 'results/filteredCounts/{sample}.genes_cells_filtered.h5'
    params:
        nmads_fractionMT = config['tools']['filter_genes_and_cells']['nmads_fractionMT'],
        nmads_NODG = config['tools']['filter_genes_and_cells']['nmads_NODG'],
        threshold_fractionMT = config['tools']['filter_genes_and_cells']['threshold_fractionMT'],
        threshold_NODG = config['tools']['filter_genes_and_cells']['threshold_NODG'],
        remove_doublets = config['tools']['filter_genes_and_cells']['remove_doublets'],
        outDir = 'results/filteredCounts/',
        lsfoutfile = 'results/filteredCounts/{sample}.filter_genes_and_cells.lsfout.log',
        lsferrfile = 'results/filteredCounts/{sample}.filter_genes_and_cells.lsferr.log',
        scratch = config['tools']['filter_genes_and_cells']['scratch'],
        mem = config['tools']['filter_genes_and_cells']['mem'],
        time = config['tools']['filter_genes_and_cells']['time'],
        genomeVersion = config['tools']['filter_genes_and_cells']['genomeVersion'],
        sample = '{sample}'
    threads:
        config['tools']['filter_genes_and_cells']['threads']
    benchmark:
        'results/filteredCounts/{sample}.filter_genes_and_cells.benchmark'
    shell:
            '{config[tools][filter_genes_and_cells][call]} ' +
            '--hdf5File {input.infile} ' +
            '--nmads_NODG {params.nmads_NODG} ' +
            '--nmads_fractionMT {params.nmads_fractionMT} ' +
            '--threshold_NODG {params.threshold_NODG} ' +
            '--threshold_fractionMT {params.threshold_fractionMT} ' +
            '--genomeVersion {params.genomeVersion} ' +
            '--doublet_barcodes {input.doublets} ' +
            '--remove_doublets {params.remove_doublets} ' +
            '--sample {params.sample} ' +
            '--outDir {params.outDir}'



# perform normalisation, cell cycle correction and other preprocessing using sctransform
rule sctransform_preprocessing:
    input:
        hdf5_file =  'results/filteredCounts/{sample}.h5',
    output:
        outfile = 'results/counts_corrected/{sample}.corrected.RDS',
        highly_variable = 'results/counts_corrected/{sample}.corrected.variable_genes.h5',
    params:
        lsfoutfile = 'results/counts_corrected/{sample}.corrected.lsfout.log',
        lsferrfile = 'results/counts_corrected/{sample}.corrected.lsferr.log',
        scratch = config['tools']['sctransform_preprocessing']['scratch'],
        mem = config['tools']['sctransform_preprocessing']['mem'],
        time = config['tools']['sctransform_preprocessing']['time'],
        sample = '{sample}',
        number_genes = config['tools']['sctransform_preprocessing']['number_genes'],
        min_var = config['tools']['sctransform_preprocessing']['min_var'],
        n_nn = config['tools']['sctransform_preprocessing']['n_nn'],
        outDir = 'results/counts_corrected/',
    threads:
        config['tools']['sctransform_preprocessing']['threads']
    benchmark:
        'results/counts_corrected/{sample}.corrected.benchmark'
    shell:
        "{config[tools][sctransform_preprocessing][call]} --inHDF5 {input.hdf5_file} --sample {params.sample} --number_genes {params.number_genes} --min_var {params.min_var} --n_nn {params.n_nn} --outdir {params.outDir} "


# perform clustering with phenograph
rule phenograph:
    input:
        infile = 'results/counts_corrected/{sample}.variable_genes.h5'
    output:
        outfile = 'results/clustering/{sample}.clusters_phenograph.csv',
        distance_matrix = 'results/clustering/{sample}.distance_matrix.tsv',
        modularity_score = 'results/clustering/{sample}.modularity_score.txt'
    params:
        n_neighbours = config['tools']['clustering']['phenograph']['n_neighbours'],
        min_cluster_size = config['tools']['clustering']['phenograph']['min_cluster_size'],
        log_normalize = config['tools']['clustering']['phenograph']['log_normalize'],
        lsfoutfile = 'results/clustering/{sample}.phenograph.lsfout.log',
        lsferrfile = 'results/clustering/{sample}.phenograph.lsferr.log',
        scratch = config['tools']['clustering']['phenograph']['scratch'],
        mem = config['tools']['clustering']['phenograph']['mem'],
        time = config['tools']['clustering']['phenograph']['time']
    threads:
        config['tools']['clustering']['phenograph']['threads']
    benchmark:
        'results/clustering/{sample}.phenograph.benchmark'
    shell:
        '{config[tools][clustering][phenograph][call]} ' +
        '--input_file {input.infile} ' +
        '--output_file {output.outfile} ' +
        '--distance_matrix {output.distance_matrix} ' +
        '--modularity_score {output.modularity_score} ' +
        '--n_neighbours {params.n_neighbours} ' +
        '--min_size {params.min_cluster_size} ' +
        '--log_normalize {params.log_normalize} ' +
        '--n_threads {threads}'


# prepare sce object in RDS file for cell type classification
rule prepare_celltyping:
    input:
        RDS_file = 'results/counts_corrected/{sample}.RDS',
        cluster = 'results/clustering/{sample}.clusters_phenograph.csv',
        distanceMatrix = 'results/clustering/{sample}.distance_matrix.tsv',
        modularity_score = 'results/clustering/{sample}.modularity_score.txt',
    output:
        outfile = 'results/prep_celltyping/{sample}.RDS'
    params:
        lsfoutfile = 'results/prep_celltyping/{sample}.prepare_celltyping.lsfout.log',
        lsferrfile = 'results/prep_celltyping/{sample}.prepare_celltyping.lsferr.log',
        scratch = config['tools']['prepare_celltyping']['scratch'],
        mem = config['tools']['prepare_celltyping']['mem'],
        time = config['tools']['prepare_celltyping']['time'],
        outputDirec = 'results/prep_celltyping/',
        sampleName = '{sample}',
    threads:
        config['tools']['prepare_celltyping']['threads']
    benchmark:
        'results/prep_celltyping/{sample}.prepare_celltyping.benchmark'
    shell:
        "{config[tools][prepare_celltyping][call]} --in_sce {input.RDS_file} --phenograph_cluster {input.cluster} --outputDirec {params.outputDirec} --sampleName {params.sampleName} --distanceMatrix {input.distanceMatrix} --modularity_score {input.modularity_score} "


# perform cell type classification
rule cell_type_classification:
    input:
        infile = 'results/prep_celltyping/{sample}.RDS',
    output:
        outfile = 'results/celltype_classification/{sample}.phenograph_celltype_association.txt',
        out_sce = 'results/celltype_classification/{sample}.RDS'
    params:
        lsfoutfile = 'results/celltype_classification/{sample}.cell_type_classification.lsfout.log',
        lsferrfile = 'results/celltype_classification/{sample}.cell_type_classification.lsferr.log',
        scratch = config['tools']['cell_type_classification']['scratch'],
        mem = config['tools']['cell_type_classification']['mem'],
        time = config['tools']['cell_type_classification']['time'],
        min_genes = config['tools']['cell_type_classification']['min_genes'],
        celltype_lists = config['resources']['celltype_lists'],
        celltype_config = config['resources']['celltype_config'],
        outputDirec = 'results/celltype_classification/',
        sampleName = '{sample}',
    threads:
        config['tools']['cell_type_classification']['threads']
    benchmark:
        'results/celltype_classification/{sample}.cell_type_classification.benchmark'
    shell:
        '{config[tools][cell_type_classification][call]} ' +
        '--SCE {input.infile} ' +
        '--celltype_lists {params.celltype_lists} ' +
        '--celltype_config {params.celltype_config} ' +
        '--sampleName {params.sampleName} ' +
        '--min_genes {params.min_genes} ' +
        '--outputDirec {params.outputDirec} '


# filter out atypical cells from sce object
rule remove_atypical:
    input:
        infile = 'results/celltype_classification/{sample}.RDS',
        cluster_table = 'results/celltype_classification/{sample}.phenograph_celltype_association.txt',
    output:
        out_sce = 'results/atypical_removed/{sample}.atypical_removed.RDS',
        out_table = 'results/atypical_removed/{sample}.atypical_removed.phenograph_celltype_association.txt'
    params:
        lsfoutfile = 'results/atypical_removed/{sample}.atypical_removed.lsfout.log',
        lsferrfile = 'results/atypical_removed/{sample}.atypical_removed.lsferr.log',
        scratch = config['tools']['remove_atypical']['scratch'],
        mem = config['tools']['remove_atypical']['mem'],
        time = config['tools']['remove_atypical']['time'],
        celltype_config = config['resources']['celltype_config'],
        outputDirec = 'results/atypical_removed/',
        sample_name = '{sample}',
        threshold_filter = config['tools']['remove_atypical']['threshold_filter'],
        min_threshold = config['tools']['remove_atypical']['min_threshold'],
        threshold_type = config['tools']['remove_atypical']['threshold_type'],
    threads:
        config['tools']['remove_atypical']['threads']
    benchmark:
        'results/atypical_removed/{sample}.atypical_removed.benchmark'
    shell:
        "{config[tools][remove_atypical][call]} --sce_in {input.infile} --cluster_table {input.cluster_table} --celltype_config {params.celltype_config} --threshold_filter {params.threshold_filter} --min_threshold {params.min_threshold} --threshold_type {params.threshold_type} --outDir {params.outputDirec} --sample_name {params.sample_name} "


# perform gsva gene set analysis
rule gsva:
    input:
        infile = 'results/atypical_removed/{sample}.RDS',
    output:
        outfile = 'results/gsva/{sample}.gsetscore_hm.png',
    params:
        lsfoutfile = 'results/gsva/{sample}.gsva.lsfout.log',
        lsferrfile = 'results/gsva/{sample}.gsva.lsferr.log',
        scratch = config['tools']['gsva']['scratch'],
        mem = config['tools']['gsva']['mem'],
        time = config['tools']['gsva']['time'],
        outputDirec = 'results/gsva/',
        sampleName = '{sample}',
        genesets = config['resources']['genesets'],
    threads:
        config['tools']['gsva']['threads']
    benchmark:
        'results/gsva/{sample}.gsva.benchmark'
    shell:
        "{config[tools][gsva][call]} --SCE {input.infile} --geneset {params.genesets} --outputDirec {params.outputDirec} --sampleName {params.sampleName} "


# generate plots about sample composition and gene expression
rule plotting:
    input:
        infile = 'results/atypical_removed/{sample}.RDS',
    output:
        outfile = 'results/plotting/{sample}.celltype_barplot.png',
    params:
        lsfoutfile = 'results/plotting/{sample}.plotting.lsfout.log',
        lsferrfile = 'results/plotting/{sample}.plotting.lsferr.log',
        scratch = config['tools']['plotting']['scratch'],
        mem = config['tools']['plotting']['mem'],
        time = config['tools']['plotting']['time'],
        outputDirec = 'results/plotting/',
        sampleName = '{sample}',
        genes_of_interest = config['resources']['priority_genes'],
        colour_config = config['resources']['colour_config'],
        use_alias = config['tools']['plotting']['use_alias']
    threads:
        config['tools']['plotting']['threads']
    benchmark:
        'results/plotting/{sample}.plotting.benchmark'
    shell:
        "{config[tools][plotting][call]} --sce_in {input.infile} --genelist {params.genes_of_interest} --outDir {params.outputDirec} --sampleName {params.sampleName} --colour_config {params.colour_config} --toggle_label {params.use_alias}"


# adapt this directory in master snake file to prevent recomputing the cohort in each analysis
### THIS SHOULD PROBABLY BE CHANGED WITH NEW STRUCTURE. CLEARIFY! TODO 

# perform assembly of the nonmalignant reference cohort
rule assemble_nonmalignant_cohort:
    input:
        inputDir_hdf5 = config['tools']['assemble_non_malignant_reference']['hdf5_dir'],
        inputDir_classify = config['tools']['assemble_non_malignant_reference']['celltype_dir']
    output:
        outfile = 'results/non_malignant_reference/nonmalignant_reference_cohort.h5'
    params:
        outDir = 'results/non_malignant_reference/',
        non_malignant_types = config['tools']['assemble_non_malignant_reference']['non_malignant_types'],
        lsfoutfile = 'results/non_malignant_reference/nonmalignant_reference_cohort.lsfout.log',
        lsferrfile = 'results/non_malignant_reference/nonmalignant_reference_cohort.lsferr.log',
        scratch = config['tools']['assemble_non_malignant_reference']['scratch'],
        mem = config['tools']['assemble_non_malignant_reference']['mem'],
        time = config['tools']['assemble_non_malignant_reference']['time']
    threads:
        config['tools']['assemble_non_malignant_reference']['threads']
    benchmark:
        'results/non_malignant_reference/nonmalignant_reference_cohort.benchmark'
    shell:
        "{config[tools][assemble_non_malignant_reference][call]} --hdf5_dir {input.inputDir_hdf5} --celltype_dir {input.inputDir_classify} --non_malignant_celltypes {params.non_malignant_types} --out_dir {params.outDir} "

# plot reference cohort tSNEs
rule plot_tSNEs_nonmalignant_cohort:
    input:
        inputHDF5 = 'results/non_malignant_reference/nonmalignant_reference_cohort.h5'
    output:
        outfile_batch = 'results/non_malignant_reference/nonmalignant_reference_cohort.tSNE_batch.png',
	outfile_ct = 'results/non_malignant_reference/nonmalignant_reference_cohort.tSNE_celltype.png'
    params:
        lsfoutfile = 'results/non_malignant_reference/plot_tSNEs_nonmalignant.lsfout.log',
        lsferrfile = 'results/non_malignant_reference/plot_tSNEs_nonmalignant.lsferr.log',
        scratch = config['tools']['plot_tSNE_nonmalignant']['scratch'],
        mem = config['tools']['plot_tSNE_nonmalignant']['mem'],
        time = config['tools']['plot_tSNE_nonmalignant']['time']
    threads:
        config['tools']['plot_tSNE_nonmalignant']['threads']
    benchmark:
        'results/non_malignant_reference/plot_tSNEs_nonmalignant.benchmark'
    shell:
        "{config[tools][plot_tSNE_nonmalignant][call]} --hdf5File {input.inputHDF5} --outFile_batch {output.outfile_batch} --outFile_ct {output.outfile_ct} "


# perform the differential expression analysis using a Wilcoxon test
rule diff_exp_genes:
    input:
        sce_in = 'results/atypical_removed/{sample}.RDS',
        cell_types = 'results/atypical_removed/{sample}.phenograph_celltype_association.txt'
    output:
        #outpath = dynamic('results/diff_exp/' + "{sample}.{clusterid}.DEgenes.tsv"),
        success = 'results/diff_exp/{sample}.diffExp_success.txt'
    params:
        sampleName = '{sample}',
        malignant = config['inputOutput']['malignant_cell_type'],
        threshold_comparison = config['tools']['diff_exp']['threshold_comparison'],
        fdr_cut = config['tools']['diff_exp']['fdr_cut'],
        fc_cut = config['tools']['diff_exp']['fc_cut'],
        mindiff2second = config['tools']['diff_exp']['mindiff2second'],
        minNumberNonMalignant = config['tools']['diff_exp']['minNumberNonMalignant'],
        lsfoutfile = 'results/diff_exp/{sample}.diff_exp.lsfout.log',
        lsferrfile = 'results/diff_exp/{sample}.diff_exp.lsferr.log',
        scratch = config['tools']['diff_exp']['scratch'],
        mem = config['tools']['diff_exp']['mem'],
        time = config['tools']['diff_exp']['time'],
        outpath = 'results/diff_exp/',
    threads:
        config['tools']['diff_exp']['threads']
    benchmark:
        'results/diff_exp/{sample}.diff_exp.benchmark'
    shell:
        "{config[tools][diff_exp][call]} --sample_data {input.sce_in} " +
        "--sampleName {params.sampleName} " +
        "--cluster_table {input.cell_types} " +
        "--malignant_tag {params.malignant} " +
        "--fdr_cut {params.fdr_cut} " +
        "--fc_cut {params.fc_cut} " +
        "--mindiff2second {params.mindiff2second} " +
        "--threshold_comparison {params.threshold_comparison} " +
        "--minNumberNonMalignant {params.minNumberNonMalignant} " +
        "--outdir {params.outpath} ; " +
        "date > {output.success} "


# give out gene expression values per cluster
rule gene_exp:
    input:
        sce_in = 'results/atypical_removed/{sample}.RDS',
    output:
        out = 'results/gene_exp/{sample}.gene_expression_per_cluster.tsv'
    params:
        sampleName = '{sample}',
        lsfoutfile = 'results/gene_exp/{sample}.gene_exp.lsfout.log',
        lsferrfile = 'results/gene_exp/{sample}.gene_exp.lsferr.log',
        scratch = config['tools']['gene_exp']['scratch'],
        mem = config['tools']['gene_exp']['mem'],
        time = config['tools']['gene_exp']['time'],
        outpath = 'results/gene_exp/',
        threshold_sample = config['tools']['gene_exp']['threshold_sample'],
        type_sample = config['tools']['gene_exp']['type_sample'],
        priority_genes = config['resources']['priority_genes'],
    threads:
        config['tools']['gene_exp']['threads']
    benchmark:
        'results/gene_exp/{sample}.gene_exp.benchmark'
    shell:
        "{config[tools][gene_exp][call]} --sce_in {input.sce_in} --priority_genes {params.priority_genes} --filtering_threshold_sample {params.threshold_sample} --filter_type_sample {params.type_sample} --outDir {params.outpath} --sample_name {params.sampleName} "


# This rule generates general quality control plots to hdf5 expression files
rule generate_qc_plots :
    input:
        infile = '{sample}.h5'
    output:
        out = '{sample}.h5.histogram_library_sizes.png'
    params:
        lsfoutfile = '{sample}.generate_qc_plots.lsfout.log',
        lsferrfile = '{sample}.generate_qc_plots.lsferr.log',
        scratch = config['tools']['generate_qc_plots']['scratch'],
        mem = config['tools']['generate_qc_plots']['mem'],
        time = config['tools']['generate_qc_plots']['time'],
    threads:
        config['tools']['generate_qc_plots']['threads']
    benchmark:
        '{sample}.generate_qc_plots.benchmark'
    shell:
        '{config[tools][generate_qc_plots][call]} --hdf5File {input.infile} '


# This rule creates a box plot comparing cell type fractions across samples
rule generate_cell_type_boxplot:
    input:
        previous_samples = config['resources']['previous_samples'],
        sample_cell_types = 'results/atypical_removed/{sample}.phenograph_celltype_association.txt'
    output:
        out = 'results/plotting/{sample}.boxplot_cell_types_cohort.png'
    params:
        sampleName = '{sample}',
        sampleName_short = config['tools']['cellranger_count']['cellranger_sampleName'],
        outDir = 'results/plotting/',
        lsfoutfile = 'results/plotting/{sample}.boxplot_cell_types_cohort.lsfout.log',
        lsferrfile = 'results/plotting/{sample}.boxplot_cell_types_cohort.lsferr.log',
        scratch = config['tools']['generate_cell_type_boxplot']['scratch'],
        mem = config['tools']['generate_cell_type_boxplot']['mem'],
        time = config['tools']['generate_cell_type_boxplot']['time'],
    threads:
        config['tools']['generate_cell_type_boxplot']['threads']
    benchmark:
        'results/plotting/{sample}.boxplot_cell_types_cohort.benchmark'
    shell:
        '{config[tools][generate_cell_type_boxplot][call]} --previous_samples {input.previous_samples} --current_sample {input.sample_cell_types} --sampleName {params.sampleName} --sampleName_short {params.sampleName_short} --outDir {params.outDir}'



# This rule integrates samples of the cohort and visualizes the integration with UMAPs
rule sample_integration:
    input:
        previous_samples = config['resources']['previous_samples_counts'],
        current_sample = 'results/atypical_removed/{sample}.RDS'
    output:
        out = 'results/plotting/{sample}.sample_integration_highlight_current.png'
    params:
        sampleName = '{sample}',
        outDir = 'results/plotting/',
        lsfoutfile = 'results/plotting/{sample}.sample_integration.lsfout.log',
        lsferrfile = 'results/plotting/{sample}.sample_integration.lsferr.log',
        scratch = config['tools']['sample_integration']['scratch'],
        mem = config['tools']['sample_integration']['mem'],
        time = config['tools']['sample_integration']['time'],
        sampleName_short = config['tools']['cellranger_count']['cellranger_sampleName'],
	colour_config = config['resources']['colour_config']
    threads:
        config['tools']['sample_integration']['threads']
    benchmark:
        'results/plotting/{sample}.sample_integration.benchmark'
    shell:
        '{config[tools][sample_integration][call]} --cohort_list {input.previous_samples} --sample_data {input.current_sample} --sampleName {params.sampleName} --sampleName_short {params.sampleName_short} --colour_config {params.colour_config} --outdir {params.outDir}'



# calculate for each cluster the number of cells it countains and the percentage of all cells
rule cellPercentInCluster:
    input:
        clusterCsv = 'results/atypical_removed/{sample}.phenograph_celltype_association.txt'
    output:
#    print(drugID)
        out = 'results/clusterpercent/{sample}.clusters_cell_count_percent.txt'
    params:
        lsfoutfile = 'results/clusterpercent/{sample}.clusterPercent.lsfout.log',
        lsferrfile =  'results/clusterpercent/{sample}.clusterPercent.lsferr.log',
        scratch = config['tools']['cellPercentInCluster']['scratch'],
        mem = config['tools']['cellPercentInCluster']['mem'],
        time = config['tools']['cellPercentInCluster']['time'],
        variousParams = config['tools']['cellPercentInCluster']['variousParams']
    threads:
        config['tools']['cellPercentInCluster']['threads']
    benchmark:
        'results/clusterpercent/{sample}.clusterPercent.benchmark'
    shell:
        '{config[tools][cellPercentInCluster][call]} --inputTable {input.clusterCsv} --outFile {output.out} {params.variousParams}'