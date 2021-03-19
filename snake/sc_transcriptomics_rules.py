if not 'CELLRANGER_IN' in globals():
    CELLRANGER_IN = INPUTDIR
if not 'CELLRANGER_OUT' in globals():
    CELLRANGER_OUT = OUTDIR + 'cellranger_run/'
    #CR_MATRIX_PATH = '/outs/filtered_gene_bc_matrices/'+config['resources']['transcriptome_code']+'/'
# cellranger call to process the raw samples
rule cellranger_count: 
    input:
        fastqs_dir = CELLRANGER_IN, #+ '{sample}',
        reference = config['resources']['reference_transcriptome']
    output:
        features_file = CELLRANGER_OUT + '{sample}.features.tsv',
        matrix_file = CELLRANGER_OUT +'{sample}.matrix.mtx',
        barcodes_file = CELLRANGER_OUT + '{sample}.barcodes.tsv'
    params:
        cr_out = CELLRANGER_OUT,
        local_cores = config['tools']['cellranger_count']['local_cores'],
        lsfoutfile = CELLRANGER_OUT + '{sample}.cellranger_count.lsfout.log',
        lsferrfile = CELLRANGER_OUT + '{sample}.cellranger_count.lsferr.log',
        scratch = config['tools']['cellranger_count']['scratch'],
        mem = config['tools']['cellranger_count']['mem'],
        time = config['tools']['cellranger_count']['time'],
        variousParams = config['tools']['cellranger_count']['variousParams'],
        cellranger_sampleName = config['tools']['cellranger_count']['cellranger_sampleName'],
        metrics_summary = CELLRANGER_OUT + '{sample}.metrics_summary.csv',
        web_summary = CELLRANGER_OUT + '{sample}.web_summary.html'
    threads:
        config['tools']['cellranger_count']['threads']
    benchmark:
        CELLRANGER_OUT + '{sample}.cellranger_count.benchmark'
    # NOTE: cellranger count function cannot specify the output directory, the output it the path you call it from.
    # Therefore, a subshell is used here.
    shell:
        '(cd {params.cr_out}; {config[tools][cellranger_count][call]} count --id={params.cellranger_sampleName} --transcriptome={input.reference} --localcores={params.local_cores} --fastqs={input.fastqs_dir} --nosecondary {params.variousParams}); gunzip {params.cr_out}{params.cellranger_sampleName}/outs/filtered_feature_bc_matrix/features.tsv.gz ; gunzip {params.cr_out}{params.cellranger_sampleName}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz ; gunzip {params.cr_out}{params.cellranger_sampleName}/outs/filtered_feature_bc_matrix/matrix.mtx.gz ; ln -s "{params.cr_out}{params.cellranger_sampleName}/outs/filtered_feature_bc_matrix/features.tsv" "{output.features_file}"; ln -s "{params.cr_out}{params.cellranger_sampleName}/outs/filtered_feature_bc_matrix/matrix.mtx" "{output.matrix_file}"; ln -s "{params.cr_out}{params.cellranger_sampleName}/outs/filtered_feature_bc_matrix/barcodes.tsv" "{output.barcodes_file}" ; ln -s "{params.cr_out}{params.cellranger_sampleName}/outs/web_summary.html" "{params.web_summary}" ; ln -s "{params.cr_out}{params.cellranger_sampleName}/outs/metrics_summary.csv" "{params.metrics_summary}"'


if not 'CREATEHD5_IN' in globals():
    CREATEHD5_IN = CELLRANGER_OUT
if not 'CREATEHD5_OUT' in globals():
    CREATEHD5_OUT = OUTDIR + 'rawCounts/'

# create hdf5 from count matrix, genes, and cell barcodes file
rule create_hdf5:
    input:
        genes_file = CREATEHD5_IN + '{sample}.features.tsv',
        matrix_file = CREATEHD5_IN +'{sample}.matrix.mtx',
        barcodes_file = CREATEHD5_IN + '{sample}.barcodes.tsv'
    output:
        outfile = CREATEHD5_OUT + '{sample}.h5'
    params:
        lsfoutfile = CREATEHD5_OUT + '{sample}.create_hd5.lsfout.log',
        lsferrfile = CREATEHD5_OUT + '{sample}.create_hd5.lsferr.log',
        scratch = config['tools']['create_hd5']['scratch'],
        mem = config['tools']['create_hd5']['mem'],
        time = config['tools']['create_hd5']['time']
    threads:
        config['tools']['create_hd5']['threads']
    benchmark:
        CREATEHD5_OUT + '{sample}.create_hd5.benchmark'
    shell:
        '{config[tools][create_hd5][call]} -g {input.genes_file} -m {input.matrix_file} -b {input.barcodes_file} -o {output.outfile}'

if not 'FILTER_GENES_CELLS_IN' in globals():
    FILTER_GENES_CELLS_IN = CREATEHD5_OUT
if not 'FILTER_GENES_CELLS_OUT' in globals():
    FILTER_GENES_CELLS_OUT = OUTDIR + 'filteredCounts/'


# identify doublets with scDblFinder
rule identify_doublets:
    input:
        infile = FILTER_GENES_CELLS_IN + '{sample}.h5'
    output:
        outfile = FILTER_GENES_CELLS_OUT + '{sample}.doublet_barcodes.txt'
    params:
        sample = '{sample}',
        outdir = FILTER_GENES_CELLS_OUT,
        lsfoutfile = FILTER_GENES_CELLS_OUT + '{sample}.identify_doublets.lsfout.log',
        lsferrfile = FILTER_GENES_CELLS_OUT + '{sample}.identify_doublets.lsferr.log',
        scratch = config['tools']['identify_doublets']['scratch'],
        mem = config['tools']['identify_doublets']['mem'],
        time = config['tools']['identify_doublets']['time']
    threads:
        config['tools']['identify_doublets']['threads']
    benchmark:
        FILTER_GENES_CELLS_OUT + '{sample}.identify_doublets.benchmark'
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
        infile = FILTER_GENES_CELLS_IN + '{sample}.h5',
        doublets = FILTER_GENES_CELLS_OUT + '{sample}.doublet_barcodes.txt'
    output:
        outfile = FILTER_GENES_CELLS_OUT + '{sample}.genes_cells_filtered.h5'
    params:
        nmads_fractionMT = config['tools']['filter_genes_and_cells']['nmads_fractionMT'],
        nmads_NODG = config['tools']['filter_genes_and_cells']['nmads_NODG'],
        threshold_fractionMT = config['tools']['filter_genes_and_cells']['threshold_fractionMT'],
        threshold_NODG = config['tools']['filter_genes_and_cells']['threshold_NODG'],
        outDir = FILTER_GENES_CELLS_OUT,
        lsfoutfile = FILTER_GENES_CELLS_OUT + '{sample}.filter_genes_and_cells.lsfout.log',
        lsferrfile = FILTER_GENES_CELLS_OUT + '{sample}.filter_genes_and_cells.lsferr.log',
        scratch = config['tools']['filter_genes_and_cells']['scratch'],
        mem = config['tools']['filter_genes_and_cells']['mem'],
        time = config['tools']['filter_genes_and_cells']['time'],
        genomeVersion = config['tools']['filter_genes_and_cells']['genomeVersion'],
        sample = '{sample}'
    threads:
        config['tools']['filter_genes_and_cells']['threads']
    benchmark:
        FILTER_GENES_CELLS_OUT + '{sample}.filter_genes_and_cells.benchmark'
    shell:
            '{config[tools][filter_genes_and_cells][call]} ' +
            '--hdf5File {input.infile} ' +
            '--nmads_NODG {params.nmads_NODG} ' +
            '--nmads_fractionMT {params.nmads_fractionMT} ' +
            '--threshold_NODG {params.threshold_NODG} ' +
            '--threshold_fractionMT {params.threshold_fractionMT} ' +
            '--genomeVersion {params.genomeVersion} ' +
            '--doublet_barcodes {input.doublets} ' +
            '--sample {params.sample} ' +
            '--outDir {params.outDir}'


if not 'SCTRANSFORM_IN' in globals():
    SCTRANSFORM_IN = FILTER_GENES_CELLS_OUT
if not 'SCTRANSFORM_OUT' in globals():
    SCTRANSFORM_OUT = OUTDIR + 'counts_corrected/'

# perform normalisation, cell cycle correction and other preprocessing using sctransform
rule sctransform_preprocessing:
    input:
        hdf5_file =  SCTRANSFORM_IN + '{sample}.h5',
    output:
        outfile = SCTRANSFORM_OUT + '{sample}.corrected.RDS',
        highly_variable = SCTRANSFORM_OUT + '{sample}.corrected.variable_genes.h5',
    params:
        lsfoutfile = SCTRANSFORM_OUT + '{sample}.corrected.lsfout.log',
        lsferrfile = SCTRANSFORM_OUT + '{sample}.corrected.lsferr.log',
        scratch = config['tools']['sctransform_preprocessing']['scratch'],
        mem = config['tools']['sctransform_preprocessing']['mem'],
        time = config['tools']['sctransform_preprocessing']['time'],
        sample = '{sample}',
        number_genes = config['tools']['sctransform_preprocessing']['number_genes'],
        min_var = config['tools']['sctransform_preprocessing']['min_var'],
        n_nn = config['tools']['sctransform_preprocessing']['n_nn'],
        outDir = SCTRANSFORM_OUT,
    threads:
        config['tools']['sctransform_preprocessing']['threads']
    benchmark:
        SCTRANSFORM_OUT + '{sample}.corrected.benchmark'
    shell:
        "{config[tools][sctransform_preprocessing][call]} --inHDF5 {input.hdf5_file} --sample {params.sample} --number_genes {params.number_genes} --min_var {params.min_var} --n_nn {params.n_nn} --outdir {params.outDir} "


if not 'PHENOGRAPH_IN' in globals():
    PHENOGRAPH_IN = SCTRANSFORM_OUT
if not 'PHENOGRAPH_OUT' in globals():
    PHENOGRAPH_OUT = OUTDIR + 'clustering/'

# perform clustering with phenograph
rule phenograph:
    input:
        infile = PHENOGRAPH_IN + '{sample}.variable_genes.h5'
    output:
        outfile = PHENOGRAPH_OUT + '{sample}.clusters_phenograph.csv',
        distance_matrix = PHENOGRAPH_OUT + '{sample}.distance_matrix.tsv',
        modularity_score = PHENOGRAPH_OUT + '{sample}.modularity_score.txt'
    params:
        n_neighbours = config['tools']['clustering']['phenograph']['n_neighbours'],
        min_cluster_size = config['tools']['clustering']['phenograph']['min_cluster_size'],
        log_normalize = config['tools']['clustering']['phenograph']['log_normalize'],
        lsfoutfile = PHENOGRAPH_OUT + '{sample}.phenograph.lsfout.log',
        lsferrfile = PHENOGRAPH_OUT + '{sample}.phenograph.lsferr.log',
        scratch = config['tools']['clustering']['phenograph']['scratch'],
        mem = config['tools']['clustering']['phenograph']['mem'],
        time = config['tools']['clustering']['phenograph']['time']
    threads:
        config['tools']['clustering']['phenograph']['threads']
    benchmark:
        PHENOGRAPH_OUT + '{sample}.phenograph.benchmark'
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


if not 'PREPARE_CELLTYPING_IN' in globals():
    PREPARE_CELLTYPING_IN = SCTRANSFORM_OUT
if not 'PREPARE_CELLTYPING_OUT' in globals():
    PREPARE_CELLTYPING_OUT = OUTDIR + 'prep_celltyping/'

# prepare sce object in RDS file for cell type classification
rule prepare_celltyping:
    input:
        RDS_file = PREPARE_CELLTYPING_IN + '{sample}.RDS',
        cluster = PHENOGRAPH_OUT + '{sample}.clusters_phenograph.csv',
        distanceMatrix = PHENOGRAPH_OUT + '{sample}.distance_matrix.tsv',
        modularity_score = PHENOGRAPH_OUT + '{sample}.modularity_score.txt',
    output:
        outfile = PREPARE_CELLTYPING_OUT + '{sample}.RDS'
    params:
        lsfoutfile = PREPARE_CELLTYPING_OUT + '{sample}.prepare_celltyping.lsfout.log',
        lsferrfile = PREPARE_CELLTYPING_OUT + '{sample}.prepare_celltyping.lsferr.log',
        scratch = config['tools']['prepare_celltyping']['scratch'],
        mem = config['tools']['prepare_celltyping']['mem'],
        time = config['tools']['prepare_celltyping']['time'],
        outputDirec = PREPARE_CELLTYPING_OUT,
        sampleName = '{sample}',
    threads:
        config['tools']['prepare_celltyping']['threads']
    benchmark:
        PREPARE_CELLTYPING_OUT + '{sample}.prepare_celltyping.benchmark'
    shell:
        "{config[tools][prepare_celltyping][call]} --in_sce {input.RDS_file} --phenograph_cluster {input.cluster} --outputDirec {params.outputDirec} --sampleName {params.sampleName} --distanceMatrix {input.distanceMatrix} --modularity_score {input.modularity_score} "


if not 'CELLTYPECLASS_IN' in globals():
    CELLTYPECLASS_IN = PREPARE_CELLTYPING_OUT
if not 'CELLTYPECLASS_OUT' in globals():
    CELLTYPECLASS_OUT = OUTDIR + 'celltype_classification/'

# perform cell type classification
rule cell_type_classification:
    input:
        infile = CELLTYPECLASS_IN + '{sample}.RDS',
    output:
        outfile = CELLTYPECLASS_OUT + '{sample}.phenograph_celltype_association.txt',
        out_sce = CELLTYPECLASS_OUT + '{sample}.RDS'
    params:
        lsfoutfile = CELLTYPECLASS_OUT + '{sample}.cell_type_classification.lsfout.log',
        lsferrfile = CELLTYPECLASS_OUT + '{sample}.cell_type_classification.lsferr.log',
        scratch = config['tools']['cell_type_classification']['scratch'],
        mem = config['tools']['cell_type_classification']['mem'],
        time = config['tools']['cell_type_classification']['time'],
        min_genes = config['tools']['cell_type_classification']['min_genes'],
        celltype_lists = config['resources']['celltype_lists'],
        celltype_config = config['resources']['celltype_config'],
        outputDirec = CELLTYPECLASS_OUT,
        sampleName = '{sample}',
    threads:
        config['tools']['cell_type_classification']['threads']
    benchmark:
        CELLTYPECLASS_OUT + '{sample}.cell_type_classification.benchmark'
    shell:
        '{config[tools][cell_type_classification][call]} ' +
        '--SCE {input.infile} ' +
        '--celltype_lists {params.celltype_lists} ' +
        '--celltype_config {params.celltype_config} ' +
        '--sampleName {params.sampleName} ' +
        '--min_genes {params.min_genes} ' +
        '--outputDirec {params.outputDirec} '


if not 'REMOVE_ATYPICAL_IN' in globals():
    REMOVE_ATYPICAL_IN = CELLTYPECLASS_OUT
if not 'REMOVE_ATYPICAL_OUT' in globals():
    REMOVE_ATYPICAL_OUT = OUTDIR + 'atypical_removed/'

# filter out atypical cells from sce object
rule remove_atypical:
    input:
        infile = REMOVE_ATYPICAL_IN + '{sample}.RDS',
        cluster_table = REMOVE_ATYPICAL_IN + '{sample}.phenograph_celltype_association.txt',
    output:
        out_sce = REMOVE_ATYPICAL_OUT + '{sample}.atypical_removed.RDS',
        out_table = REMOVE_ATYPICAL_OUT + '{sample}.atypical_removed.phenograph_celltype_association.txt'
    params:
        lsfoutfile = REMOVE_ATYPICAL_OUT + '{sample}.atypical_removed.lsfout.log',
        lsferrfile = REMOVE_ATYPICAL_OUT + '{sample}.atypical_removed.lsferr.log',
        scratch = config['tools']['remove_atypical']['scratch'],
        mem = config['tools']['remove_atypical']['mem'],
        time = config['tools']['remove_atypical']['time'],
        celltype_config = config['resources']['celltype_config'],
        outputDirec = REMOVE_ATYPICAL_OUT,
        sample_name = '{sample}',
        threshold_filter = config['tools']['remove_atypical']['threshold_filter'],
        min_threshold = config['tools']['remove_atypical']['min_threshold'],
        threshold_type = config['tools']['remove_atypical']['threshold_type'],
    threads:
        config['tools']['remove_atypical']['threads']
    benchmark:
        REMOVE_ATYPICAL_OUT + '{sample}.atypical_removed.benchmark'
    shell:
        "{config[tools][remove_atypical][call]} --sce_in {input.infile} --cluster_table {input.cluster_table} --celltype_config {params.celltype_config} --threshold_filter {params.threshold_filter} --min_threshold {params.min_threshold} --threshold_type {params.threshold_type} --outDir {params.outputDirec} --sample_name {params.sample_name} "


if not 'GSVA_IN' in globals():
    GSVA_IN = REMOVE_ATYPICAL_OUT
if not 'GSVA_OUT' in globals():
    GSVA_OUT = OUTDIR + 'gsva/'

# perform gsva gene set analysis
rule gsva:
    input:
        infile = GSVA_IN + '{sample}.RDS',
    output:
        outfile = GSVA_OUT + '{sample}.gsetscore_hm.png',
    params:
        lsfoutfile = GSVA_OUT + '{sample}.gsva.lsfout.log',
        lsferrfile = GSVA_OUT + '{sample}.gsva.lsferr.log',
        scratch = config['tools']['gsva']['scratch'],
        mem = config['tools']['gsva']['mem'],
        time = config['tools']['gsva']['time'],
        outputDirec = GSVA_OUT,
        sampleName = '{sample}',
        genesets = config['resources']['genesets'],
    threads:
        config['tools']['gsva']['threads']
    benchmark:
        GSVA_OUT + '{sample}.gsva.benchmark'
    shell:
        "{config[tools][gsva][call]} --SCE {input.infile} --geneset {params.genesets} --outputDirec {params.outputDirec} --sampleName {params.sampleName} "


if not 'PLOTTING_IN' in globals():
    PLOTTING_IN = REMOVE_ATYPICAL_OUT
if not 'PLOTTING_OUT' in globals():
    PLOTTING_OUT = OUTDIR + 'plotting/'

# generate plots about sample composition and gene expression
rule plotting:
    input:
        infile = PLOTTING_IN + '{sample}.RDS',
    output:
        outfile = PLOTTING_OUT + '{sample}.celltype_barplot.png',
    params:
        lsfoutfile = PLOTTING_OUT + '{sample}.plotting.lsfout.log',
        lsferrfile = PLOTTING_OUT + '{sample}.plotting.lsferr.log',
        scratch = config['tools']['plotting']['scratch'],
        mem = config['tools']['plotting']['mem'],
        time = config['tools']['plotting']['time'],
        outputDirec = PLOTTING_OUT,
        sampleName = '{sample}',
        genes_of_interest = config['resources']['priority_genes'],
        colour_config = config['resources']['colour_config'],
        use_alias = config['tools']['plotting']['use_alias']
    threads:
        config['tools']['plotting']['threads']
    benchmark:
        PLOTTING_OUT + '{sample}.plotting.benchmark'
    shell:
        "{config[tools][plotting][call]} --sce_in {input.infile} --genelist {params.genes_of_interest} --outDir {params.outputDirec} --sampleName {params.sampleName} --colour_config {params.colour_config} --toggle_label {params.use_alias}"


# adapt this directory in master snake file to prevent recomputing the cohort in each analysis
if not 'ASSEMBLE_NONMALIG_OUT' in globals():
    ASSEMBLE_NONMALIG_OUT = OUTDIR + 'non_malignant_reference/'

# perform assembly of the nonmalignant reference cohort
rule assemble_nonmalignant_cohort:
    input:
        inputDir_hdf5 = config['tools']['assemble_non_malignant_reference']['hdf5_dir'],
        inputDir_classify = config['tools']['assemble_non_malignant_reference']['celltype_dir']
    output:
        outfile = ASSEMBLE_NONMALIG_OUT + 'nonmalignant_reference_cohort.h5'
    params:
        outDir = ASSEMBLE_NONMALIG_OUT,
        non_malignant_types = config['tools']['assemble_non_malignant_reference']['non_malignant_types'],
        lsfoutfile = ASSEMBLE_NONMALIG_OUT + 'nonmalignant_reference_cohort.lsfout.log',
        lsferrfile = ASSEMBLE_NONMALIG_OUT + 'nonmalignant_reference_cohort.lsferr.log',
        scratch = config['tools']['assemble_non_malignant_reference']['scratch'],
        mem = config['tools']['assemble_non_malignant_reference']['mem'],
        time = config['tools']['assemble_non_malignant_reference']['time']
    threads:
        config['tools']['assemble_non_malignant_reference']['threads']
    benchmark:
        ASSEMBLE_NONMALIG_OUT + 'nonmalignant_reference_cohort.benchmark'
    shell:
        "{config[tools][assemble_non_malignant_reference][call]} --hdf5_dir {input.inputDir_hdf5} --celltype_dir {input.inputDir_classify} --non_malignant_celltypes {params.non_malignant_types} --out_dir {params.outDir} "

# plot reference cohort tSNEs
rule plot_tSNEs_nonmalignant_cohort:
    input:
        inputHDF5 = ASSEMBLE_NONMALIG_OUT + 'nonmalignant_reference_cohort.h5'
    output:
        outfile_batch = ASSEMBLE_NONMALIG_OUT + 'nonmalignant_reference_cohort.tSNE_batch.png',
	outfile_ct = ASSEMBLE_NONMALIG_OUT + 'nonmalignant_reference_cohort.tSNE_celltype.png'
    params:
        lsfoutfile = ASSEMBLE_NONMALIG_OUT + 'plot_tSNEs_nonmalignant.lsfout.log',
        lsferrfile = ASSEMBLE_NONMALIG_OUT + 'plot_tSNEs_nonmalignant.lsferr.log',
        scratch = config['tools']['plot_tSNE_nonmalignant']['scratch'],
        mem = config['tools']['plot_tSNE_nonmalignant']['mem'],
        time = config['tools']['plot_tSNE_nonmalignant']['time']
    threads:
        config['tools']['plot_tSNE_nonmalignant']['threads']
    benchmark:
        ASSEMBLE_NONMALIG_OUT + 'plot_tSNEs_nonmalignant.benchmark'
    shell:
        "{config[tools][plot_tSNE_nonmalignant][call]} --hdf5File {input.inputHDF5} --outFile_batch {output.outfile_batch} --outFile_ct {output.outfile_ct} "


if not 'DIFF_EXP_IN' in globals():
    DIFF_EXP_IN = REMOVE_ATYPICAL_OUT
if not 'DIFF_EXP_OUT' in globals():
    DIFF_EXP_OUT = OUTDIR + 'diff_exp/'

# perform the differential expression analysis using a Wilcoxon test
rule diff_exp_genes:
    input:
        sce_in = DIFF_EXP_IN + '{sample}.RDS',
        cell_types = DIFF_EXP_IN + '{sample}.phenograph_celltype_association.txt'
    output:
        #outpath = dynamic(DIFF_EXP_OUT + "{sample}.{clusterid}.DEgenes.tsv"),
        success = DIFF_EXP_OUT + "{sample}.diffExp_success.txt"
    params:
        sampleName = '{sample}',
        malignant = config['inputOutput']['malignant_cell_type'],
        threshold_comparison = config['tools']['diff_exp']['threshold_comparison'],
        fdr_cut = config['tools']['diff_exp']['fdr_cut'],
        fc_cut = config['tools']['diff_exp']['fc_cut'],
        mindiff2second = config['tools']['diff_exp']['mindiff2second'],
        minNumberNonMalignant = config['tools']['diff_exp']['minNumberNonMalignant'],
        lsfoutfile = DIFF_EXP_OUT + '{sample}.diff_exp.lsfout.log',
        lsferrfile = DIFF_EXP_OUT + '{sample}.diff_exp.lsferr.log',
        scratch = config['tools']['diff_exp']['scratch'],
        mem = config['tools']['diff_exp']['mem'],
        time = config['tools']['diff_exp']['time'],
        outpath = DIFF_EXP_OUT,
    threads:
        config['tools']['diff_exp']['threads']
    benchmark:
        DIFF_EXP_OUT + '{sample}.diff_exp.benchmark'
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


if not 'GENE_EXP_IN' in globals():
    GENE_EXP_IN = REMOVE_ATYPICAL_OUT
if not 'GENE_EXP_OUT' in globals():
    GENE_EXP_OUT = OUTDIR + 'gene_exp/'

# give out gene expression values per cluster
rule gene_exp:
    input:
        sce_in = GENE_EXP_IN + '{sample}.RDS',
    output:
        out = GENE_EXP_OUT + "{sample}.gene_expression_per_cluster.tsv"
    params:
        sampleName = '{sample}',
        lsfoutfile = GENE_EXP_OUT + '{sample}.gene_exp.lsfout.log',
        lsferrfile = GENE_EXP_OUT + '{sample}.gene_exp.lsferr.log',
        scratch = config['tools']['gene_exp']['scratch'],
        mem = config['tools']['gene_exp']['mem'],
        time = config['tools']['gene_exp']['time'],
        outpath = GENE_EXP_OUT,
        threshold_sample = config['tools']['gene_exp']['threshold_sample'],
        type_sample = config['tools']['gene_exp']['type_sample'],
        priority_genes = config['resources']['priority_genes'],
    threads:
        config['tools']['gene_exp']['threads']
    benchmark:
        GENE_EXP_OUT + '{sample}.gene_exp.benchmark'
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

if not 'BOXPLOT_IN' in globals():
	BOXPLOT_IN = REMOVE_ATYPICAL_OUT
if not 'BOXPLOT_OUT' in globals():
	BOXPLOT_OUT = PLOTTING_OUT

# This rule creates a box plot comparing cell type fractions across samples
rule generate_cell_type_boxplot:
    input:
        previous_samples = config['resources']['previous_samples'],
        sample_cell_types = BOXPLOT_IN + '{sample}.phenograph_celltype_association.txt'
    output:
        out = BOXPLOT_OUT + '{sample}.boxplot_cell_types_cohort.png'
    params:
        sampleName = '{sample}',
        sampleName_short = config['tools']['cellranger_count']['cellranger_sampleName'],
        outDir = BOXPLOT_OUT,
        lsfoutfile = BOXPLOT_OUT + '{sample}.boxplot_cell_types_cohort.lsfout.log',
        lsferrfile = BOXPLOT_OUT + '{sample}.boxplot_cell_types_cohort.lsferr.log',
        scratch = config['tools']['generate_cell_type_boxplot']['scratch'],
        mem = config['tools']['generate_cell_type_boxplot']['mem'],
        time = config['tools']['generate_cell_type_boxplot']['time'],
    threads:
        config['tools']['generate_cell_type_boxplot']['threads']
    benchmark:
        BOXPLOT_OUT + '{sample}.boxplot_cell_types_cohort.benchmark'
    shell:
        '{config[tools][generate_cell_type_boxplot][call]} --previous_samples {input.previous_samples} --current_sample {input.sample_cell_types} --sampleName {params.sampleName} --sampleName_short {params.sampleName_short} --outDir {params.outDir}'


if not 'INTEGRATION_IN' in globals():
	INTEGRATION_IN = REMOVE_ATYPICAL_OUT
if not 'INTEGRATION_OUT' in globals():
	INTEGRATION_OUT = PLOTTING_OUT

# This rule integrates samples of the cohort and visualizes the integration with UMAPs
rule sample_integration:
    input:
        previous_samples = config['resources']['previous_samples_counts'],
        current_sample = INTEGRATION_IN + '{sample}.RDS'
    output:
        out = INTEGRATION_OUT + '{sample}.sample_integration_highlight_current.png'
    params:
        sampleName = '{sample}',
        outDir = INTEGRATION_OUT,
        lsfoutfile = INTEGRATION_OUT + '{sample}.sample_integration.lsfout.log',
        lsferrfile = INTEGRATION_OUT + '{sample}.sample_integration.lsferr.log',
        scratch = config['tools']['sample_integration']['scratch'],
        mem = config['tools']['sample_integration']['mem'],
        time = config['tools']['sample_integration']['time'],
        sampleName_short = config['tools']['cellranger_count']['cellranger_sampleName'],
	colour_config = config['resources']['colour_config']
    threads:
        config['tools']['sample_integration']['threads']
    benchmark:
        INTEGRATION_OUT + '{sample}.sample_integration.benchmark'
    shell:
        '{config[tools][sample_integration][call]} --cohort_list {input.previous_samples} --sample_data {input.current_sample} --sampleName {params.sampleName} --sampleName_short {params.sampleName_short} --colour_config {params.colour_config} --outdir {params.outDir}'


if not 'PERCENTAGE_IN' in globals():
    PERCENTAGE_IN = REMOVE_ATYPICAL_OUT
if not 'PERCENTAGE_OUT' in globals():
    PERCENTAGE_OUT = OUTDIR + 'clusterpercent/'

# calculate for each cluster the number of cells it countains and the percentage of all cells
rule cellPercentInCluster:
    input:
        clusterCsv = PERCENTAGE_IN + '{sample}.phenograph_celltype_association.txt'
    output:
#    print(drugID)
        out = PERCENTAGE_OUT + '{sample}.clusters_cell_count_percent.txt'
    params:
        lsfoutfile = PERCENTAGE_OUT  + '{sample}.clusterPercent.lsfout.log',
        lsferrfile =  PERCENTAGE_OUT  + '{sample}.clusterPercent.lsferr.log',
        scratch = config['tools']['cellPercentInCluster']['scratch'],
        mem = config['tools']['cellPercentInCluster']['mem'],
        time = config['tools']['cellPercentInCluster']['time'],
        variousParams = config['tools']['cellPercentInCluster']['variousParams']
    threads:
        config['tools']['cellPercentInCluster']['threads']
    benchmark:
        PERCENTAGE_OUT + '{sample}.clusterPercent.benchmark'
    shell:
        '{config[tools][cellPercentInCluster][call]} --inputTable {input.clusterCsv} --outFile {output.out} {params.variousParams}'
