# Rules of the clinical part of the scAmpi pipeline

# Parse csv file with expression levels and extract differentially expressed genes based on filter criteria
rule parse_filter_DE_genes:
    input:
        tsv = 'results/diff_exp_analysis/{sample}/{sample}.{i}.DEgenes.tsv',
    output:
        out = 'results/parse_diff_exp/{sample}.{i}.txt',
    params:
        variousParams = config['tools']['parse_filter_DE_genes']['variousParams'],
        custom_script = workflow.source_path("../scripts/parse_filter_DE_genes.py")
    conda:
        '../envs/parse_filter_DE_genes.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/parse_diff_exp/benchmark/{sample}.{i}.parse.benchmark'
    shell:
        'python {params.custom_script} '
        '{input.tsv} '
        '{output.out} '
        '{params.variousParams}'


# query identified variants at dgidb
rule query_dgidb:
    input:
        infile = 'results/parse_diff_exp/{sample}.{i}.txt'
    output:
        outfile = 'results/query_dgidb/{sample}.{i}.dgidb.txt',
        outfileCompleteTable = 'results/query_dgidb/{sample}.{i}.dgidb.txt.CompleteTable.txt',
        outfileGeneCategory = 'results/query_dgidb/{sample}.{i}.dgidb.txt.GeneCategories.txt'
    params:
        minDatabaseNum = config['tools']['query_dgidb']['minDatabaseNum'],
        colName_genes = config['tools']['query_dgidb']['colName_genes'],
        custom_script = workflow.source_path("../scripts/query_dgidb.R"),
    conda:
        '../envs/query_dgidb.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['medium'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['medium']
    benchmark:
        'databaseQuery/benchmark/{sample}.{i}.dgidbQuery.benchmark'
    shell:
         'Rscript {params.custom_script} '
         '{input.infile} '
         '{output.outfile} '
         '{params.minDatabaseNum} '
         '{params.colName_genes}'


# clinical trials query
rule query_clinical_trials:
    input:
        infile = 'results/query_dgidb/{sample}.{i}.dgidb.txt.CompleteTable.txt',
        downloadSuccess = 'results/clinical_trials/downloadSuccess.txt',
    output:
        outfile = 'results/clinical_trials/{sample}.{i}.dgidb.txt.CompleteTable.ClinicalTrials.allDrugs.txt'
    params:
        cancerType = config['tools']['download_clinical_trials']['cancerType'],
        whiteList = config['tools']['query_clinical_trials']['whiteList'],
        blackList = config['tools']['query_clinical_trials']['blackList'],
        outDirec = 'results/clinical_trials/',
        custom_script = workflow.source_path("../scripts/queryClinicalTrials.py"),
    conda:
        '../envs/query_clinical_trials.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['medium'],
    threads:
        config['computingResources']['threads']['medium']
    benchmark:
        'results/clinical_trials/benchmark/{sample}.{i}.clinicalTrialsQuery.benchmark'
    shell:
        'python {params.custom_script} '
        '{input.infile} '
        '{output.outfile} '
        '{params.outDirec}/{params.cancerType}_clinicalTrials/ '
        '"{params.whiteList}" '
        '"{params.blackList}" '


# download the clinical trials necessary for the query
rule download_clinical_trials:
    output:
        outfile = 'results/clinical_trials/downloadSuccess.txt'
    params:
        cancerType = config['tools']['download_clinical_trials']['cancerType'],
        outDirec = 'results/clinical_trials/',
#    conda:
#        '../envs/query_dgidb.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['medium'],
    threads:
        config['computingResources']['threads']['medium']
    benchmark:
        'results/clinical_trials/benchmark/downloadClinicalTrials.benchmark'
    shell:
        ('wget "https://clinicaltrials.gov/search?term={params.cancerType}&studyxml=true" '
        '-O {params.outDirec}/{params.cancerType}_clinicalTrials.zip ; '
        'unzip {params.outDirec}/{params.cancerType}_clinicalTrials.zip '
        '-d {params.outDirec}/{params.cancerType}_clinicalTrials ; '
        'touch {params.outDirec}/downloadSuccess.txt' )


# Combine different database queries, annotate input table with clinical information
rule annotate_DE_clinical_info:
    input:
        infile = 'results/parse_diff_exp/{sample}.{i}.txt',
        pathwayDB = config['resources']['pathwayDB'],
        inDGIDB = 'results/query_dgidb/{sample}.{i}.dgidb.txt.GeneCategories.txt',
        inClinicalTrials = 'results/clinical_trials/{sample}.{i}.dgidb.txt.CompleteTable.ClinicalTrials.allDrugs.txt'
    output:
        outTable = 'results/clinical_annotation/{sample}.{i}.clinicalAnnotation.txt',
        outTable_dgdidbIndependent = 'results/clinical_annotation/{sample}.{i}.clinicalAnnotation.txt_dgidbIndependent.txt'
    params:
        variousParams = config['tools']['annotate_DE_clinical_info']['variousParams'],
        custom_script = workflow.source_path("../scripts/annotate_DE_clinical_info.py"),
    conda:
        '../envs/annotate_DE_clinical_info.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/clinical_annotation/benchmark/{sample}.{i}.annotate_DE_clinical_info.benchmark'
    shell:
        'python {params.custom_script} '
        '--inputTable {input.infile} '
        '--outFile {output.outTable} '
        '--pathwayDB {input.pathwayDB} '
        '--dgidb_categ {input.inDGIDB} '
        '--clinTrials {input.inClinicalTrials} '
        '{params.variousParams} '


# query identified expression in civic
rule query_civic:
    input:
        infile = 'results/clinical_annotation/{sample}.{i}.clinicalAnnotation.txt'
    output:
        outfile = 'results/query_civic/{sample}.{i}.clinicalAnnotation.civic.txt'
    params:
        cancerType = config['tools']['query_civic']['cancerType'],
        blackList = config['tools']['query_civic']['blackList'],
        highLevel = config['tools']['query_civic']['highLevel'],
        colName_gene = config['tools']['query_civic']['colName_gene'],
        colName_logFC = config['tools']['query_civic']['colName_logFC'],
        strictExpression = config['tools']['query_civic']['strictExpression'],
        custom_script = workflow.source_path("../scripts/query_civic_expr.py"),
    conda:
        '../envs/query_civic.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/query_civic/benchmark/{sample}.{i}.query_civic.benchmark'
    shell:
        'python {params.custom_script} '
        '--inputTable {input.infile} '
        '--outFile {output.outfile} '
        '--cancerTypeList "{params.cancerType}" '
        '--blackList "{params.blackList}" '
        '--highLevelList "{params.highLevel}" '
        '--colName_gene {params.colName_gene} '
        '--colName_logFC {params.colName_logFC} '
        '--strictExpression {params.strictExpression}'


# gene set enrichment on the DE results
rule gene_set_enrichment:
    input:
        infile = 'results/diff_exp_analysis/{sample}/{sample}.{i}.DEgenes.tsv'
    output:
        outfile = 'results/gene_set_enrichment/{sample}.{i}.enrichedGeneSets.txt'
    params:
        geneSetDB = config['resources']['genesets'],
        variousParams = config['tools']['gene_set_enrichment']['variousParams'],
        custom_script = workflow.source_path("../scripts/gene_set_enrichment.R"),
    conda:
        '../envs/gene_set_enrichment.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/gene_set_enrichment/benchmark/{sample}.{i}.gene_set_enrichment.benchmark'
    shell:
        'Rscript {params.custom_script} '
        '{input.infile} '
        '{output.outfile} '
        '{params.geneSetDB} '
        '{params.variousParams}'


# gene set enrichment on the DE results malignant vs. malignant
rule gene_set_enrichment_mal_vs_mal:
    input:
        infile = 'results/diff_exp_analysis/{sample}/vs_other_malignant/{sample}.DEmalignant.{i}.DEgenes.tsv'
    output:
        outfile = 'results/gene_set_enrichment/vs_other_malignant/{sample}.DEmalignant.{i}.enrichedGeneSets.txt'
    params:
        geneSetDB = config['resources']['genesets'],
        variousParams = config['tools']['gene_set_enrichment']['variousParams'],
        custom_script = workflow.source_path("../scripts/gene_set_enrichment.R"),
    conda:
        '../envs/gene_set_enrichment.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/gene_set_enrichment/vs_other_malignant/benchmark/{sample}.DEmalignant.{i}.gene_set_enrichment.benchmark'
    shell:
        'Rscript {params.custom_script} '
        '{input.infile} '
        '{output.outfile} '
        '{params.geneSetDB} '
        '{params.variousParams}'


def getGeneSetHeatmapFiles(wildcards):
#    if "DEmalignant" in wildcards.sample:
#        if len(CLUSTER_IDS_MALIGNANT) == 0:
#            return [GENESETANALYSIS_OUT]
#        return expand(GENESETANALYSIS_OUT + wildcards.sample + '.{clusterid}.enrichedGeneSets.txt', clusterid = CLUSTER_IDS_MALIGNANT)
    checkpoint_output = checkpoints.diff_exp_analysis.get(**wildcards).output[0]
    return expand('results/gene_set_enrichment/{sample}.{i}.enrichedGeneSets.txt',
#        return expand("results/{result_outdir}/{sample}.{i}.{suffix}.txt",
            i = glob_wildcards(os.path.join(checkpoint_output, "{sample,[^/]+}.{i,[^/]+}.DEgenes.tsv")).i,
            sample = wildcards.sample)

def getGeneSetHeatmapFilesMalignant(wildcards):
#    if "DEmalignant" in wildcards.sample:
#        if len(CLUSTER_IDS_MALIGNANT) == 0:
#            return [GENESETANALYSIS_OUT]
#        return expand(GENESETANALYSIS_OUT + wildcards.sample + '.{clusterid}.enrichedGeneSets.txt', clusterid = CLUSTER_IDS_MALIGNANT)
    checkpoint_output = checkpoints.diff_exp_analysis.get(**wildcards).output[0]
    return expand('results/gene_set_enrichment/vs_other_malignant/{sample}.DEmalignant.{i}.enrichedGeneSets.txt',
#        return expand("results/{result_outdir}/{sample}.{i}.{suffix}.txt",
            i = glob_wildcards(os.path.join(checkpoint_output, "vs_other_malignant/{sample,[^/]+}.DEmalignant.{i,[^/]+}.DEgenes.tsv")).i,
            sample = wildcards.sample)


## plot heat map for gene set enrichment
## NOTE: empty inputs cause snakemake to crash, thus we implemented the workaround
## that returns GENESETANALYSIS_OUT in case of an empty list of enrichment analyses comparing tumor sub-clones
rule plot_gene_set_enrichment:
    input:
        inDir = getGeneSetHeatmapFiles
        #inDir = expand(GENESETANALYSIS_OUT + '{{sample}}.{clusterid}.enrichedGeneSets.txt', clusterid = CLUSTER_IDS)
    output:
        outfile = 'results/gene_set_enrichment/{sample}.heatmap_enrichment.png'
    params:
        comparison_direc = 'results/gene_set_enrichment/',
        custom_script = workflow.source_path("../scripts/plot_genesets_heatmap.R"),
    conda:
        '../envs/plot_gene_set_enrichment.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/gene_set_enrichment/benchmark/{sample}.plot_gene_set_enrichment.benchmark'
    shell:
        'if [ "{input.inDir}" != "{params.comparison_direc}" ] ; '
        'then echo "test1" ; '
        'Rscript {params.custom_script} {output.outfile} {input.inDir} ; '
        'else touch {output.outfile} ; '
        'fi'


## plot heat map for gene set enrichment
## NOTE: empty inputs cause snakemake to crash, thus we implemented the workaround
## that returns GENESETANALYSIS_OUT in case of an empty list of enrichment analyses comparing tumor sub-clones
rule plot_gene_set_enrichment_mal_vs_mal:
    input:
        inDir = getGeneSetHeatmapFilesMalignant
        #inDir = expand(GENESETANALYSIS_OUT + '{{sample}}.{clusterid}.enrichedGeneSets.txt', clusterid = CLUSTER_IDS)
    output:
        outfile = 'results/gene_set_enrichment/vs_other_malignant/{sample}.DEmalignant.heatmap_enrichment.png'
    params:
        comparison_direc = 'results/gene_set_enrichment/vs_other_malignant/',
        custom_script = workflow.source_path("../scripts/plot_genesets_heatmap.R"),
    conda:
        '../envs/plot_gene_set_enrichment.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/gene_set_enrichment/vs_other_malignant/benchmark/{sample}.DEmalignant.plot_gene_set_enrichment.benchmark'
    shell:
        'if [ "{input.inDir}" != "{params.comparison_direc}" ] ; '
        'then echo "test1" ; '
        'Rscript {params.custom_script} {output.outfile} {input.inDir} ; '
        'else touch {output.outfile} ; '
        'fi'


# parse the *.dgidb.txt.CompleteTable.ClinicalTrials.txt files of all clusters
# and generate a table that shows for each Drug the clusters that can be targeted by this drug
# and the weight of the drug for calculating minimum set cover
def get_parse_for_minSetCover_input(wildcards):
    checkpoint_output = checkpoints.diff_exp_analysis.get(**wildcards).output[0]
    return expand('results/clinical_trials/{sample}.{i}.dgidb.txt.CompleteTable.ClinicalTrials.{{type}}.txt',
            i = glob_wildcards(os.path.join(checkpoint_output, "{sample,[^/]+}.{i,[^/]+}.DEgenes.tsv")).i,
            sample = wildcards.sample)

rule parse_for_minSetCover:
    input:
#        infiles = expand(PARSETRIALSTABLE_IN + '{{sample}}.{i}.dgidb.txt.CompleteTable.ClinicalTrials.{{type}}.txt', clusterid = CLUSTER_IDS),
        infiles = get_parse_for_minSetCover_input,
        drugList = config['resources']['drugList']
    output:
        out = 'results/drug_combination/{sample}.drugToCluster.{type}.txt'
    params:
        colName_clinTrial = config['tools']['parse_for_minSetCover']['colName_clinTrial'],
        colName_DGIDB_score = config['tools']['parse_for_minSetCover']['colName_DGIDB_score'],
        custom_script = workflow.source_path("../scripts/parse_for_minSetCover.py"),
    conda:
        '../envs/parse_for_minSetCover.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/drug_combination/benchmark/{sample}.{type}.parse_for_minSetCover.benchmark'
    shell:
        'python {params.custom_script} '
        '--inFiles {input.infiles} '
        '--outFile {output.out} '
        '--colName_clinTrial {params.colName_clinTrial} '
        '--colName_DGIDB_Score {params.colName_DGIDB_score} '
        '--drug_list {input.drugList}'


# calculate for each cluster the number of cells it countains and the percentage of all cells
rule cell_percent_in_cluster:
    input:
        clusterCsv = 'results/atypical_removed/{sample}.atypical_removed.phenograph_celltype_association.txt'
    output:
        out = 'results/clustering/{sample}.clusters_cell_count_percent.txt'
    params:
        variousParams = config['tools']['cell_percent_in_cluster']['variousParams'],
        custom_script = workflow.source_path("../scripts/cell_percent_in_cluster.py"),
    conda:
        '../envs/cell_percent_in_cluster.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['medium'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['medium']
    benchmark:
        'results/clustering/benchmark/{sample}.clusterPercent.benchmark'
    shell:
        'python {params.custom_script} '
        '--inputTable {input.clusterCsv} '
        '--outFile {output.out} '
        '{params.variousParams}'



# deduce minimum set cover to find drug (combination) that targets all clusters via interaction with DE gene of this cluster
rule find_minSetCover:
    input:
        infile = 'results/drug_combination/{sample}.drugToCluster.{type}.txt',
        percTable = 'results/clustering/{sample}.clusters_cell_count_percent.txt'
    output:
        out = 'results/drug_combination/{sample}.drugCombination.{type}.txt'
    params:
        variousParams = config['tools']['find_minSetCover']['variousParams'],
        custom_script = workflow.source_path("../scripts/find_minSetCover.py"),
    conda:
        '../envs/find_minSetCover.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/drug_combination/benchmark/{sample}.{type}.find_minSetCover.benchmark'
    shell:
        'python {params.custom_script} '
        '--input {input.infile} '
        '--outFile {output.out} '
        '--percentageTable {input.percTable} '
        '{params.variousParams}'


# filter the drug-gene interaction results for drugs that are included in TP melanoma clinical list of drugs
rule filter_drugs:
    input:
        infile = 'results/clinical_trials/{sample}.{i}.dgidb.txt.CompleteTable.ClinicalTrials.allDrugs.txt',
        drugList = config['resources']['drugList']
    output:
        out = 'results/clinical_trials/{sample}.{i}.dgidb.txt.CompleteTable.ClinicalTrials.filteredDrugs.txt'
    params:
        custom_script = workflow.source_path("../scripts/filter_drugs.py"),
    conda:
        '../envs/filter_drugs.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/clinical_trials/benchmark/{sample}.{i}.filter_drugs.benchmark'
    shell:
        'python {params.custom_script} '
        '--inFile {input.infile} '
        '--outFile {output.out} '
        '--drugList {input.drugList}'


# this rule is to preprocess the output of the rule parse_for_minSetCover of type drug,clusters,weight (tab separated)
# to have the necessary input format for the UpSetR package.
# That is drug,cluster1,cluster2..clustern (tab separated) with 1 or 0 indicating if a drug targets a DE of a cluster
rule preprocess_upsetr_plot:
    input:
        infile = 'results/drug_combination/{sample}.drugToCluster.{type}.txt'
    output:
        out = 'results/upsetr_plot/{sample}.drugToCluster.{type}.processedForUpSetR.txt'
    params:
        custom_script = workflow.source_path("../scripts/preprocess_upsetr_plot.py"),
    conda:
        '../envs/preprocess_upsetr_plot.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/upsetr_plot/benchmark/{sample}.{type}.upsetr_plot.benchmark'
    shell:
        'python {params.custom_script} '
        '--inFile {input.infile} '
        '--outFile {output.out}'


# this rule generates a plot with UpSetR (comparably to venn diagramm)
# The plot displays intersections of the drug sets that target DE genes in the clusters
rule plot_upsetr:
    input:
        infile = 'results/upsetr_plot/{sample}.drugToCluster.{type}.processedForUpSetR.txt'
    output:
        out = 'results/upsetr_plot/{sample}.drugToCluster.{type}.vennplot.png'
    params:
        variousParams = config['tools']['plot_upsetr']['variousParams'],
        custom_script = workflow.source_path("../scripts/plot_upsetr.R"),
    conda:
        '../envs/plot_upsetr.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/upsetr_plot/benchmark/{sample}.{type}.plot_upsetr.benchmark'
    shell:
        'Rscript {params.custom_script} '
        '--inFile {input.infile} '
        '--outFile {output.out} '
        '{params.variousParams}'


# This rule generates a table with a full list of all clinically relevant drugs and the information
# for each cluster if a drug-gene interaction was found in dgidb between the respective drug and
# any differentially expressed gene of the cluster.
rule get_full_druglist_to_subclones:
    input:
        infile = 'results/drug_combination/{sample}.drugToCluster.allDrugs.txt',
    output:
        out = 'results/drug_combination/{sample}.full_druglist_to_subclones.txt'
    params:
        drugList = config['resources']['drugList'],
        custom_script = workflow.source_path("../scripts/get_full_druglist_to_subclones_assignm.py"),
    conda:
        '../envs/get_full_druglist_to_subclones.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['low'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['low']
    benchmark:
        'results/drug_combination/benchmark/{sample}.full_druglist_to_subclones.benchmark'
    shell:
        'python {params.custom_script} '
        '--in_drugToCluster {input.infile} '
        '--in_drugList {params.drugList} '
        '--outFile {output.out} '


# check whether all civic queries are finished
def getCivicQueryResults(wildcards):
    checkpoint_output = checkpoints.diff_exp_analysis.get(**wildcards).output[0]
#    return expand(PLOT_DRUGS_IN + wildcards.sample + '.{i}.clinicalAnnotation.civic.txt',
    return expand('results/query_civic/{sample}.{i}.clinicalAnnotation.civic.txt',
        sample = wildcards.sample,
        i = glob_wildcards(os.path.join(checkpoint_output, "{sample,[^/]+}.{i,[^/]+}.DEgenes.tsv")).i)


# this rule generates a UMAP plot that shows drug prediction on the tumor clones
rule plot_drug_prediction:
    input:
        rdsFile = 'results/atypical_removed/{sample}.atypical_removed.RDS',
        inFiles = getCivicQueryResults,
        drugList = config['resources']['drugList'],
        drugCombis = config['resources']['drugCombinations'],
        civicDict = config['resources']['civicDict']
    output:
        out = 'results/plot_drug_prediction/{sample}.drug_prediction_umap.png'
    params:
        sampleName = '{sample}',
        inputDir = 'results/query_civic/',
        outputDirec = 'results/plot_drug_prediction/',
        variousParams = config['tools']['plot_drug_prediction']['variousParams'],
        custom_script = workflow.source_path("../scripts/show_drugPrediction_on_clones.R"),
    conda:
        '../envs/plot_drug_prediction.yaml'
    resources:
        mem_mb = config['computingResources']['mem']['medium'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['medium']
    benchmark:
        'results/plot_drug_prediction/benchmark/{sample}.plot_drug_prediction.benchmark'
    shell:
        'Rscript {params.custom_script} '
        '--SCE {input.rdsFile} '
        '--drugPredDir {params.inputDir} '
        '--outputDirec {params.outputDirec} '
        '--drugList {input.drugList} '
        '--combiList {input.drugCombis} '
        '--civicDict {input.civicDict} '
        '--sampleName {params.sampleName} '
        '{params.variousParams} '



# define input for aggregate rule
# restrain the wildcards (the global restraint alone does not work)
def aggregate_input(suffix, result_outdir):
    """Define a checkpoint compatible function that generates filenames."""
    def tmp(wildcards):
        checkpoint_output = checkpoints.diff_exp_analysis.get(**wildcards).output[0]
        return expand("results/{result_outdir}/{sample}.{i}.{suffix}.txt",
            result_outdir = result_outdir,
            suffix = suffix,
            sample = wildcards.sample,
            # restrain wildcards to not consider subdirectories
            # do this here because `glob_wildcards` does not
            # respect the global wildcard_constraints directive
            i = glob_wildcards(os.path.join(checkpoint_output, "{sample,[^/]+}.{i,[^/]+}.DEgenes.tsv")).i)
    return tmp


# define what rules after the checkpoint differential expression should be performed.
# the aggregate rule is triggered by the final scampi rule
rule aggregate:
    input:
        aggregate_input(config['tools']['aggregate']['suffix'],
                        config['tools']['aggregate']['result_outdir'])
    output:
        'results/{outdir}/{sample}.aggregated.txt'
    resources:
        mem_mb = config['computingResources']['mem']['medium'],
        time_min = config['computingResources']['time']['low'],
    threads:
        config['computingResources']['threads']['medium']
    shell:
        'touch {output}'

