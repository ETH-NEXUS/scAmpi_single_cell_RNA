# Rules of the clinical part of the scAmpi pipeline

# Parse csv file with expression levels and extract differentially expressed genes based on filter criteria
rule parse_filter_DE_genes:
    input:
        tsv = 'results/diff_exp_analysis/{sample}/{sample}.{i}.DEgenes.tsv',
    output:
        out = 'results/parse_diff_exp/{sample}.{i}.txt',
    params:
        variousParams = config['tools']['parse_filter_DE_genes']['variousParams']
    conda:
        '../envs/parse_filter_DE_genes.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/parse_diff_exp/benchmark/{sample}.{i}.parse.benchmark'
    shell:
        'python workflow/scripts/parse_filter_DE_genes.py '
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
        colName_genes = config['tools']['query_dgidb']['colName_genes']
    conda:
        '../envs/query_dgidb.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'databaseQuery/benchmark/{sample}.{i}.dgidbQuery.benchmark'
    shell:
         'Rscript workflow/scripts/query_dgidb.R '
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
        outDirec = 'results/clinical_trials/'
    conda:
        '../envs/query_clinical_trials.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/clinical_trials/benchmark/{sample}.{i}.clinicalTrialsQuery.benchmark'
    shell:
        'python workflow/scripts/queryClinicalTrials.py '
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
        outDirec = 'results/clinical_trials/'
#    conda:
#        '../envs/query_dgidb.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
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
        variousParams = config['tools']['annotate_DE_clinical_info']['variousParams']
    conda:
        '../envs/annotate_DE_clinical_info.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/clinical_annotation/benchmark/{sample}.{i}.annotate_DE_clinical_info.benchmark'
    shell:
        'python workflow/scripts/annotate_DE_clinical_info.py '
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
        strictExpression = config['tools']['query_civic']['strictExpression']
    conda:
        '../envs/query_civic.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/query_civic/benchmark/{sample}.{i}.query_civic.benchmark'
    shell:
        'python workflow/scripts/query_civic_expr.py '
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
        variousParams = config['tools']['gene_set_enrichment']['variousParams']
    conda:
        '../envs/gene_set_enrichment.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/gene_set_enrichment/benchmark/{sample}.{i}.gene_set_enrichment.benchmark'
    shell:
        'Rscript workflow/scripts/gene_set_enrichment.R '
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
        variousParams = config['tools']['gene_set_enrichment']['variousParams']
    conda:
        '../envs/gene_set_enrichment.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/gene_set_enrichment/vs_other_malignant/benchmark/{sample}.DEmalignant.{i}.gene_set_enrichment.benchmark'
    shell:
        'Rscript workflow/scripts/gene_set_enrichment.R '
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
            i = glob_wildcards(os.path.join(checkpoint_output, "{sample,[^/]+}.{i,[^/]+}.DEgenes.tsv")).i,
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
        comparison_direc = 'results/gene_set_enrichment/'
    conda:
        '../envs/plot_gene_set_enrichment.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/gene_set_enrichment/benchmark/{sample}.plot_gene_set_enrichment.benchmark'
    shell:
        'if [ "{input.inDir}" != "{params.comparison_direc}" ] ; '
        'then echo "test1" ; '
        'Rscript workflow/scripts/plot_genesets_heatmap.R {output.outfile} {input.inDir} ; '
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
        comparison_direc = 'results/gene_set_enrichment/vs_other_malignant/'
    conda:
        '../envs/plot_gene_set_enrichment.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/gene_set_enrichment/vs_other_malignant/benchmark/{sample}.DEmalignant.plot_gene_set_enrichment.benchmark'
    shell:
        'if [ "{input.inDir}" != "{params.comparison_direc}" ] ; '
        'then echo "test1" ; '
        'Rscript workflow/scripts/plot_genesets_heatmap.R {output.outfile} {input.inDir} ; '
        'else touch {output.outfile} ; '
        'fi'


#if not 'DRUGCOMBINATION' in globals():
#    DRUGCOMBINATION = OUTDIR + 'drugCombination/'
#
#if not 'PARSETRIALSTABLE_IN' in globals():
#    PARSETRIALSTABLE_IN = CLINICALTRIALS_OUT
#if not 'PARSETRIALSTABLE_OUT' in globals():
#    PARSETRIALSTABLE_OUT = DRUGCOMBINATION
#
## parse the *.dgidb.txt.CompleteTable.ClinicalTrials.txt files of all clusters
## and generate a table that shows for each Drug the clusters that can be targeted by this drug
## and the weight of the drug for calculating minimum set cover
#
#rule parseDgidbTrialsTable_for_minSetCover:
#    input:
#        infiles = expand(PARSETRIALSTABLE_IN + '{{sample}}.{clusterid}.dgidb.txt.CompleteTable.ClinicalTrials.{{type}}.txt', clusterid = CLUSTER_IDS),
#        drugList = config['resources']['drugList']
#    output:
#        out = PARSETRIALSTABLE_OUT + '{sample}.drugToCluster.{type}.txt'
#    params:
#        colName_clinTrial = config['tools']['parseDgidbTrialsTable_for_minSetCover']['colName_clinTrial'],
#        colName_DGIDB_score = config['tools']['parseDgidbTrialsTable_for_minSetCover']['colName_DGIDB_score'],
#        lsfoutfile = PARSETRIALSTABLE_OUT + '{sample}.drugToCluster.{type}.lsfout.log',
#        lsferrfile = PARSETRIALSTABLE_OUT + '{sample}.drugToCluster.{type}.lsferr.log',
#        scratch = config['tools']['parseDgidbTrialsTable_for_minSetCover']['scratch'],
#        mem = config['tools']['parseDgidbTrialsTable_for_minSetCover']['mem'],
#        time = config['tools']['parseDgidbTrialsTable_for_minSetCover']['time']
#    threads:
#        config['tools']['parseDgidbTrialsTable_for_minSetCover']['threads']
#    benchmark:
#        PARSETRIALSTABLE_OUT + '{sample}.drugToCluster.{type}.benchmark'
#    shell:
#        '{config[tools][parseDgidbTrialsTable_for_minSetCover][call]} --inFiles {input.infiles} --outFile {output.out} --colName_clinTrial {params.colName_clinTrial} --colName_DGIDB_Score {params.colName_DGIDB_score} --drug_list {input.drugList}'
#
#if not 'MINSETCOVER_IN' in globals():
#    MINSETCOVER_IN = PARSETRIALSTABLE_OUT
#if not 'MINSETCOVER_OUT' in globals():
#    MINSETCOVER_OUT = DRUGCOMBINATION
#
## deduce minimum set cover to find drug (combination) that targets all clusters via interaction with DE gene of this cluster
#
#rule findminSetCover:
#    input:
#        infile = MINSETCOVER_IN + '{sample}.drugToCluster.{type}.txt',
#        percTable = PERCENTAGE_OUT + '{sample}.clusters_cell_count_percent.txt'
#    output:
#        out = MINSETCOVER_OUT + '{sample}.drugCombination.{type}.txt'
#    params:
#        lsfoutfile = MINSETCOVER_OUT + '{sample}.drugCombination.{type}.lsfout.log',
#        lsferrfile = MINSETCOVER_OUT + '{sample}.drugCombination.{type}.lsferr.log',
#        scratch = config['tools']['findminSetCover']['scratch'],
#        mem = config['tools']['findminSetCover']['mem'],
#        time = config['tools']['findminSetCover']['time'],
#        variousParams = config['tools']['findminSetCover']['variousParams']
#    threads:
#        config['tools']['findminSetCover']['threads']
#    benchmark:
#        MINSETCOVER_OUT + '{sample}.drugCombination.{type}.benchmark'
#    shell:
#        '{config[tools][findminSetCover][call]} --input {input.infile} --outFile {output.out} --percentageTable {input.percTable} {params.variousParams}'
#
#
#if not 'FILTERDRUGS_IN' in globals():
#    FILTERDRUGS_IN = CLINICALTRIALS_OUT
#if not 'FILTERDRUGS_OUT' in globals():
#    FILTERDRUGS_OUT = CLINICALTRIALS_OUT
## filter the drug-gene interaction results for drugs that are included in TP melanoma clinical list of drugs
#rule filterDrugs:
#    input:
#        infile = FILTERDRUGS_IN + '{sample}.{clusterid}.dgidb.txt.CompleteTable.ClinicalTrials.allDrugs.txt',
#        drugList = config['resources']['drugList']
#    output:
#        out = FILTERDRUGS_OUT + '{sample}.{clusterid}.dgidb.txt.CompleteTable.ClinicalTrials.filteredDrugs.txt'
#    params:
#        lsfoutfile = FILTERDRUGS_OUT + '{sample}.{clusterid}.filterDrugs.lsfout.log',
#        lsferrfile = FILTERDRUGS_OUT + '{sample}.{clusterid}.filterDrugs.lsferr.log',
#        scratch = config['tools']['filterDrugs']['scratch'],
#        mem = config['tools']['filterDrugs']['mem'],
#        time = config['tools']['filterDrugs']['time'],
#    threads:
#        config['tools']['filterDrugs']['threads']
#    benchmark:
#        FILTERDRUGS_OUT + '{sample}.{clusterid}.filterDrugs.benchmark'
#    shell:
#        '{config[tools][filterDrugs][call]} --inFile {input.infile} --outFile {output.out} --drugList {input.drugList}'
#
#
#if not 'PREPROCESSUPSETR_IN' in globals():
#    PREPROCESSUPSETR_IN = DRUGCOMBINATION
#if not 'PREPROCESSUPSETR_OUT' in globals():
#    PREPROCESSUPSETR_OUT = DRUGCOMBINATION
#
## this rule is to preprocess the output of the rule parseDgidbTrialsTable_for_minSetCover of type drug,clusters,weight (tab separated)
## to have the necessary input format for the UpSetR package, which is drug,cluster1,cluster2..clustern (tab separated) with 1 or 0 indicating if a drug targets a DE of a cluster
#rule preprocessForUpSetR_venn:
#    input:
#        infile = PREPROCESSUPSETR_IN + '{sample}.drugToCluster.{type}.txt'
#    output:
#        out = PREPROCESSUPSETR_OUT + '{sample}.drugToCluster.{type}.processedForUpSetR.txt'
#    params:
#        lsfoutfile = PREPROCESSUPSETR_OUT + '{sample}.drugToCluster.{type}.processedForUpSetR.lsfout.log',
#        lsferrfile = PREPROCESSUPSETR_OUT + '{sample}.drugToCluster.{type}.processedForUpSetR.lsferr.log',
#        scratch = config['tools']['preprocessForUpSetR']['scratch'],
#        mem = config['tools']['preprocessForUpSetR']['mem'],
#        time = config['tools']['preprocessForUpSetR']['time']
#    threads:
#        config['tools']['preprocessForUpSetR']['threads']
#    benchmark:
#        PREPROCESSUPSETR_OUT + '{sample}.drugToCluster.{type}.processedForUpSetR.benchmark'
#    shell:
#        '{config[tools][preprocessForUpSetR][call]} --inFile {input.infile} --outFile {output.out}'
#
#
#if not 'PLOTUPSETR_IN' in globals():
#    PLOTUPSETR_IN = DRUGCOMBINATION
#if not 'PLOTUPSETR_OUT' in globals():
#    PLOTUPSETR_OUT = DRUGCOMBINATION
#
## this rule generates a plot with UpSetR (comparably to venn diagramm)
## The plot displays intersections of the drug sets that target DE genes in the clusters
#rule plotUpSetR_venn:
#    input:
#        infile = PLOTUPSETR_IN + '{sample}.drugToCluster.{type}.processedForUpSetR.txt'
#    output:
#        out = PLOTUPSETR_OUT + '{sample}.drugToCluster.{type}.vennplot.png'
#    params:
#        lsfoutfile = PLOTUPSETR_OUT + '{sample}.drugToCluster.{type}.vennplot.lsfout.log',
#        lsferrfile = PLOTUPSETR_OUT + '{sample}.drugToCluster.{type}.vennplot.lsferr.log',
#        scratch = config['tools']['plotUpSetR']['scratch'],
#        mem = config['tools']['plotUpSetR']['mem'],
#        time = config['tools']['plotUpSetR']['time'],
#        variousParams = config['tools']['plotUpSetR']['variousParams']
#    threads:
#        config['tools']['plotUpSetR']['threads']
#    benchmark:
#        PLOTUPSETR_OUT + '{sample}.drugToCluster.{type}.vennplot.benchmark'
#    shell:
#        '{config[tools][plotUpSetR][call]} --inFile {input.infile} --outFile {output.out} {params.variousParams}'
#
#
#if not 'FULL_DRUGLIST_TO_SUBCLONES_IN' in globals():
#    FULL_DRUGLIST_TO_SUBCLONES_IN = PARSETRIALSTABLE_OUT
#if not 'FULL_DRUGLIST_TO_SUBCLONES_OUT' in globals():
#    FULL_DRUGLIST_TO_SUBCLONES_OUT = PARSETRIALSTABLE_OUT
#
## This rule generates a table with a full list of all clinically relevant drugs and the information
## for each cluster if a drug gene interaction was found in dgidb between the respective drug and 
## any differentially expressed gene of the cluster.
#
#rule get_full_druglist_to_subclones:
#    input:
#        infile = FULL_DRUGLIST_TO_SUBCLONES_IN + '{sample}.drugToCluster.allDrugs.txt',
#    output:
#        out = FULL_DRUGLIST_TO_SUBCLONES_OUT + '{sample}.full_druglist_to_subclones.txt'
#    params:
#        lsfoutfile = FULL_DRUGLIST_TO_SUBCLONES_OUT + '{sample}.full_druglist_to_subclones.lsfout.log',
#        lsferrfile = FULL_DRUGLIST_TO_SUBCLONES_OUT + '{sample}.full_druglist_to_subclones.lsferr.log',
#        scratch = config['tools']['get_full_druglist_to_subclones']['scratch'],
#        mem = config['tools']['get_full_druglist_to_subclones']['mem'],
#        time = config['tools']['get_full_druglist_to_subclones']['time'],
#        drugList = config['resources']['drugList']
#    threads:
#        config['tools']['get_full_druglist_to_subclones']['threads']
#    benchmark:
#        FULL_DRUGLIST_TO_SUBCLONES_OUT + '{sample}.full_druglist_to_subclones.benchmark'
#    shell:
#        '{config[tools][get_full_druglist_to_subclones][call]} --in_drugToCluster {input.infile} --in_drugList {params.drugList} --outFile {output.out} '
#
#
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
        variousParams = config['tools']['plot_drug_prediction']['variousParams']
    conda:
        '../envs/plot_drug_prediction.yaml'
    resources:
        mem_mb = config['computingResources']['mediumRequirements']['mem'],
        time_min = config['computingResources']['mediumRequirements']['time'],
    threads:
        config['computingResources']['mediumRequirements']['threads']
    benchmark:
        'results/plot_drug_prediction/benchmark/{sample}.plot_drug_prediction.benchmark'
    shell:
        'Rscript workflow/scripts/show_drugPrediction_on_clones.R '
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
    shell:
        'touch {output}'


###   TODO

# define input for aggregate rule
# restrain the wildcards (the global restraint alone does not work)
#def aggregate_malignant_vs_malignant(suffix, result_outdir):
#    """Define a checkpoint compatible function that generates filenames."""
#    def tmp(wildcards):
#        checkpoint_output = checkpoints.diff_exp_analysis.get(**wildcards).output[0]
#        return expand("results/{result_outdir}/{sample}.{i}.{suffix}.txt",
#            result_outdir = result_outdir,
#            suffix = suffix,
#            sample = wildcards.sample,
#            # restrain wildcards to not consider subdirectories
#            # do this here because `glob_wildcards` does not
#            # respect the global wildcard_constraints directive
#            i = glob_wildcards(os.path.join(checkpoint_output, "{sample,[^/]+}.{i,[^/]+}.DEgenes.tsv")).i)
#    return tmp

# define aggregated malignant vs. malignant gene set enrichment
# the aggregate rule is triggered by the final scampi rule
#rule aggregate:
#    input:
#        aggregate_input(config['tools']['aggregate']['suffix'],
#        config['tools']['aggregate']['result_outdir'])
#    output:
#        'results/{outdir}/{sample}.aggregated.txt'
#    shell:
#        'touch {output}'
