
if not 'PROCESSDE_IN' in globals():
    PROCESSDE_IN = DIFF_EXP_OUT
if not 'PROCESSDE_OUT' in globals():
    PROCESSDE_OUT = OUTDIR + 'genes_DE/'

# Parse csv file with expression levels and extract differentially expressed genes based on filter criteria
rule parseAndFilter_DEgenes:
    input:
        tsv = PROCESSDE_IN + '{sample}.{clusterid}.DEgenes.tsv',
    output:
        out = PROCESSDE_OUT + '{sample}.{clusterid}.txt'
    params:
        lsfoutfile = PROCESSDE_OUT + '{sample}.{clusterid}.parse.lsfout.log',
        lsferrfile = PROCESSDE_OUT + '{sample}.{clusterid}.parse.lsferr.log',
        scratch = config['tools']['parseAndFilter_DEgenes']['scratch'],
        mem = config['tools']['parseAndFilter_DEgenes']['mem'],
        time = config['tools']['parseAndFilter_DEgenes']['time'],
        variousParams = config['tools']['parseAndFilter_DEgenes']['variousParams']
    threads:
        config['tools']['parseAndFilter_DEgenes']['threads']
    benchmark:
        PROCESSDE_OUT + '{sample}.{clusterid}.parse.benchmark'
    shell:
        '{config[tools][parseAndFilter_DEgenes][call]} {input.tsv} {output.out} {params.variousParams}'

if not 'DGIDB_IN' in globals():
    DGIDB_IN = OUTDIR + 'databaseQuery/'
if not 'DGIDB_OUT' in globals():
    DGIDB_OUT = OUTDIR + 'databaseQuery/'

# query identified variants at dgidb
rule dgidbQuery:
    input:
        infile = DGIDB_IN + '{sample}.{clusterid}.txt'
    output:
        outfile = DGIDB_OUT + '{sample}.{clusterid}.dgidb.txt',
        outfileCompleteTable = DGIDB_OUT + '{sample}.{clusterid}.dgidb.txt.CompleteTable.txt',
        outfileGeneCategory = DGIDB_OUT + '{sample}.{clusterid}.dgidb.txt.GeneCategories.txt'
    params:
        lsfoutfile = DGIDB_OUT + '{sample}.{clusterid}.dgidbQuery.lsfout.log',
        lsferrfile = DGIDB_OUT + '{sample}.{clusterid}.dgidbQuery.lsferr.log',
        scratch = config['tools']['queryDGIDB']['scratch'],
        mem = config['tools']['queryDGIDB']['mem'],
        time = config['tools']['queryDGIDB']['time'],
        minDatabaseNum = config['tools']['queryDGIDB']['minDatabaseNum'],
        colName_genes = config['tools']['queryDGIDB']['colName_genes']
    threads:
        config['tools']['queryDGIDB']['threads']
    benchmark:
        DGIDB_OUT + '{sample}.{clusterid}.dgidbQuery.benchmark'
    shell:
         '{config[tools][queryDGIDB][call]} {input.infile} {output.outfile} {params.minDatabaseNum} {params.colName_genes}'


if not 'CLINICALTRIALS_IN' in globals():
    CLINICALTRIALS_IN = DGIDB_OUT
if not 'CLINICALTRIALS_OUT' in globals():
    CLINICALTRIALS_OUT = OUTDIR + 'clinicalTrials/'

# clnical trials query
rule clinicalTrialsQuery:
    input:
        infile = CLINICALTRIALS_IN + '{sample}.{clusterid}.dgidb.txt.CompleteTable.txt',
        downloadSuccess = CLINICALTRIALS_OUT + 'downloadSuccess.txt',
        #clinicalTrialsFolder = {config['resources']['general']['clinicalTrialsFolder']}
    output:
        outfile = CLINICALTRIALS_OUT + '{sample}.{clusterid}.dgidb.txt.CompleteTable.ClinicalTrials.allDrugs.txt'
    params:
        lsfoutfile = CLINICALTRIALS_OUT + '{sample}.{clusterid}.clinicalTrialsQuery.lsfout.log',
        lsferrfile = CLINICALTRIALS_OUT + '{sample}.{clusterid}.clinicalTrialsQuery.lsferr.log',
        scratch = config['tools']['queryClinicalTrials']['scratch'],
        mem = config['tools']['queryClinicalTrials']['mem'],
        time = config['tools']['queryClinicalTrials']['time'],
        cancerType = config['tools']['downloadClinicalTrials']['cancerType'],
        whiteList = config['tools']['queryClinicalTrials']['whiteList'],
        blackList = config['tools']['queryClinicalTrials']['blackList'],
        outDirec = CLINICALTRIALS_OUT
    threads:
        config['tools']['queryClinicalTrials']['threads']
    benchmark:
        CLINICALTRIALS_OUT + '{sample}.{clusterid}.clinicalTrialsQuery.benchmark'
    shell:
        '{config[tools][queryClinicalTrials][call]} {input.infile} {output.outfile} {params.outDirec}/{params.cancerType}_clinicalTrials/ "{params.whiteList}" "{params.blackList}"'
        #'{config[tools][queryClinicalTrials][call]} {input.infile} {output.outfile} {input.clinicalTrialsFolder}/'


# download the clinical trials necessary for the query
rule downloadClinicalTrials:
    output:
        outfile = CLINICALTRIALS_OUT + 'downloadSuccess.txt'
    params:
        lsfoutfile = CLINICALTRIALS_OUT + 'downloadClinicalTrials.lsfout.log',
        lsferrfile = CLINICALTRIALS_OUT + 'downloadClinicalTrials.lsferr.log',
        scratch = config['tools']['downloadClinicalTrials']['scratch'],
        mem = config['tools']['downloadClinicalTrials']['mem'],
        time = config['tools']['downloadClinicalTrials']['time'],
        cancerType = config['tools']['downloadClinicalTrials']['cancerType'],
        outDirec = CLINICALTRIALS_OUT
    threads:
        config['tools']['downloadClinicalTrials']['threads']
    benchmark:
        CLINICALTRIALS_OUT + 'downloadClinicalTrials.benchmark'
    shell:
        ('wget "https://clinicaltrials.gov/search?term={params.cancerType}&studyxml=true" ' + 
        '-O {params.outDirec}/{params.cancerType}_clinicalTrials.zip ; ' + 
        'unzip {params.outDirec}/{params.cancerType}_clinicalTrials.zip -d {params.outDirec}/{params.cancerType}_clinicalTrials ; ' +
        'touch {params.outDirec}/downloadSuccess.txt')


if not 'ANNOTATECLINICAL_IN' in globals():
    ANNOTATECLINICAL_IN = PROCESSDE_OUT
if not 'ANNOTATECLINICAL_OUT' in globals():
    ANNOTATECLINICAL_OUT = OUTDIR + 'clinicalAnnotation/'

# Combine different database queries, annotate input table with clinical information
rule annotate_DE_clinicalInformation:
    input:
        infile = ANNOTATECLINICAL_IN + '{sample}.{clusterid}.txt',
        pathwayDB = config['resources']['pathwayDB'],
        inDGIDB = DGIDB_OUT + '{sample}.{clusterid}.dgidb.txt.GeneCategories.txt',
        inClinicalTrials = CLINICALTRIALS_OUT + '{sample}.{clusterid}.dgidb.txt.CompleteTable.ClinicalTrials.allDrugs.txt'
    output:
        outTable = ANNOTATECLINICAL_OUT + '{sample}.{clusterid}.clinicalAnnotation.txt',
        outTable_dgdidbIndependent = ANNOTATECLINICAL_OUT + '{sample}.{clusterid}.clinicalAnnotation.txt_dgidbIndependent.txt'
    params:
        lsfoutfile = ANNOTATECLINICAL_OUT + '{sample}.{clusterid}.annotateClinical.lsfout.log',
        lsferrfile = ANNOTATECLINICAL_OUT + '{sample}.{clusterid}.annotateClinical.lsferr.log',
        scratch = config['tools']['annotateClinical']['scratch'],
        mem = config['tools']['annotateClinical']['mem'],
        time = config['tools']['annotateClinical']['time'],
        variousParams = config['tools']['annotateClinical']['variousParams']
    threads:
        config['tools']['annotateClinical']['threads']
    benchmark:
        ANNOTATECLINICAL_OUT + '{sample}.{clusterid}.annotateClinical.benchmark'
    shell:
        '{config[tools][annotateClinical][call]} --inputTable {input.infile} --outFile {output.outTable} --pathwayDB {input.pathwayDB} --dgidb_categ {input.inDGIDB} --clinTrials {input.inClinicalTrials} {params.variousParams}'


if not 'CIVIC_IN' in globals():
    CIVIC_IN = ANNOTATECLINICAL_OUT
if not 'CIVIC_OUT' in globals():
    CIVIC_OUT = ANNOTATECLINICAL_OUT

# query identified expression in civic
rule queryCIVIC:
    input:
        infile = CIVIC_IN +  '{sample}.{clusterid}.clinicalAnnotation.txt'
    output:
        outfile = CIVIC_OUT + '{sample}.{clusterid}.clinicalAnnotation.civic.txt'
    params:
        lsfoutfile = CIVIC_OUT + '{sample}.{clusterid}.queryCIVIC.lsfout.log',
        lsferrfile = CIVIC_OUT + '{sample}.{clusterid}.queryCIVIC.lsferr.log',
        scratch = config['tools']['queryCIVIC']['scratch'],
        mem = config['tools']['queryCIVIC']['mem'],
        time = config['tools']['queryCIVIC']['time'],
        cancerType = config['tools']['queryCIVIC']['cancerType'],
        blackList = config['tools']['queryCIVIC']['blackList'],
        highLevel = config['tools']['queryCIVIC']['highLevel'],
        colName_gene = config['tools']['queryCIVIC']['colName_gene'],
        colName_logFC = config['tools']['queryCIVIC']['colName_logFC'],
        strictExpression = config['tools']['queryCIVIC']['strictExpression']
    threads:
        config['tools']['queryCIVIC']['threads']
    benchmark:
        CIVIC_OUT + '{sample}.{clusterid}.queryCIVIC.benchmark'
    shell:
        '{config[tools][queryCIVIC][call]} --inputTable {input.infile} --outFile {output.outfile} --cancerTypeList "{params.cancerType}" --blackList "{params.blackList}" --highLevelList "{params.highLevel}" --colName_gene {params.colName_gene} --colName_logFC {params.colName_logFC} --strictExpression {params.strictExpression}'


if not 'GENESETANALYSIS_IN' in globals():
    GENESETANALYSIS_IN = DIFF_EXP_OUT
if not 'GENESETANALYSIS_OUT' in globals():
    GENESETANALYSIS_OUT = OUTDIR + 'geneSetAnalysis/'

# clnical trials query
rule geneSetEnrichment:
    input:
        infile = GENESETANALYSIS_IN + '{sample}.{clusterid}.DEgenes.tsv'
    output:
        outfile = GENESETANALYSIS_OUT + '{sample}.{clusterid}.enrichedGeneSets.txt'
    params:
        lsfoutfile = GENESETANALYSIS_OUT + '{sample}.{clusterid}.geneSetAnalysis.lsfout.log',
        lsferrfile = GENESETANALYSIS_OUT + '{sample}.{clusterid}.geneSetAnalysis.lsferr.log',
        scratch = config['tools']['geneSetEnrichment']['scratch'],
        mem = config['tools']['geneSetEnrichment']['mem'],
        time = config['tools']['geneSetEnrichment']['time'],
        geneSetDB = config['resources']['genesets'],
        variousParams = config['tools']['geneSetEnrichment']['variousParams']
    threads:
        config['tools']['geneSetEnrichment']['threads']
    benchmark:
        GENESETANALYSIS_OUT + '{sample}.{clusterid}.geneSetAnalysis.benchmark'
    shell:
        '{config[tools][geneSetEnrichment][call]} {input.infile} {output.outfile} {params.geneSetDB} {params.variousParams}'

def getGeneSetHeatmapFiles(wildcards):
    if "DEmalignant" in wildcards.sample:
        if len(CLUSTER_IDS_MALIGNANT) == 0:
            return [GENESETANALYSIS_OUT]
        return expand(GENESETANALYSIS_OUT + wildcards.sample + '.{clusterid}.enrichedGeneSets.txt', clusterid = CLUSTER_IDS_MALIGNANT)
    return expand(GENESETANALYSIS_OUT + wildcards.sample + '.{clusterid}.enrichedGeneSets.txt', clusterid = CLUSTER_IDS)


# plot heat map for gene set enrichment
# NOTE: empty inputs cause snakemake to crash, thus we implemented the workaround that returns GENESETANALYSIS_OUT in case of an empty list of enrichment analyses comparing tumor sub clones
rule plot_gene_set_enrichment:
    input:
        inDir = getGeneSetHeatmapFiles
        #inDir = expand(GENESETANALYSIS_OUT + '{{sample}}.{clusterid}.enrichedGeneSets.txt', clusterid = CLUSTER_IDS)
    output:
        outfile = GENESETANALYSIS_OUT + '{sample}.heatmap_enrichment.png'
    params:
        lsfoutfile = GENESETANALYSIS_OUT + '{sample}.heatmap_enrichment.lsfout.log',
        lsferrfile = GENESETANALYSIS_OUT + '{sample}.heatmap_enrichment.lsferr.log',
        scratch = config['tools']['plotGeneSetEnrichment']['scratch'],
        mem = config['tools']['plotGeneSetEnrichment']['mem'],
        time = config['tools']['plotGeneSetEnrichment']['time'],
        variousParams = config['tools']['plotGeneSetEnrichment']['variousParams'],
        comparison_direc = GENESETANALYSIS_OUT
    threads:
        config['tools']['plotGeneSetEnrichment']['threads']
    benchmark:
        GENESETANALYSIS_OUT + '{sample}.heatmap_enrichment.benchmark'
    shell:
        'if [ "{input.inDir}" != "{params.comparison_direc}" ] ; then echo "test1" ; {config[tools][plotGeneSetEnrichment][call]} {output.outfile} {input.inDir} ; else touch {output.outfile} ; fi'


if not 'DRUGCOMBINATION' in globals():
    DRUGCOMBINATION = OUTDIR + 'drugCombination/'

if not 'PARSETRIALSTABLE_IN' in globals():
    PARSETRIALSTABLE_IN = CLINICALTRIALS_OUT
if not 'PARSETRIALSTABLE_OUT' in globals():
    PARSETRIALSTABLE_OUT = DRUGCOMBINATION

# parse the *.dgidb.txt.CompleteTable.ClinicalTrials.txt files of all clusters
# and generate a table that shows for each Drug the clusters that can be targeted by this drug
# and the weight of the drug for calculating minimum set cover

rule parseDgidbTrialsTable_for_minSetCover:
    input:
        infiles = expand(PARSETRIALSTABLE_IN + '{{sample}}.{clusterid}.dgidb.txt.CompleteTable.ClinicalTrials.{{type}}.txt', clusterid = CLUSTER_IDS),
        drugList = config['resources']['drugList']
    output:
        out = PARSETRIALSTABLE_OUT + '{sample}.drugToCluster.{type}.txt'
    params:
        colName_clinTrial = config['tools']['parseDgidbTrialsTable_for_minSetCover']['colName_clinTrial'],
        colName_DGIDB_score = config['tools']['parseDgidbTrialsTable_for_minSetCover']['colName_DGIDB_score'],
        lsfoutfile = PARSETRIALSTABLE_OUT + '{sample}.drugToCluster.{type}.lsfout.log',
        lsferrfile = PARSETRIALSTABLE_OUT + '{sample}.drugToCluster.{type}.lsferr.log',
        scratch = config['tools']['parseDgidbTrialsTable_for_minSetCover']['scratch'],
        mem = config['tools']['parseDgidbTrialsTable_for_minSetCover']['mem'],
        time = config['tools']['parseDgidbTrialsTable_for_minSetCover']['time']
    threads:
        config['tools']['parseDgidbTrialsTable_for_minSetCover']['threads']
    benchmark:
        PARSETRIALSTABLE_OUT + '{sample}.drugToCluster.{type}.benchmark'
    shell:
        '{config[tools][parseDgidbTrialsTable_for_minSetCover][call]} --inFiles {input.infiles} --outFile {output.out} --colName_clinTrial {params.colName_clinTrial} --colName_DGIDB_Score {params.colName_DGIDB_score} --drug_list {input.drugList}'

if not 'MINSETCOVER_IN' in globals():
    MINSETCOVER_IN = PARSETRIALSTABLE_OUT
if not 'MINSETCOVER_OUT' in globals():
    MINSETCOVER_OUT = DRUGCOMBINATION

# deduce minimum set cover to find drug (combination) that targets all clusters via interaction with DE gene of this cluster

rule findminSetCover:
    input:
        infile = MINSETCOVER_IN + '{sample}.drugToCluster.{type}.txt',
        percTable = PERCENTAGE_OUT + '{sample}.clusters_cell_count_percent.txt'
    output:
        out = MINSETCOVER_OUT + '{sample}.drugCombination.{type}.txt'
    params:
        lsfoutfile = MINSETCOVER_OUT + '{sample}.drugCombination.{type}.lsfout.log',
        lsferrfile = MINSETCOVER_OUT + '{sample}.drugCombination.{type}.lsferr.log',
        scratch = config['tools']['findminSetCover']['scratch'],
        mem = config['tools']['findminSetCover']['mem'],
        time = config['tools']['findminSetCover']['time'],
        variousParams = config['tools']['findminSetCover']['variousParams']
    threads:
        config['tools']['findminSetCover']['threads']
    benchmark:
        MINSETCOVER_OUT + '{sample}.drugCombination.{type}.benchmark'
    shell:
        '{config[tools][findminSetCover][call]} --input {input.infile} --outFile {output.out} --percentageTable {input.percTable} {params.variousParams}'


if not 'FILTERDRUGS_IN' in globals():
    FILTERDRUGS_IN = CLINICALTRIALS_OUT
if not 'FILTERDRUGS_OUT' in globals():
    FILTERDRUGS_OUT = CLINICALTRIALS_OUT
# filter the drug-gene interaction results for drugs that are included in TP melanoma clinical list of drugs
rule filterDrugs:
    input:
        infile = FILTERDRUGS_IN + '{sample}.{clusterid}.dgidb.txt.CompleteTable.ClinicalTrials.allDrugs.txt',
        drugList = config['resources']['drugList']
    output:
        out = FILTERDRUGS_OUT + '{sample}.{clusterid}.dgidb.txt.CompleteTable.ClinicalTrials.filteredDrugs.txt'
    params:
        lsfoutfile = FILTERDRUGS_OUT + '{sample}.{clusterid}.filterDrugs.lsfout.log',
        lsferrfile = FILTERDRUGS_OUT + '{sample}.{clusterid}.filterDrugs.lsferr.log',
        scratch = config['tools']['filterDrugs']['scratch'],
        mem = config['tools']['filterDrugs']['mem'],
        time = config['tools']['filterDrugs']['time'],
    threads:
        config['tools']['filterDrugs']['threads']
    benchmark:
        FILTERDRUGS_OUT + '{sample}.{clusterid}.filterDrugs.benchmark'
    shell:
        '{config[tools][filterDrugs][call]} --inFile {input.infile} --outFile {output.out} --drugList {input.drugList}'


if not 'PREPROCESSUPSETR_IN' in globals():
    PREPROCESSUPSETR_IN = DRUGCOMBINATION
if not 'PREPROCESSUPSETR_OUT' in globals():
    PREPROCESSUPSETR_OUT = DRUGCOMBINATION

# this rule is to preprocess the output of the rule parseDgidbTrialsTable_for_minSetCover of type drug,clusters,weight (tab separated)
# to have the necessary input format for the UpSetR package, which is drug,cluster1,cluster2..clustern (tab separated) with 1 or 0 indicating if a drug targets a DE of a cluster
rule preprocessForUpSetR_venn:
    input:
        infile = PREPROCESSUPSETR_IN + '{sample}.drugToCluster.{type}.txt'
    output:
        out = PREPROCESSUPSETR_OUT + '{sample}.drugToCluster.{type}.processedForUpSetR.txt'
    params:
        lsfoutfile = PREPROCESSUPSETR_OUT + '{sample}.drugToCluster.{type}.processedForUpSetR.lsfout.log',
        lsferrfile = PREPROCESSUPSETR_OUT + '{sample}.drugToCluster.{type}.processedForUpSetR.lsferr.log',
        scratch = config['tools']['preprocessForUpSetR']['scratch'],
        mem = config['tools']['preprocessForUpSetR']['mem'],
        time = config['tools']['preprocessForUpSetR']['time']
    threads:
        config['tools']['preprocessForUpSetR']['threads']
    benchmark:
        PREPROCESSUPSETR_OUT + '{sample}.drugToCluster.{type}.processedForUpSetR.benchmark'
    shell:
        '{config[tools][preprocessForUpSetR][call]} --inFile {input.infile} --outFile {output.out}'


if not 'PLOTUPSETR_IN' in globals():
    PLOTUPSETR_IN = DRUGCOMBINATION
if not 'PLOTUPSETR_OUT' in globals():
    PLOTUPSETR_OUT = DRUGCOMBINATION

# this rule generates a plot with UpSetR (comparably to venn diagramm)
# The plot displays intersections of the drug sets that target DE genes in the clusters
rule plotUpSetR_venn:
    input:
        infile = PLOTUPSETR_IN + '{sample}.drugToCluster.{type}.processedForUpSetR.txt'
    output:
        out = PLOTUPSETR_OUT + '{sample}.drugToCluster.{type}.vennplot.png'
    params:
        lsfoutfile = PLOTUPSETR_OUT + '{sample}.drugToCluster.{type}.vennplot.lsfout.log',
        lsferrfile = PLOTUPSETR_OUT + '{sample}.drugToCluster.{type}.vennplot.lsferr.log',
        scratch = config['tools']['plotUpSetR']['scratch'],
        mem = config['tools']['plotUpSetR']['mem'],
        time = config['tools']['plotUpSetR']['time'],
        variousParams = config['tools']['plotUpSetR']['variousParams']
    threads:
        config['tools']['plotUpSetR']['threads']
    benchmark:
        PLOTUPSETR_OUT + '{sample}.drugToCluster.{type}.vennplot.benchmark'
    shell:
        '{config[tools][plotUpSetR][call]} --inFile {input.infile} --outFile {output.out} {params.variousParams}'


if not 'FULL_DRUGLIST_TO_SUBCLONES_IN' in globals():
    FULL_DRUGLIST_TO_SUBCLONES_IN = PARSETRIALSTABLE_OUT
if not 'FULL_DRUGLIST_TO_SUBCLONES_OUT' in globals():
    FULL_DRUGLIST_TO_SUBCLONES_OUT = PARSETRIALSTABLE_OUT

# This rule generates a table with a full list of all clinically relevant drugs and the information
# for each cluster if a drug gene interaction was found in dgidb between the respective drug and 
# any differentially expressed gene of the cluster.

rule get_full_druglist_to_subclones:
    input:
        infile = FULL_DRUGLIST_TO_SUBCLONES_IN + '{sample}.drugToCluster.allDrugs.txt',
    output:
        out = FULL_DRUGLIST_TO_SUBCLONES_OUT + '{sample}.full_druglist_to_subclones.txt'
    params:
        lsfoutfile = FULL_DRUGLIST_TO_SUBCLONES_OUT + '{sample}.full_druglist_to_subclones.lsfout.log',
        lsferrfile = FULL_DRUGLIST_TO_SUBCLONES_OUT + '{sample}.full_druglist_to_subclones.lsferr.log',
        scratch = config['tools']['get_full_druglist_to_subclones']['scratch'],
        mem = config['tools']['get_full_druglist_to_subclones']['mem'],
        time = config['tools']['get_full_druglist_to_subclones']['time'],
        drugList = config['resources']['drugList']
    threads:
        config['tools']['get_full_druglist_to_subclones']['threads']
    benchmark:
        FULL_DRUGLIST_TO_SUBCLONES_OUT + '{sample}.full_druglist_to_subclones.benchmark'
    shell:
        '{config[tools][get_full_druglist_to_subclones][call]} --in_drugToCluster {input.infile} --in_drugList {params.drugList} --outFile {output.out} '


if not 'PLOT_DRUGS_IN' in globals():
    PLOT_DRUGS_IN = CIVIC_OUT
if not 'PLOT_DRUGS_OUT' in globals():
    PLOT_DRUGS_OUT = DRUGCOMBINATION

# check whether all civic queries are finished
def getCivicQueryResults(wildcards):
    return expand(PLOT_DRUGS_IN + wildcards.sample + '.{clusterid}.clinicalAnnotation.civic.txt', clusterid = CLUSTER_IDS)

# this rule generates a UMAP plot that shows drug prediction on the tumor clones
rule plot_drug_prediction:
    input:
        rdsFile = REMOVE_ATYPICAL_OUT + '{sample}.RDS',
        inFiles = getCivicQueryResults,
        drugList = config['resources']['drugList'],
        drugCombis = config['resources']['drugCombinations'],
        civicDict = config['resources']['civicDict']
    output:
        out = PLOT_DRUGS_OUT + '{sample}.drug_prediction_umap.png'
    params:
        lsfoutfile = PLOT_DRUGS_OUT + '{sample}.drug_prediction_umap.lsfout.log',
        lsferrfile = PLOT_DRUGS_OUT + '{sample}.drug_prediction_umap.lsferr.log',
        scratch = config['tools']['show_drugPrediction_on_clones']['scratch'],
        sampleName = '{sample}',
        inputDir = PLOT_DRUGS_IN,
        outputDirec = PLOT_DRUGS_OUT,
        mem = config['tools']['show_drugPrediction_on_clones']['mem'],
        time = config['tools']['show_drugPrediction_on_clones']['time'],
        variousParams = config['tools']['show_drugPrediction_on_clones']['variousParams']
    threads:
        config['tools']['show_drugPrediction_on_clones']['threads']
    benchmark:
        PLOT_DRUGS_OUT + '{sample}.drug_prediction_umap.benchmark'
    shell:
        '{config[tools][show_drugPrediction_on_clones][call]} --SCE {input.rdsFile} --drugPredDir {params.inputDir} --outputDirec {params.outputDirec} --drugList {input.drugList} --combiList {input.drugCombis} --civicDict {input.civicDict} --sampleName {params.sampleName} {params.variousParams}'

