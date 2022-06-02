from snakemake.utils import validate
from snakemake.utils import min_version
import pandas as pd
import os
from glob import glob


# do not allow subdirectories as part of 'sample' or 'i' wildcards
wildcard_constraints:
    sample="[^/]+",
    i="[^/]+",

# Define minimum Snakemake version
min_version("6.12.1")

# Include Config file
configfile: "config/config.yaml"

# This file includes common functions used in the pipeline
include: "rules/misc_snake.smk"
# Include rules
include: "rules/scAmpi_clinical_rules.smk"


# include local rules
localrules:
    clinical_only,

# infer cluster IDs from input files
SAMPLE,IDS, = glob_wildcards("results/parse_diff_exp/{sample}.{id}.txt")

# defines output of a full clinical run
rule clinical_only:
    input:
        # clinical annotation
        expand("results/clinical_annotation/{sample}.{clusterid}.clinicalAnnotation.txt", sample = sample_ids, clusterid = IDS),
        expand("results/query_civic/{sample}.{clusterid}.clinicalAnnotation.civic.txt", sample = sample_ids, clusterid = IDS),
        # parse_for_minSetCover (aggregation rule)
        expand("results/drug_combination/{sample}.drugToCluster.allDrugs.txt", sample = sample_ids),
        expand("results/drug_combination/{sample}.drugToCluster.filteredDrugs.txt", sample = sample_ids),
        # preprocess for upsetR plot
        expand("results/upsetr_plot/{sample}.drugToCluster.allDrugs.processedForUpSetR.txt", sample = sample_ids),
        expand("results/upsetr_plot/{sample}.drugToCluster.filteredDrugs.processedForUpSetR.txt", sample = sample_ids),
        # plot upset plot
        expand("results/upsetr_plot/{sample}.drugToCluster.allDrugs.vennplot.png", sample = sample_ids),
        expand("results/upsetr_plot/{sample}.drugToCluster.filteredDrugs.vennplot.png", sample = sample_ids),
        # druglist_to_subclones
        expand("results/drug_combination/{sample}.full_druglist_to_subclones.txt", sample = sample_ids),
    output:
        "results/complete.txt",
    shell:
        "date > {output}"


### Adaptation specific for clinical-only run

# parse the *.dgidb.txt.CompleteTable.ClinicalTrials.txt files of all clusters
# and generate table that shows for each Drug the clusters that can be targeted by this drug 
# and the weight of the drug for calculating minimum set cover
def get_parse_for_minSetCover_input_clinical_only(wildcards):
    return expand('results/clinical_trials/{sample}.{i}.dgidb.txt.CompleteTable.ClinicalTrials.{{type}}.txt',
            i = IDS,
            sample = wildcards.sample)

# replace the cellranger step that is the default in this pipeline with the open source starsolo
ruleorder: parse_for_minSetCover_clinical_only > parse_for_minSetCover
use rule parse_for_minSetCover as parse_for_minSetCover_clinical_only with:
    input:
        infiles = get_parse_for_minSetCover_input_clinical_only,
        drugList = config['resources']['drugList']
