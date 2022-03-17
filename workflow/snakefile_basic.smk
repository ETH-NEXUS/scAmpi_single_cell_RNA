from snakemake.utils import min_version
from snakemake.utils import validate
import pandas as pd
import os
from glob import glob

# do not allow subdirectories as part of `sample` or `i` wildcards
wildcard_constraints:
    sample = "[^/]+",
    i = "[^/]+"

# Define minimum Snakemake version
min_version("6.12.1")

# Include Config file
configfile: "config/config.yaml"

# Include report functionality
#report: "../report/workflow.rst"

# This file includes common functions used in the pipeline
include: "rules/misc_snake.smk"

# Include rules
include: "rules/scAmpi_basic_rules.smk"
include: "rules/scAmpi_clinical_rules.smk"

# include local rules
localrules: all,check_output,scAmpi_basic

# final rule of pipeline
rule all:
    input:
        expand('results/finished/{sample}.complete.txt', sample = sample_ids)
    output:
        'results/complete.txt'
    shell:
        'date > {output}'


# input function for local rule `check_output` in snakefile_basic.smk
# retrieves the info from the config file if only the basic part, or both, basic and clinical, should be run.
def define_output(wildcards):
    if config['inputOutput']['basic_only']:
        return 'results/finished/{sample}.scAmpi_basic.txt'
    else:
        return expand('results/finished/{sample}.scAmpi_{part}.txt', part = ['basic', 'clinical'], sample = wildcards.sample)

# rule that checks if both parts of scampi, or the basic part only should be run
rule check_output:
    input:
        define_output
    output:
        'results/finished/{sample}.complete.txt'

# defines output of scampi basic
rule scAmpi_basic:
    input:
#        'results/cellranger_run/{sample}.features.tsv',
#        'results/counts_raw/{sample}.h5',
#        'results/counts_filtered/{sample}.doublet_barcodes.txt',
#        'results/counts_raw/{sample}.h5.histogram_library_sizes.png',
#        'results/counts_filtered/{sample}.genes_cells_filtered.h5.histogram_library_sizes.png',
#        'results/counts_corrected/{sample}.corrected.RDS',
#        'results/clustering/{sample}.clusters_phenograph.csv',
        'results/atypical_removed/{sample}.atypical_removed.RDS',
        'results/clustering/{sample}.clusters_cell_count_percent.txt',
        'results/gene_exp/{sample}.gene_expression_per_cluster.tsv',
        'results/plotting/{sample}.celltype_barplot.png',
        'results/gsva/{sample}.gsetscore_hm.png',
        'results/diff_exp_analysis/{sample}/',
    output:
        'results/finished/{sample}.scAmpi_basic.txt'
    shell:
        'date > {output}'


# input function for local rule `clinical_mode` in snakefile_basic.smk
def count_clusters(wildcards):
    checkpoint_output = checkpoints.diff_exp_analysis.get(**wildcards).output[0]
    all_clusters = expand('results/diff_exp_analysis/{sample}/{sample}_{i}.DEgenes.tsv',
    i = glob_wildcards(os.path.join(checkpoint_output, "{sample,[^/]+}.{i,[^/]+}.DEgenes.tsv")).i,
    sample = wildcards.sample)
    all_count = len(all_clusters)

    malignant_clusters = expand('results/diff_exp_analysis/{sample}/vs_other_malignant/{sample}.DEmalignant.{i}.DEgenes.tsv',
    i = glob_wildcards(os.path.join(checkpoint_output, "{sample,[^/]+}.DEmalignant.{i,[^/]+}.DEgenes.tsv")).i,
    sample = wildcards.sample)
    malignant_count = len(malignant_clusters)

    # get list of all clusters
    if  all_count > 0:
        return expand('results/finished/{sample}.clinical_full.txt',
        sample = wildcards.sample)
    elif malignant_count > 0:
        return expand('results/finished/{sample}.clinical_malignant_only.txt',
        sample = wildcards.sample)
    else:
        return expand('results/finished/{sample}.clinical_nonmalignant.txt',
        sample = wildcards.sample)

# rule that checks if the full clinical workflow can be run, or only a subset of the rules.
# this depends on the number of clusters found in the sample and cannot be determined beforehand.
rule clinical_mode:
    input:
        count_clusters
    output:
        'results/finished/{sample}.scAmpi_clinical.txt'
    shell:
        'echo {input.my_input}'


# defines output of a full clinical run
rule scAmpi_basic:
rule clinical_full:
    input:
        # trigger clinical part of the pipeline
        'results/aggregated/{sample}.aggregated.txt',
        # plot_drug_prediction is also aggregation rule (as is aggregate)
        'results/plot_drug_prediction/{sample}.drug_prediction_umap.png',
        # plot gene set enrichment heatmap (is also aggregation rule)
        'results/gene_set_enrichment/{sample}.heatmap_enrichment.png',
        'results/gene_set_enrichment/vs_other_malignant/{sample}.DEmalignant.heatmap_enrichment.png',
        # parse_for_minSetCover (is also aggregation rule)
        'results/drug_combination/{sample}.drugToCluster.allDrugs.txt',
        'results/drug_combination/{sample}.drugToCluster.filteredDrugs.txt',
        # preprocess for upsetR plot
        'results/upsetr_plot/{sample}.drugToCluster.allDrugs.processedForUpSetR.txt',
        'results/upsetr_plot/{sample}.drugToCluster.filteredDrugs.processedForUpSetR.txt',
        # plot upset plot
        'results/upsetr_plot/{sample}.drugToCluster.allDrugs.vennplot.png',
        'results/upsetr_plot/{sample}.drugToCluster.filteredDrugs.vennplot.png',
        # find minSetCover
        'results/drug_combination/{sample}.drugCombination.allDrugs.txt',
        'results/drug_combination/{sample}.drugCombination.filteredDrugs.txt',
        # druglist_to_subclones
        'results/drug_combination/{sample}.full_druglist_to_subclones.txt',
    output:
        'results/finished/{sample}.clinical_full.txt'
    shell:
        'date > {output}'


# defines output of a reduced clinical run.
# this is triggered if either no malignant or no non-malignant cells are found in the sample.
rule clinical_malignant_only:
    input:
        # druglist_to_subclones
        'results/drug_combination/{sample}.full_druglist_to_subclones.txt',
    output:
        'results/finished/{sample}.clinical_malignant_only.txt'
    shell:
        'date > {output}'


# defines output of a reduced clinical run.
# this is triggered if either no malignant or no non-malignant cells are found in the sample.
rule clinical_nonmalignant:
    input:
        # druglist_to_subclones
        'results/drug_combination/{sample}.full_druglist_to_subclones.txt',
    output:
        'results/finished/{sample}.clinical_nonmalignant.txt'
    shell:
        'date > {output}'
