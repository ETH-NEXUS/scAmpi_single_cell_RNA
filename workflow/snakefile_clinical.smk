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
localrules: all

# final rule of pipeline
rule all:
    input:
        expand("results/finished/{sample}.scAmpi_clinical.txt", sample = sample_ids)
    output:
        'results/complete.txt'
    shell:
        'date > {output}'


# defines output of scampi basic
rule scAmpi_basic:
    input:
#        'results/cellranger_run/{sample}.features.tsv',
#        "results/counts_raw/{sample}.h5",
#        "results/counts_filtered/{sample}.doublet_barcodes.txt",
        "results/counts_raw/{sample}.h5.histogram_library_sizes.png",
        "results/counts_filtered/{sample}.genes_cells_filtered.h5.histogram_library_sizes.png",
#        "results/counts_corrected/{sample}.corrected.RDS",
#        "results/clustering/{sample}.clusters_phenograph.csv",
#        "results/atypical_removed/{sample}.atypical_removed.RDS",
        "results/gene_exp/{sample}.gene_expression_per_cluster.tsv",
        "results/plotting/{sample}.celltype_barplot.png",
        "results/gsva/{sample}.gsetscore_hm.png",
        "results/diff_exp_analysis/{sample}/",
    output:
        "results/finished/{sample}.scAmpi_basic.txt",
    resources:
        mem_mb=config["computingResources"]["mem"]["low"],
        time_min=config["computingResources"]["time"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    shell:
        "date > {output}"


# rule that checks if the full clinical workflow can be run, or only a subset of the rules.
# this depends on the number of clusters found in the sample and cannot be determined beforehand.
rule clinical_mode:
    input:
        count_clusters,
    output:
        "results/finished/{sample}.scAmpi_clinical.txt",
    resources:
        mem_mb=config["computingResources"]["mem"]["low"],
        time_min=config["computingResources"]["time"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    shell:
        "echo {input}"


# defines output of a full clinical run
rule clinical_full:
    input:
        # trigger basic part of scAmpi
        "results/finished/{sample}.scAmpi_basic.txt",

        # trigger clinical part of scAmpi
        "results/aggregated/{sample}.aggregated.txt",
        # plot_drug_prediction is also aggregation rule (as is aggregate)
        "results/plot_drug_prediction/{sample}.drug_prediction_umap.png",
        # plot gene set enrichment heatmap (is also aggregation rule)
        "results/gene_set_enrichment/{sample}.heatmap_enrichment.png",
        "results/gene_set_enrichment/vs_other_malignant/{sample}.DEmalignant.heatmap_enrichment.png",
        # parse_for_minSetCover (is also aggregation rule)
        "results/drug_combination/{sample}.drugToCluster.allDrugs.txt",
        "results/drug_combination/{sample}.drugToCluster.filteredDrugs.txt",
        # preprocess for upsetR plot
        "results/upsetr_plot/{sample}.drugToCluster.allDrugs.processedForUpSetR.txt",
        "results/upsetr_plot/{sample}.drugToCluster.filteredDrugs.processedForUpSetR.txt",
        # plot upset plot
        "results/upsetr_plot/{sample}.drugToCluster.allDrugs.vennplot.png",
        "results/upsetr_plot/{sample}.drugToCluster.filteredDrugs.vennplot.png",
        # find minSetCover
        "results/drug_combination/{sample}.drugCombination.allDrugs.txt",
        "results/drug_combination/{sample}.drugCombination.filteredDrugs.txt",
        # druglist_to_subclones
        "results/drug_combination/{sample}.full_druglist_to_subclones.txt",
    output:
        "results/finished/{sample}.clinical_full.txt",
    resources:
        mem_mb=config["computingResources"]["mem"]["low"],
        time_min=config["computingResources"]["time"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    shell:
        "date > {output}"


# defines output of a reduced clinical run.
# this is triggered if either no malignant or no non-malignant cells are found in the sample.
rule clinical_malignant_only:
    input:
        # trigger basic part of scAmpi
        "results/finished/{sample}.scAmpi_basic.txt",

        # trigger reduced clinical part of scAmpi
        # plot gene set enrichment heatmap (is also aggregation rule)
        "results/gene_set_enrichment/vs_other_malignant/{sample}.DEmalignant.heatmap_enrichment.png",
    output:
        "results/finished/{sample}.clinical_malignant_only.txt",
    resources:
        mem_mb=config["computingResources"]["mem"]["low"],
        time_min=config["computingResources"]["time"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    shell:
        "date > {output}"


# defines output of a reduced clinical run.
# this is triggered if either no malignant or no non-malignant cells are found in the sample.
rule clinical_nonmalignant:
    input:
        # trigger basic part of scAmpi
        "results/finished/{sample}.scAmpi_basic.txt",
        # no steps of scAmpi clinical can be run
    output:
        "results/finished/{sample}.clinical_nonmalignant.txt",
    resources:
        mem_mb=config["computingResources"]["mem"]["low"],
        time_min=config["computingResources"]["time"]["low"],
    threads: config["computingResources"]["threads"]["low"]
    shell:
        "date > {output}"
