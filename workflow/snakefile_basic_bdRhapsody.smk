from snakemake.utils import min_version
from snakemake.utils import validate
import pandas as pd
import os
from glob import glob


# Environment variable specifying where the singularity image is located
# can be set as: export SINGULARITY_CACHE=<path to where our singularity images are on customapps>
#envvars:
#    "SINGULARITY_CACHE"



# do not allow subdirectories as part of `sample` or `i` wildcards
wildcard_constraints:
    sample="[^/]+",
    i="[^/]+",


# Define minimum Snakemake version
min_version("6.12.1")


# Include Config file
configfile: "config/config.yaml"


# Include report functionality
# report: "../report/workflow.rst"


# This file includes common functions used in the pipeline
include: "rules/misc_snake.smk"
# Include rules
include: "rules/scAmpi_basic_rules.smk"
include: "rules/bdRhapsody_rules.smk"


# include local rules
localrules:
    scAmpi_basic,


# final rule of pipeline
# defines output of scampi basic
rule scAmpi_basic:
    input:
        expand("results/counts_raw/{sample}.h5", sample=sample_ids),
        expand(
            "results/counts_filtered/{sample}.doublet_barcodes.txt", sample=sample_ids
        ),
        expand(
            "results/counts_raw/{sample}.raw.histogram_library_sizes.png",
            sample=sample_ids,
        ),
        expand(
            "results/counts_filtered/{sample}.genes_cells_filtered.histogram_library_sizes.png",
            sample=sample_ids,
        ),
        expand("results/counts_corrected/{sample}.corrected.RDS", sample=sample_ids),
        expand("results/clustering/{sample}.clusters_phenograph.csv", sample=sample_ids),
# Removing atypical cells should not be standard in many situations
       # expand(
       #     "results/atypical_removed/{sample}.atypical_removed.RDS", sample=sample_ids
       # ),
        expand(
            "results/gene_exp/{sample}.gene_expression_per_cluster.tsv",
            sample=sample_ids,
        ),
        expand("results/plotting/{sample}.celltype_barplot.png", sample=sample_ids),
        "results/celltyping/cells_per_celltype_and_sample.txt",

        # both gsva and de analysis have not yet been tested and used on mouse bd rhapsody data.
#        expand("results/gsva/{sample}.gsetscore_hm.png", sample=sample_ids),
#        expand("results/diff_exp_analysis/{sample}/", sample=sample_ids),
    output:
        "results/complete.txt",
    shell:
        "date > {output}"


ruleorder: create_hdf5_bdrhapsody > create_hdf5


use rule create_hdf5 as create_hdf5_bdrhapsody with:
    input:
        genes_file="results/bdr_matrix/{sample}.features.tsv",
        matrix_file="results/bdr_matrix/{sample}.matrix.mtx",
        barcodes_file="results/bdr_matrix/{sample}.barcodes.tsv",
    conda:
        "envs/create_hdf5.yaml"
