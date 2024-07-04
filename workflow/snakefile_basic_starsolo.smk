from snakemake.utils import min_version
from snakemake.utils import validate
import pandas as pd
import os
from glob import glob


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
include: "rules/starsolo_rules.smk"


# include local rules
localrules:
    scAmpi_basic,


# final rule of pipeline
# defines output of scampi basic
rule scAmpi_basic:
    input:
        expand("results/starsolo/{sample}.features.tsv", sample=sample_ids),
        expand("results/counts_raw/{sample}.h5", sample=sample_ids),
        expand(
            "results/identify_doublets/{sample}.doublet_barcodes.txt", sample=sample_ids
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
        expand(
            "results/atypical_removed/{sample}.atypical_removed.RDS", sample=sample_ids
        ),
        expand(
            "results/gene_exp/{sample}.gene_expression_per_cluster.tsv",
            sample=sample_ids,
        ),
        expand("results/plotting/{sample}.celltype_barplot.png", sample=sample_ids),
        expand("results/gsva/{sample}.gsetscore_hm.png", sample=sample_ids),
        expand("results/diff_exp_analysis/{sample}/", sample=sample_ids),
    output:
        "results/complete.txt",
    shell:
        "date > {output}"


# replace the cellranger step that is the default in this pipeline with the open source starsolo
ruleorder: create_hdf5_starsolo > create_hdf5


use rule create_hdf5 as create_hdf5_starsolo with:
    input:
        genes_file="results/starsolo/{sample}.features.tsv",
        matrix_file="results/starsolo/{sample}.matrix.mtx",
        barcodes_file="results/starsolo/{sample}.barcodes.tsv",
    conda:
        "envs/create_hdf5.yaml"
