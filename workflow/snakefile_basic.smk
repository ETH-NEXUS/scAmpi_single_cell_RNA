import os
import glob
import sys
import datetime

# Include Config file
configfile: "config/config.yaml"


# This function adapts the config object to include full path information
include: "rules/misc_snake.smk"

config["inputOutput"]["sample_map"]

# old framework
SAMPLEMAPPING = config['inputOutput']['sample_map']


# Include the rules
include: "rules/scAmpi_basic_rules.smk"

# This rule defines which files should be created
localrules: scAmpi_basic
rule scAmpi_basic:
    input:
        expand('results/cellranger_run/{sample}.features.tsv', sample = getSampleNames()),
        expand('results/counts_raw/{sample}.h5', sample = getSampleNames()),
        expand('results/counts_raw/{sample}.h5.histogram_library_sizes.png', sample = getSampleNames()),
        expand('results/counts_filtered/{sample}.genes_cells_filtered.h5.histogram_library_sizes.png', sample = getSampleNames()),
        expand('results/counts_corrected/{sample}.genes_cells_filtered.corrected.RDS', sample = getSampleNames()),
        expand('results/clustering/{sample}.genes_cells_filtered.corrected.clusters_phenograph.csv', sample = getSampleNames()),
        expand('results/atypical_removed/{sample}.genes_cells_filtered.corrected.atypical_removed.RDS', sample = getSampleNames()),
        expand('results/clustering/{sample}.genes_cells_filtered.corrected.atypical_removed.clusters_cell_count_percent.txt', sample = getSampleNames()),
        expand('results/diff_exp/{sample}.genes_cells_filtered.corrected.atypical_removed.diffExp_success.txt', sample = getSampleNames()),
        expand('results/gene_exp/{sample}.genes_cells_filtered.corrected.atypical_removed.gene_expression_per_cluster.tsv', sample = getSampleNames()),
        expand('results/plotting/{sample}.genes_cells_filtered.corrected.atypical_removed.celltype_barplot.png', sample = getSampleNames()),
        expand('results/gsva/{sample}.genes_cells_filtered.corrected.atypical_removed.gsetscore_hm.png', sample = getSampleNames()),
    output:
        'results/complete_scAmpi_basic.txt'
    params:
        lsfoutfile = 'results/complete_scAmpi_basic.lsfout.log',
        lsferrfile = 'results/complete_scAmpi_basic.lsferr.log',
        mem = '1000',
        scratch = '1000',
        time = '1'
    benchmark:
        'results/complete_scAmpi_basic.txt.benchmark'
    shell:
        'date > {output}'
