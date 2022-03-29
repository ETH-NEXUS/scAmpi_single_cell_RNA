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
include: "rules/final_rules.smk"

# include local rules
localrules: all

# final rule of pipeline
rule all:
    input:
        expand('results/finished/{sample}.complete.txt', sample = sample_ids)
    output:
        'results/complete.txt'
    shell:
        'date > {output}'
