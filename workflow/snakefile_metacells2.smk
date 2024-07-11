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
min_version("8.14.0")


# Include Config file
configfile: "config/config.yaml"


# Include report functionality
report: "report/workflow.rst"


# This file includes common functions used in the pipeline
include: "rules/misc_snake.smk"
# Include rules
include: "rules/scAmpi_basic_rules.smk"
include: "rules/scAmpi_metacells2_rules.smk"


# include local rules
localrules:
    scAmpi_metacells2,


# final rule of pipeline
# defines output of scampi metacells
rule scAmpi_metacells2:
    input:       
        expand("results/plotting/{sample}.celltype_barplot.png", sample=sample_ids),
        expand("results/gsva/{sample}.gsetscore_hm.png", sample=sample_ids),
        expand("results/diff_exp_analysis/{sample}/", sample=sample_ids),
        expand("results/plotting/{sample}_metacells2.celltype_barplot.png", sample=sample_ids),
        expand("workflow/report/rules/metacells2/{sample}_celltyping_summary.rst", sample=sample_ids),
        expand("results/gsva/{sample}_metacells2.gsetscore_hm.png", sample=sample_ids),
        expand("results/gsva/{sample}.gsetscore_hm.png", sample=sample_ids),
   