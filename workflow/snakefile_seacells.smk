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
report: "report/workflow.rst"


# make sure seacell specific rules are applied if sample wildcard is followed by _seacells in the file names
ruleorder: sctransform_preprocessing_filtered_seacells > sctransform_preprocessing


# This file includes common functions used in the pipeline
include: "rules/misc_snake.smk"
# Include rules
include: "rules/scAmpi_basic_rules.smk"
include: "rules/scAmpi_seacells_rules.smk"


# include local rules
localrules:
    scAmpi_seacells,


# final rule of pipeline
# defines output of scampi seacells
rule scAmpi_seacells:
    input:
        #expand("results/plotting/{sample}.celltype_barplot.png", sample=sample_ids),
        expand(
            "results/plotting/{sample}_seacells.celltype_barplot.png",
            sample=sample_ids,
        ),
        expand("results/gsva/{sample}.gsetscore_hm.png", sample=sample_ids),
        #expand("results/diff_exp_analysis/{sample}/", sample=sample_ids),
        expand(
            "workflow/report/rules/seacells/{sample}_celltyping_summary.rst",
            sample=sample_ids,
        ),
        expand("results/gsva/{sample}_seacells.gsetscore_hm.png", sample=sample_ids),
        "results/metacells/seacells_celltype_purity.png",
        expand("results/celltype_gsva_c6/{sample}_seacells_celltype_GSVA.tsv", sample=sample_ids),
        expand("results/celltype_gsva_c6/{sample}_celltype_GSVA.tsv", sample=sample_ids),
        
