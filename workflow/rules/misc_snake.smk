import os
from os.path import join, dirname, splitext, isdir, basename
from glob import glob
import sys
import inspect
import copy
import pandas as pd
from snakemake.utils import validate


#########################################
###  Check config file for missing values
#########################################


# Define config object
# it behaves like a dict, but prints a helpfull error message if something is missing
class ConfigDict(dict):
    def __init__(self, *args, parent_path=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.parent_path = parent_path or ["Config"]
    def __getitem__(self, key):
        try:
            # logger.debug(f"get {self.full_path} -> {key}")
            value = super().__getitem__(key)
        except KeyError:
            logger.error(f"KeyError: Did not find '{key}' in {self.full_path}")
            raise
        return value
    def __setitem__(self, key, value):
        if isinstance(value, dict):
            # Ensure that nested dicts are also wrapped in ConfigDict.
            value = ConfigDict.from_dict(value, parent_path=self.parent_path + [key])
        super().__setitem__(key, value)
    @property 
    def full_path(self):
        return ' -> '.join(map(str, self.parent_path ))
    @classmethod
    def from_dict(cls, dct, parent_path=None):
        """Convert a regular dictionary to ConfigDict, recursively."""
        parent_path = parent_path or ["Config"]
        obj = cls(parent_path=parent_path)        
        for key, value in dct.items():
            obj[key] = value        
        return obj


# Check with the above class definitions if the config file contains all necessary values
config = ConfigDict.from_dict(config)


cr_version = config["tools"]["cellranger_count"]["version"]
assert cr_version in ["7.1.0", "8.0.1"], f"Unsuppoerted cellranger version {cr_version}"

#################################################
### import sample map and retrieve sample names
#################################################

sample_table = pd.read_table(config["inputOutput"]["sample_map"], sep="\t", header=0)
sample_ids = set(sample_table["sample"])
# if samples specified in config, process only those
if "samples" in config:
    if isinstance(config["samples"], str):
        config["samples"] = [config["samples"]]
    assert all(
        sa in sample_ids for sa in config["samples"]
    ), "specified samples in config which are not present in sample table"
    assert sample_ids, "empty sample ids list specified"
    sample_ids = set(config["samples"])
    sample_table = sample_table.query("sample in @sample_ids")
# for the reporting, we need the samples in the config dict.
config["samples"] = sample_ids

logger.info(f"processing {len(sample_ids)} samples: {', '.join(sample_ids)}")

fq_files = {sa: [] for sa in sample_ids}
link_names = {}
fastq_dir = config["inputOutput"]["input_fastqs"]
if not isdir(fastq_dir):
    logger.error(f"Cannot find fastq directory {fastq_dir}")
elif "file_stem" in sample_table:
    for _, row in sample_table.iterrows():
        sa = row["sample"]
        fs = row["file_stem"]
        if not isdir(join(fastq_dir, sa)):
            found = [
                join(fastq_dir, f)
                for f in os.listdir(fastq_dir)
                if f.startswith(fs) and f.endswith(".fastq.gz")
            ]
        else:
            found = [
                join(fastq_dir, sa, f)
                for f in os.listdir(join(fastq_dir, sa))
                if f.startswith(fs) and f.endswith(".fastq.gz")
            ]
        fq_files[sa].extend(found)
        for fq in found:
            link_names[fq] = basename(fq).replace(fs, sa)
            logger.debug(f"{fs=}, {sa=}, {link_names[fq]=}")
    fq_link_dict = {v: k for k, v in link_names.items()}
else:  # folder structure is already fastq_dir/sample/sample_S1_L001_I1_001.fastq.gz
    for _, row in sample_table.iterrows():
        sa = row["sample"]
        if isdir(join(fastq_dir, sa)):
            fq_files[sa].extend(
                [
                    join(fastq_dir, sa, f)
                    for f in os.listdir(join(fastq_dir, sa))
                    if f.startswith(sa) and f.endswith(".fastq.gz")
                ]
            )
for sa, fq_list in fq_files.items():
    logger.info(f"found {len(fq_list)} fastq files for sample {sa}")


validate(sample_table, "../schema/sample_map.schema.yaml")


#################################################
### Resources
#################################################


max_mem_mb=config["computingResources"]["max_mem_mb"]


#################################################
### helper functions for rules
#################################################


# input function for local rule `clinical_mode` in snakefile.smk
# With this function one of the three clinical output rules are triggered, depending on the number of clusters found.
def count_clusters(wildcards):
    checkpoint_output = checkpoints.diff_exp_analysis.get(**wildcards).output[0]
    all_clusters = expand(
        "results/diff_exp_analysis/{sample}/{sample}_{i}.DEgenes.tsv",
        i=glob_wildcards(
            os.path.join(checkpoint_output, "{sample,[^/]+}.{i,[^/]+}.DEgenes.tsv")
        ).i,
        sample=wildcards.sample,
    )

    all_count = len(all_clusters)

    malignant_clusters = expand(
        "results/diff_exp_analysis/{sample}/vs_other_malignant/{sample}.DEmalignant.{i}.DEgenes.tsv",
        i=glob_wildcards(
            os.path.join(
                checkpoint_output,
                "vs_other_malignant/{sample,[^/]+}.DEmalignant.{i,[^/]+}.DEgenes.tsv",
            )
        ).i,
        sample=wildcards.sample,
    )

    malignant_count = len(malignant_clusters)

    # get list of all clusters
    if all_count > 0:
        # then there is at least one malignant cluster and at least two non-malignant as is requested by the DE script to run
        return expand(
            "results/finished/{sample}.clinical_full.txt", sample=wildcards.sample
        )
    elif malignant_count > 0:
        # there are no non-malignant clusters but several malignant
        return expand(
            "results/finished/{sample}.clinical_malignant_only.txt",
            sample=wildcards.sample,
        )
    else:
        # there are no malignant clusters and therefore no DE in this sample (or in total only one cluster)
        return expand(
            "results/finished/{sample}.clinical_nonmalignant.txt",
            sample=wildcards.sample,
        )


def get_extension(filename):
    """
    Extracts the extension from a filename.
    Handles compound extensions like '.fq.gz' or '.fastq.gz'.
    """
    _, ext = splitext(filename)  # Split off the last extension
    if ext == ".gz":
        _, ext2 = splitext(filename[:-3])  # Split again if the file is gzipped
        ext = ext2 + ext
    return ext


def get_fastq_links(wildcards):
    sa = wildcards.sample
    fq_list = fq_files[sa]
    if not link_names:
        return fq_list
    return [join("results/input_fastq", sa, link_names[fq]) for fq in fq_list]


def get_params_remove_atypical_cells(wildcard, key):
    what = "cells"
    if "_seacells" in wildcard["sample"] or "_metacells2" in wildcard["sample"]:
        what = "metacells"
    return config["tools"]["remove_atypical_cells"][what][key]
