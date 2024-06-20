import os
import sys
import inspect
import copy
import pandas as pd
from snakemake.utils import validate

# import sample map and retrieve sample names
samples_table = pd.read_table(config["inputOutput"]["sample_map"], header=0)
sample_map = samples_table.set_index("sample_name", drop=False)
sample_ids = samples_table["sample_name"].tolist()
# if samples specified in config, process only those
if "samples" in config:
    if isinstance(config["samples"], str):
        config["samples"]=[config["samples"]]
    assert all(sa in sample_ids for sa in config["samples"]),"specified samples in config which are not present in sample table"
    assert sample_ids, "empty sample ids list specified"
    sample_ids=config["samples"]
    sample_map=sample_map.loc[sample_ids]
# for the reporting, we need the samples in the config dict.
config["samples"]=sample_ids

file_stem_dict=dict(sample_map['file_stem'])
sample_name_dict={v:k for k,v in file_stem_dict.items()}

validate(sample_map, "../schema/sample_map.schema.yaml")


#########################################
###  Check config file for missing values
#########################################

# Define config object
# it behaves like a dict, but prints a helpfull error message if something is missing
class Config:
    def __init__(self, kwargs, name="Config"):
        self.__name = name
        self.__members = {}
        for key, value in kwargs.items():
            if isinstance(value, dict):
                self.__members[key] = Config(kwargs=value, name=f'{self.__name}->{key}')
            else:
                self.__members[key] = value

    def print_error(self, key):
        msg=f'You have not specified "{key}" for "{self.__name}"'
        logger.error(f"""
                ===============================================
                {msg}
                ===============================================
                """)
        return f'config entry "{key}" not found'
    def __getitem__(self, key):
        if key in self.__members:
            return self.__members[key]
        else:
            raise KeyError(self.print_error(key))
    def __setitem__(self, key, value):
        self.__members[key]=value
# Check with the above class definitions if the config file contains all necessary values
config = Config(config)


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

def get_symlink_names(wildcards):
    sa=wildcards.sample
    in_path = config['inputOutput']['input_fastqs']
    in_files = [f for f in os.listdir(in_path) if f.startswith(wildcards.sample)]
    targets = [f.replace(file_stem_dict[sa],sa, 1) for f in in_files ]
    return targets

def get_input_fastq(wildcards):
    fn=wildcards.link_filename
    sa=next(sa for sa in sample_ids if fn.startswith(sa))
    return fn.replace(sa,file_stem_dict[sa], 1)
    
