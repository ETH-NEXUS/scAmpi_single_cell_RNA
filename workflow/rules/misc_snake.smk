import os.path
import sys
import inspect
import copy
import pandas as pd
from snakemake.utils import validate

# import sample map and retrieve sample names
samples_table = pd.read_table(config["inputOutput"]["sample_map"], header = 0)
samples = samples_table.set_index("sample", drop = False)
sample_ids = samples_table["sample"].tolist()
validate(samples, "../schema/sample_map.schema.yaml")


#########################################
###  Check config file for missing values
#########################################

fail_instantly = False

# define error object and message
class Error(object):
    def __init__(self, key, name):
        self.__key = key
        self.__name = name

    def __add__(self, other):
        return self

    def __call__(self, wildcards=None):
        sys.exit(
            """
            ===============================================
            You have not specified '{}' for '{}'
            ===============================================
            """.format(self.__key, self.__name))

    def __getitem__(self, value):
        return Error(key=self.__key, name=self.__name)

# define config object
class Config(object):
    def __init__(self, kwargs, name='Config'):
        self.__name = name
        self.__members = {}
        for (key, value) in kwargs.items():
            if isinstance(value, dict):
                self.__members[key] = Config(kwargs=value, name=key)
            else:
                self.__members[key] = value

    def __getitem__(self, key):
        if key in self.__members:
            return self.__members[key]
        else:
            if fail_instantly:
                sys.exit(
                    """
                    ===============================================
                    You have not specified '{}' for '{}'
                    ===============================================
                    """.format(key, self.__name))
            else:
                return Error(key=key, name=self.__name)

# check with the above class definitions if the config file contains all necessary values
config = Config(config)



###   input function for local rules in snakefile.smk

# input function for local rule `check_output` in snakefile_basic.smk
# retrieves the info from the config file if only the basic part, or both, basic and clinical, should be run.
def define_output(wildcards):
    if config['inputOutput']['basic_only']:
        return 'results/finished/{sample}.scAmpi_basic.txt'
    else:
        return expand('results/finished/{sample}.scAmpi_{part}.txt', part = ['basic', 'clinical'], sample = wildcards.sample)


# input function for local rule `clinical_mode` in snakefile.smk
def count_clusters(wildcards):
    checkpoint_output = checkpoints.diff_exp_analysis.get(**wildcards).output[0]
    all_clusters = expand('results/diff_exp_analysis/{sample}/{sample}_{i}.DEgenes.tsv',
    i = glob_wildcards(os.path.join(checkpoint_output, "{sample,[^/]+}.{i,[^/]+}.DEgenes.tsv")).i,
    sample = wildcards.sample)
    all_count = len(all_clusters)

    malignant_clusters = expand('results/diff_exp_analysis/{sample}/vs_other_malignant/{sample}.DEmalignant.{i}.DEgenes.tsv',
    i = glob_wildcards(os.path.join(checkpoint_output, "vs_other_malignant/{sample,[^/]+}.DEmalignant.{i,[^/]+}.DEgenes.tsv")).i,
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
