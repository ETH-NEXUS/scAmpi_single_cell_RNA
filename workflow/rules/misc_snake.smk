import os.path
import sys
import inspect
import copy
import pandas as pd
from snakemake.utils import validate

# import sample map to infer sample names
samples_table = pd.read_table(config["inputOutput"]["sample_map"], header = 0)
samples = samples_table.set_index("sample", drop = False)
sample_ids = samples_table["sample"].tolist()
validate(samples, "../schema/sample_map.schema.yaml")


### The following functions are used in the index hopping removal framework

# Generate sample-specific root paths for all samples in the sample map
# The root path of each sample is specified in the sixth column of the sample map
# The cellranger sample name for each sample is specified in the third column of the sample map
def getSampleRootPaths():
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    output = []
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                if len(lineSplit) != 6:
                    raise ValueError("Only %s fields were provided in the sample map where 6 were expected!" %(str(len(lineSplit))))
                sample = lineSplit[1].strip()
                crsample = lineSplit[2].strip()
                rootdir = lineSplit[5].strip()
                if sample in sampleMap.keys():
                    raise ValueError("Sample '%s' is not unique in the sample map!" %(sample))
                sampleMap[sample] = {}
                # NOTE: it is assumed that the provided root path already ends with "/"
                if not rootdir.endswith('/'):
                    raise ValueError("Path to root directory in the sample map must end with '/': %s" %(str(rootdir)))
                # Match each sample to its corresponding cellranger sample name
                sampleMap[sample][crsample] = str(rootdir)

    # Generate sample-specific paths
    # Rootdir is determined by the given sample, while the remaining path structure is fixed
    for sample in sampleMap.keys():
        for crsample in sampleMap[sample].keys():
            rootdir = sampleMap[sample][crsample]
            # NOTE: it is assumed that the provided root path already ends with "/"
            output.append(rootdir + crsample + '/singlecell_rna/analysis/cellranger_run/' + sample)
            # eg. Mysample-L_scR_v1.6
            # where:    rootdir = sample_dir/ (depends on the sample type)
            #           crsample = Mysample-L (short version of sample name)
            #           /singlecell_rna/analysis/cellranger_run/ (fixed)
            #           add actual sample name at the end (to be used as a filename tag) = Mysample-L_scR_v1.6
    return output


# Retrieve the cellranger sample name corresponding to a given sample in the sample map
# The cellranger name for each sample is specified in the third column of the sample map
def getCellRangerName(wildcards):
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    sampleMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                sample = lineSplit[1].strip()
                crsample = lineSplit[2].strip()
                if sample in sampleMap.keys():
                    raise ValueError("Sample '%s' is not unique in the sample map!" %(wildcards.sample))
                sampleMap[sample] = crsample

    if wildcards.sample not in sampleMap.keys():
        raise ValueError("Sample '%s' not found in the sample map!" %(wildcards.sample))
    return sampleMap[wildcards.sample]
