import os.path
import sys
import inspect
import copy

fail_instantly = False

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

#class Config(object):
#    def __init__(self, kwargs, name='Config'):
#        self.__name = name
#        self.__members = {}
#        for (key, value) in kwargs.items():
#            if isinstance(value, dict):
#                self.__members[key] = Config(kwargs=value, name=key)
#            else:
#                self.__members[key] = value
#    
#    def __getitem__(self, key):
#        if key in self.__members:
#            return self.__members[key]
#        else:
#            if fail_instantly:
#                sys.exit(
#                    """
#                    ===============================================
#                    You have not specified '{}' for '{}'
#                    ===============================================
#                    """.format(key, self.__name))
#            else:
#                return Error(key=key, name=self.__name)

#config = Config(config)

def getSampleNames():
    output = [] #[samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
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
                    sample = lineSplit[1]
                    if not (sample in output):
                        output.append(sample)
    return output

def getExperimentNames():
    output = [] #[samplename.replace(FASTQDIR,'').replace('/','')for samplename in glob.glob(FASTQDIR + '*/')]
    if output == []:
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
                    exp = lineSplit[0]
                    if not (exp in output):
                        output.append(exp)
    return output

def getSampleNamesFromExperimentNames(wildcards):
    if not 'SAMPLEMAPPING' in globals():
        return ['NOMAPPINGFILE']
    try:
        open(SAMPLEMAPPING, "r")
    except IOError:
        return ['NOMAPPINGFILE']
    expMap = dict()
    with open(SAMPLEMAPPING, "r") as f:
        for line in f:
            if line.strip() != "":
                lineSplit = line.strip().split()
                exp = lineSplit[0]
                sample = lineSplit[1]
                sampleType = lineSplit[2]
                tpoint = lineSplit[3]
                if exp not in expMap.keys():
                    expMap[exp] = []
                expMap[exp].append(sample)
    return expMap[wildcards.experiment]

def checkFilesAgainstSampleNames(files, sampleNames):
    finalFiles = []
    for f in files:
        for name in sampleNames:
            if name + "/" == f[0:len(name+"/")]:
                finalFiles.append(f)

    return finalFiles

def getSingleFastqFiles(SAMPLENAMES):
    files = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
    if files == []:
        files = [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq')]

    return checkFilesAgainstSampleNames(files, SAMPLENAMES)

    #return [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq.gz')]
    #return [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/SINGLEEND/*.fastq')]

def getPairedFastqFiles(SAMPLENAMES):
    files = [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq.gz')]
    if files == []:
        files = [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq')]
   
    return checkFilesAgainstSampleNames(files, SAMPLENAMES)

    #return [file.replace(FASTQDIR, '').replace('.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq.gz')]
    #return [file.replace(FASTQDIR, '').replace('.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*R[12].fastq')]

def getPairedFastqFilesWithoutR(SAMPLENAMES):
    files = [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
    if files == []:
        files = [file.replace(FASTQDIR, '').replace('_R1.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq')]

    return checkFilesAgainstSampleNames(files, SAMPLENAMES)

    #return [file.replace(FASTQDIR, '').replace('_R1.fastq.gz','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq.gz')]
    #return [file.replace(FASTQDIR, '').replace('_R1.fastq','')for file in glob.glob(FASTQDIR + '*/PAIREDEND/*_R1.fastq')]

def getNormalTumorFiles():
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
                exp = lineSplit[0]
                sample = lineSplit[1]
                sampleType = lineSplit[2]
                tpoint = lineSplit[3]
                if exp not in sampleMap.keys():
                    sampleMap[exp] = dict()
                if tpoint not in sampleMap[exp].keys():
                    sampleMap[exp][tpoint] = dict()
                if sampleType not in sampleMap[exp][tpoint].keys():
                    sampleMap[exp][tpoint][sampleType] = []
                sampleMap[exp][tpoint][sampleType].append(sample)
    for expKey, expValue in sampleMap.items():
        for tpointKey, tpointValue in expValue.items():
            if 'T' in tpointValue and 'N' in tpointValue:
                for sampleTumor in tpointValue['T']:
                    for sampleNormal in tpointValue['N']:
                        output.append(sampleTumor + '_vs_' + sampleNormal)
    return output


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
