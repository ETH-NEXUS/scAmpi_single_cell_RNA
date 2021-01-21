#!/usr/bin/env python

import argparse
import subprocess
import json


import sys
print(sys.version)

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample_name", required=True, help="the sample name")
parser.add_argument("-p","--path", required=True, help="/the singlecell_rna/ path")
parser.add_argument("-c", "--config", required=True, help="path to the config file")
parser.add_argument("-t", "--cancer_type", required=True, help="cancer type, e.g. melanoma")
parser.add_argument("-n", "--cellranger_sample_name", required=True, help="Shortened sample name for cellranger, without dots"),
parser.add_argument("-f", "--fastqc", required=True, help="Path to fastqc files"),
parser.add_argument("-r", "--seqRunName", required=True, help="Name tag of the sequencing run generating the fastq files")

args = parser.parse_args()


def write_sample_map(s_name):
    sample_map = "1\t" + s_name + "\tN\t1"

    with open(path+'auto_sample_map.tsv', 'w') as f:
        f.write(sample_map)
    subprocess.call(['chmod','775',path+'auto_sample_map.tsv'])

def add_slash(path):
    if path[-1] != '/':
        path = path + '/'
    return path

def write_run_sh(s_name, path):
    print(path)

    run_sh = "bsub -J " + s_name + " -e " + path + "scTP_snake.err -o " + path + "scTP_snake.out -W 23:59 -R \"rusage[mem=5000]\" \" export HDF5_USE_FILE_LOCKING=\"FALSE\"; /cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/snakemake --notemp --latency-wait 60 -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/all/scTumorProfiler/snake/snake_scTranscriptomics_master.snake --configfile " + path + "config_scTranscriptomics.json --cluster \'bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}\' -j 100 -p -k ; /cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/snakemake --notemp --latency-wait 60 -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/all/scTumorProfiler/snake/snake_clinical_master.snake --configfile " + path + "config_scTranscriptomics.json --cluster \'bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}\' -j 100 -p -k\" "

    with open(path+'auto_run.sh', 'w') as f:
        f.write(run_sh)

    subprocess.call(['chmod','775',path+'auto_run.sh'])


def write_scT_run_sh(s_name, path):
    print(path)

    run_sh = "bsub -J " + s_name + " -e " + path + "scTP_snake.err -o " + path + "scTP_snake.out -W 23:59 -R \"rusage[mem=5000]\" \" export HDF5_USE_FILE_LOCKING=\"FALSE\"; /cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/snakemake --notemp --latency-wait 60 -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/all/scTumorProfiler/snake/snake_scTranscriptomics_master.snake --configfile " + path + "config_scTranscriptomics.json --cluster \'bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}\' -j 100 -p -k\" "

    with open(path+'auto_run_sctranscriptomics.sh', 'w') as f:
        f.write(run_sh)

    subprocess.call(['chmod','775',path+'auto_run_sctranscriptomics.sh'])


def write_clinical_run_sh(s_name, path):
    print(path)

    run_sh = "bsub -J " + s_name + " -e " + path + "scTP_snake.err -o " + path + "scTP_snake.out -W 23:59 -R \"rusage[mem=5000]\" \" export HDF5_USE_FILE_LOCKING=\"FALSE\"; /cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/snakemake --notemp --latency-wait 60 -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/all/scTumorProfiler/snake/snake_clinical_master.snake --configfile " + path + "config_scTranscriptomics.json --cluster \'bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}\' -j 100 -p -k\" "

    with open(path+'auto_run_clinical.sh', 'w') as f:
        f.write(run_sh)

    subprocess.call(['chmod','775',path+'auto_run_clinical.sh'])



def write_dryrun_sh(path):

    print(path)

    dryrun_sh = "/cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/snakemake --notemp --latency-wait 60 -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/all/scTumorProfiler/snake/snake_scTranscriptomics_master.snake \
    --configfile " + path + "config_scTranscriptomics.json --cluster \'bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}\' -j 100 -p -k -n"
    
    with open(path+'auto_dryrun_sctranscriptomics.sh', 'w') as f:
        f.write(dryrun_sh)
    subprocess.call(['chmod','775',path+'auto_dryrun_sctranscriptomics.sh'])


def write_clinical_dryrun_sh(path):

    print(path)

    dryrun_sh = "/cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/snakemake --notemp --latency-wait 60 -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/all/scTumorProfiler/snake/snake_clinical_master.snake \
    --configfile " + path + "config_scTranscriptomics.json --cluster \'bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}\' -j 100 -p -k -n"
    
    with open(path+'auto_dryrun_clinical.sh', 'w') as f:
        f.write(dryrun_sh)
    subprocess.call(['chmod','775',path+'auto_dryrun_clinical.sh'])


def write_config_json(path, config, cancer_type, cellranger_sample_name,fastqc_path,seqRunName):
    
    print(path)

    with open(config) as f:
        data = json.load(f)

    data['inputOutput']['input_fastqs'] = path
    data['inputOutput']['input_fastqc'] = fastqc_path
    data['inputOutput']['analysis_output_dir'] = path
    data['inputOutput']['analysis_temp_dir'] = path + "snake_temp/"
    data['inputOutput']['sample_map'] = path + "auto_sample_map.tsv"
    data['inputOutput']['malignant_cell_type'] = cancer_type
    data['inputOutput']['sequencing_runName'] = seqRunName
    data['tools']['cellranger_count']['cellranger_sampleName'] = cellranger_sample_name

    # import pdb
    # pdb.set_trace()

    with open(path+'config_scTranscriptomics.json', 'w') as fp:
        json.dump(data, fp, indent=4, separators=(',',':'))
    #  subprocess.call('cat config_scTranscriptomics.json | python -m json.tool > config_scTranscriptomics.json', shell=True)

path = add_slash(args.path)
config = args.config
cancer_type = args.cancer_type
cellranger_sample_name = args.cellranger_sample_name
s_name = args.sample_name
fastqc_path = args.fastqc
seqRunName = args.seqRunName

write_config_json(path, config, cancer_type, cellranger_sample_name, fastqc_path,seqRunName)
write_sample_map(s_name)
write_run_sh(s_name, path)
write_scT_run_sh(s_name, path)
write_clinical_run_sh(s_name, path)
write_dryrun_sh(path)
write_clinical_dryrun_sh(path)

