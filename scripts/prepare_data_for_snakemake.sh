#!/bin/bash

# read in fastqs from given directory
# create symbolic link in new folder in output directory
# example call: 
# sh prepare_fastqs_for_snakemake.sh /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/robustness/M140905_rep1_10X_RNA_6000c_180502_NS500318_0464_AHVNHYBGX5/singlecell_rna/openbis/ /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/robustness/M140905_rep1_10X_RNA_6000c_180502_NS500318_0464_AHVNHYBGX5/singlecell_rna/fastqs/

inFastqDirec=$1*
outFastqDirec=$2

for f in $inFastqDirec.fastq.gz
do
    path=$(readlink -f $f)
    #echo $path
    fName=$(basename "$f")
    echo $fName
    newPath=$outFastqDirec$fName
    ln -s "$path" "$newPath"
done

