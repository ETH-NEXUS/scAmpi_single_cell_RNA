#!/bin/bash

# read in fastqs from given directory
# create symbolic link in new folder in output directory
# example call: 
# sh prepare_fastqs_for_snakemake.sh singlecell_rna/openbis/ singlecell_rna/fastqs/

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

