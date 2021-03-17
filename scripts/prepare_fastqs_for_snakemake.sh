#!/bin/bash

# read in fastqs from given directory
# create symbolic link in new folder in output directory
# example call: 
# sh prepare_fastqs_for_snakemake.sh singlecell_rna/openbis/ singlecell_rna/fastqs/

inFastqDirec=$1*
outFastqDirec=$2
#fastqOutNameTag=$3

for f in $inFastqDirec.fastq.gz
do
    path=$(readlink -f $f)
    #echo $path
    fName=$(basename "$f")
    echo $fName
    #ending=$(echo $fName | awk -F 'L00' '{print $2}')
    #echo $ending
    #nameParts=$(echo $fName | tr "_L00" "\n")
    #newPath=$outFastqDirec$fastqOutNameTag
    #echo $newPath
    #symFile=$newPath'_L00'$ending
    symFile=$outFastqDirec$fName
    echo $symFile
    #symCall=$('ln -s "$path" "$symFile"')
    #echo $symCall
    ln -s "$path" "$symFile"
done

