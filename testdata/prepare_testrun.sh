#!/bin/bash

#################################################
#### File name: prepare_testrun.sh
#### Author: Franziska Singer
#### Date created: March 2021
####
#### Create folders for test run + download data
####
#### Input $1: location of scAmpi git repository
#### Input $2: location of the cellranger reference
##################################################

cellranger_ref=$1

echo "Directory with cellranger reference: $cellranger_ref"


# replace reference directory path in testdata/config.yaml

# temporarily adapt IFS to preserve leading whitespaces
OLD_IFS="$IFS"
IFS=

while read -r line; do
    echo ${line//refdata-cellranger-GRCh38-3.0.0/${cellranger_ref}}
    done < testdata/config.yaml > testdata/config_temp.yaml
    mv testdata/config_temp.yaml testdata/config.yaml

# restore IFS to default
IFS="$OLD_IFS"

# download and extract example data from the 10xGenomics webpage

mkdir -p fastqs
cd fastqs
#wget https://cg.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar
#tar -xvf 5k_pbmc_v3_fastqs.tar
#rm 5k_pbmc_v3_fastqs.tar
#mv 5k_pbmc_v3_fastqs/* .
#rm -r 5k_pbmc_v3_fastqs/
