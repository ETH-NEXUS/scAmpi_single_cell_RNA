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

scAmpi_path=$1
cellranger_ref=$2
testrun_dir=$(pwd)

echo "Testrun directory: $testrun_dir"
echo "scAmpi git directory: $scAmpi_path"
echo "Directory with cellranger reference: $cellranger_ref"

# Create folders for testdata and analysis output

mkdir -p fastqs
mkdir -p analysis
mkdir -p snake_temp
mkdir -p snake_files

# copy example config file and sample map to the testrun directory and replace placeholders with correct path

cp ${scAmpi_path}/testdata/config_scAmpi_testdata.yaml snake_files/
cp ${scAmpi_path}/testdata/sample_map_testdata.tsv snake_files/

# replace test directory path
# temporarily adapt IFS to preserve leading whitespaces
OLD_IFS="$IFS"
IFS=
while read -r line; do
    echo ${line//testdir/${testrun_dir}}
    done < snake_files/config_scAmpi_testdata.yaml > snake_files/config_scAmpi_testdata_temp.yaml
    mv snake_files/config_scAmpi_testdata_temp.yaml snake_files/config_scAmpi_testdata.yaml

# replace reference directory path
while read -r line; do
    echo ${line//reference_dir/${cellranger_ref}}
    done < snake_files/config_scAmpi_testdata.yaml > snake_files/config_scAmpi_testdata_temp.yaml
    mv snake_files/config_scAmpi_testdata_temp.yaml snake_files/config_scAmpi_testdata.yaml

# replace scAmpi git path
while read -r line; do
    echo ${line//path_to_scAmpi_git/${scAmpi_path}}
    done < snake_files/config_scAmpi_testdata.yaml > snake_files/config_scAmpi_testdata_temp.yaml
    mv snake_files/config_scAmpi_testdata_temp.yaml snake_files/config_scAmpi_testdata.yaml

# restore IFS to default
IFS="$OLD_IFS"


# go to snake_files and prepare dry run and full run scripts
cd ${testrun_dir}/snake_files/
echo "snakemake -s ${scAmpi_path}/snake/snake_scAmpi_basic_master.snake --configfile config_scAmpi_testdata.yaml -n" > dryrun_scAmpi.sh
chmod +x dryrun_scAmpi.sh

echo "bsub -J scAmpi_test -eo ${testrun_dir}/snake_files/scAmpi_test.err -oo ${testrun_dir}/snake_files/scAmpi_test.out -W 23:59 \"snakemake --latency-wait 60 -s ${scAmpi_path}/snake/snake_scAmpi_basic_master.snake --configfile config_scAmpi_testdata.yaml --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R \"rusage[mem={params.mem},scratch={params.scratch}]\" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 10 -p -k\"" > run_scAmpi.sh
chmod +x run_scAmpi.sh

echo "snakemake -s ${scAmpi_path}/snake/snake_scAmpi_basic_master.snake --configfile config_scAmpi_testdata.yaml -j 10 -p -k" > run_scAmpi_local.sh
chmod +x run_scAmpi_local.sh


# download and extract example data from the 10xGenomics webpage
cd ${testrun_dir}/fastqs/
wget https://cg.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar
tar -xvf 5k_pbmc_v3_fastqs.tar
rm 5k_pbmc_v3_fastqs.tar
mv 5k_pbmc_v3_fastqs/* .
rm -r 5k_pbmc_v3_fastqs/
