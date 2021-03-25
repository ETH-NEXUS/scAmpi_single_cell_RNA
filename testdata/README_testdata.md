# Guidelines to run testdata with scAmpi

The following instructions provide a short guideline to download example data and execute the scAmpi baic scRNA-seq analysis workflow. 

#, (ii) install necessary software, (iii) prepare the config file, and (iv) run the basic scAmpi workflow.

In the following, we refer to your test directory as `testdir` and to the directory of the scAmpi git repository as `path_git_scAmpi`. Note that the testdata are human PBMCs, so you need to have a version of the cellranger human reference data files available. (For download instructions from cellranger, refer to the [10xGenomics Webpage](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation).
In the following we assume that the reference files are located in `reference_dir`.

(i) Activate environment and export cellranger path

For instructions on how to install the necessary software to run scAmpi, please refer to the README.md file in the scAmpi git repository.

In order to do a test run, you need to activate the conda environment and you need to export the path to the cellranger installation to you path. The following example calls assume that the conda environemt was installed in you HOME idrectory and that cellranger (version 6.0.0) was installed in "cellranger_dir". Note: if you do not want to add cellranger to your PATH, you can also provide the full path to the cellranger installation in the config file at the cellranger_run config block.

```
> conda avtivate scAmpi_scRNA
> export PATH=/cellranger_dir/cellranger-6.0.0:$PATH
```

(ii) Run script to prepare test directory
The script `prepare_testrun.sh` will create the fastqs and analysis directories, as well as downloads the testdata, and prepares the config file for the test run:

```
> cd testdir
> sh [path_git_scAmpi]/testdata/prepare_testrun.sh [path_git_scAmpi] [reference_dir]
```

### (i) Create sub folder and download example data

Create a folder for the fastq files, the analysis, and the genome reference in `testdir`:

```
> mkdir fastqs
> mkdir analysis
> mkdir reference
```

Download and untar 10X Genomics PBMC example data:

```
> cd fastqs
> wget https://cg.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar
> tar -xvf 5k_pbmc_v3_fastqs.tar
```

The testdata files should look like this:

```
> ls 5k_pbmc_v3_fastqs/
5k_pbmc_v3_S1_L001_I1_001.fastq.gz
5k_pbmc_v3_S1_L001_R1_001.fastq.gz
5k_pbmc_v3_S1_L001_R2_001.fastq.gz
5k_pbmc_v3_S1_L002_I1_001.fastq.gz
5k_pbmc_v3_S1_L002_R1_001.fastq.gz
5k_pbmc_v3_S1_L002_R2_001.fastq.gz
5k_pbmc_v3_S1_L003_I1_001.fastq.gz
5k_pbmc_v3_S1_L003_R1_001.fastq.gz
5k_pbmc_v3_S1_L003_R2_001.fastq.gz
5k_pbmc_v3_S1_L004_I1_001.fastq.gz
5k_pbmc_v3_S1_L004_R1_001.fastq.gz
5k_pbmc_v3_S1_L004_R2_001.fastq.gz
```

Download the 10x Genomics reference sequence (e.g. hg38)

```
> cd reference
> 
```

### (ii) Prepare installation

You can skip this section if you already installed all software necessary for scAmpi and activated the conda environment.

In order to run scAmpi, you need to install the related conda environment `scAmpi_scRNA_conda_env.yml` provided in the sub folder "envs".

```
> conda env create -f scAmpi_scRNA_conda_env.yml --name scAmpi_scRNA
```

Activate the environment:
```
> conda activate scAmpi_scRNA
```

Install Phenogrpah:
```
> pip install PhenoGraph
```

Install cellranger:
follow the steps on the 10xGenomics webpage:
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

NOTE:
*Temporary* addition to accomodate for sctransform bug:
Open an R session and install the development version of sctransform:
```
remotes::install_github("ChristophH/sctransform@develop")
```

### (iii) Prepare config file

All resources required to run scAmpi on the testdata have been included in the scAmpi repository. The config file provided in the sub folder testdata (`git_scAmpi/testdata/config_scAmpi_testdata.json`) already contains the relative paths to the respective resources. 

The relative paths in the config file need to be adapted to point to the paths of `git_scAmpi` and `testdir`.

### (iv) Run scAmpi on testdata






