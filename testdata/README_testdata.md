# Guidelines to run testdata with scAmpi

The following instructions provide a quick guideline to (i) download example data, (ii) install necessary software, (iii) prepare the config file, and (iv) run the basic scAmpi workflow.

In the following, we refer to your test directory as `testdir` and to the directory of the scAmpi git repository as `git_scAmpi`.

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






