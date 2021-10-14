# scAmpi - Single Cell Analysis mRNA pipeline

#### General overview

This scAmpi workflow is organized into two main parts: the `scAmpi_basic` part and the `scAmpi_clinical` part, which can be run independently. scAmpi_basic includes general scRNA processing steps, such as mapping, QC, normalisation, unsupervised clustering, cell type classification, and DE analysis.

scAmpi_clinial includes the search for disease relevant drug targets for differentially expressed genes. Note that the clinical part is only applied if at least one cluster identified in your sample is indicated as a diseased ("malignant") cell type.


#### Installation instructions

scAmpi provides a yml file to enable installing most software used in the default workflow as a conda environment (yaml file `scAmpi_scRNA_conda_env.yml` provided in the subfolder "envs"). Additionally, you find one concrete example of the conda environment used to run scAmpi with detailed package versions (`scAmpi_scRNA_detailed.yml`).

Example (This will install the conda environment in your home. Note that the conda may take a while to install all software packages):
```
>> conda env create -f scAmpi_scRNA_conda_env.yml --name scAmpi_scRNA
```

To activate the environment, type:
```
>> conda activate scAmpi_scRNA
```

Additionally required installations that are not available via conda:
- [Phenograph](https://github.com/dpeerlab/phenograph):
To install, first activate the conda environment and then use pip install:

```
>> pip install PhenoGraph
```

- [CIViCpy](https://github.com/griffithlab/civicpy): 
To install, first activate the conda environment and then use pip install:

```
>> pip install civicpy
```

- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger): Follow the instructions on the 10xGenomics installation support page to install cellranger and to include the cellranger binary to your path.
Webpage: [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)


*Temporary*
For the normalisation step the R package `sctransform` is used. As there is a known bug in the latest release (`0.3.2` - 2020-12-16), at the moment it is recommended to use the development version. For more information please have a look at the github repository of [sctransform](https://github.com/ChristophH/sctransform).
To install, activate the scAmpi environment, open an R session, and type:

```
>> remotes::install_github("ChristophH/sctransform@develop")
```

#### Example data

For a test run the freely available 10X Genomics data from PBMC cells can be used. A step by step guideline and example config file are provided in the directory `testdata`.


#### Before running the pipeline

Note: some steps of the scAmpi workflow perform online queries. Please make sure that this is possible on your computing system, e.g. by loading the respective modules to enable the proxy connection. (Most systems will have this enabled per default).

Before running the pipeline the `config` file needs to be adapted to contain the input and output paths for the intended analysis. Those are provided in the first section (`inputOutput`) of the config file. In addition to input and output paths, further resource information must be provided in the section `resources`. This information is primarily specifying input required for the cell type classification and the genomic reference used for the cellranger mapping. An example config file ready for adaptation, as well as a brief description of the relevant config blocks, is provided in the directory `config`.

Further, a "sample_map" must be provided, a tab delimited text file that lists all samples that should be analysed (one row per sample).
The sample map contains four columns (refer to example below): column 1: experiment ID, column 2: sample name (this ID will be used to name files and identify the sample throughout the pipeline); column 3: sample_info (free-text), e.g. tumor/normal/tissue_type; column 4: time point.
Am example file ready for adpatation is provided in the directory 'snake'.

Example:
```
1	SAMPLE-1_scR	lung	1
2	SAMPLE-2_scR	skin	1
```


#### scAmpi_basic part

Example call:

```
snakemake -s snake_scAmpi_basic_master.snake --configfile config_scAmpi.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k
```

Note that the section that follows parameter `--cluster` denotes the cluster specific notation to indicate memory and timr ressources for a job. Please adapt according to the respective job scheduling system used.


#### scAmpi_clinical part

Example call:

```
snakemake -s snake_scAmpi_clinical_master.snake --configfile config_scTranscriptomics.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k
```

Note that the section that follows parameter `--cluster` denotes the cluster specific notation to indicate memory and timr ressources for a job. Please adapt according to the respective job scheduling system used.

#### A note on using CIViC

The CIViC query implemented in scAmpi makes use of an offline cache file of the CIViC database. The cache is retrieved with the initial installation of the scAmpi software. Afterwards, users have to manually update the cache file if they want to use a new version. 
To update the cache file, load the scAmpi environemnt and open an R session.
Then type:
```
>> civic.update_cache()
```
