# scAmpi - Single Cell Analysis mRNA pipeline

#### General overview

This scAmpi workflow is organized into two main parts: the `scAmpi_basic` part and the `scAmpi_clinical` part, which can be run independently. scAmpi_basic includes general scRNA processing steps, such as mapping, QC, normalisation, unsupervised clustering, cell type classification, and DE analysis.

scAmpi_clinial includes the search for disease relevant drug targets for differentially expressed genes. Note that the clinical part is only applied if at least one cluster identified in your sample is indicated as a diseased ("malignant") cell type.


#### Installation instructions

scAmpi provides a yml file to enable installing most software used in the default workflow as a conda environment (yaml file `scAmpi_scRNA_conda_env.yml` provided in the sub folder "envs").

Example (This will install the conda environment in your home. Note that the installation may take a while.):
```
> conda env create -f scAmpi_scRNA_conda_env.yml --name scAmpi_scRNA
```
To activate the environment, type:
```
> conda activate scAmpi_scRNA
```

Additionally required installations that are not available via conda:
- [Phenograph](https://github.com/dpeerlab/phenograph):
with activated conda environment:
```
> pip install PhenoGraph
```
- [CIViCpy](https://github.com/griffithlab/civicpy): 
with activated conda environment:
```
> pip install civicpy
```

- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger): Follow the instructions on the 10xGenomics support page and include the cellranger binary to your path.
Webpage: [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)

3: download local cache of the [CIViC](https://civicdb.org) database using functionality from `civicpy`.
Webpage: https://docs.civicpy.org/en/latest/install.html

First-time download of the CIViCpy cache:

    >> from civicpy import civic
    >> civic.load_cache(on_stale='ignore')

Note that updates to this file are not handled by scAmpi and is the user who should do this using:

    >> civic.update_cache()


*Temporary*
For the normalisation step the R package `sctransform` is used. As there is a known bug in the latest release, at the moment either an older package version (e.g. 0.2) or the development version should be installed. The development version can be installed within R using `remotes::install_github("ChristophH/sctransform@develop")`.


#### Example data

For a test run the freely available 10X Genomics data from PBMC cells can be used. Please find an example config file and in the directory `testdata`


#### Before running the pipeline

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
snakemake --notemp --latency-wait 60 -s snake_scAmpi_basic_master.snake --configfile config_scAmpi.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k
```

Note that the section that follows parameter `--cluster` denotes the cluster specific notation to indicate memory and timr ressources for a job. Please adapt according to the respective job scheduling system used.


#### scAmpi_clinical part

Example call:

```
snakemake --notemp --latency-wait 60 -s snake_scAmpi_clinical_master.snake --configfile config_scTranscriptomics.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k
```

Note that the section that follows parameter `--cluster` denotes the cluster specific notation to indicate memory and timr ressources for a job. Please adapt according to the respective job scheduling system used.

