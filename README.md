# scAmpi - Single Cell Analysis mRNA pipeline

#### General overview

This scAmpi workflow is organized into two main parts: the `scAmpi_basic` part and the `scAmpi_clinical` part, which can be run independently. scAmpi_basic includes general scRNA processing steps, such as mapping, QC, normalisation, unsupervised clustering, cell type classification, and DE analysis.

scAmpi_clinial includes the search for disease relevant drug targets for differentially expressed genes. Note that the clinical part is only applied if at least one cluster identified in your sample is indicated as a diseased ("malignant") cell type.

![README_rulegraph](https://user-images.githubusercontent.com/38692323/175028270-2ac20406-720d-4941-bfb9-e924e5f65759.png)


#### Installation instructions

scAmpi follows the best practices of the Snakemake workflow manager in providing the software needed to run the pipeline in per-rule conda environments. Those environmnents are specified in the `envs/` directory in yaml files that are named `{rule_name}.yaml`. The easiest way to install and use the software is by running Snakemake with the `--use-conda` parameter. Snakemake will try to find the environments of the yaml files the rules point to, and install them if they are not already available. The directory for installing the conda environments can be specified with the `--conda-prefix` parameter.

1. Make sure `snakemake` is in your PATH.
   Follow the instructions on how to install `snakemake` [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
2. Install all conda environments the workflow needs before running an analysis with

```bash
snakemake --use-conda --conda-create-envs-only --conda-prefix /my/directory/for/conda/envs/ -s workflow/snakefile_basic.smk --configfile config/config.yaml
```

- `--use-conda` instructs snakemake to utilize the `conda:` directive in the rules
- `--conda-create-envs-only` specifies that only the installation of conda environments is triggered, not the analysis of the samples.
- *(optional):* with `--conda-prefix /my/directory/for/conda/envs/` a directory for the installation of the conda environments can be specified.

#### Installations of tools for initial read mapping and counting
For the read mapping and UMI counting step scAmpi offers pre-defined rules for using either Cellranger or STARsolo. Both tools are not available for installation via conda and need to be installed separately. Only one of the tools needs to be installed, depending on the method of choice.
- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger): Follow the instructions on the 10xGenomics installation support page to install cellranger and to include it into the PATH.
Webpage: [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)
- [STAR](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) as open source alternative to Cellranger. For installation, follow the instructions in the excellent STAR documentation and include it in your PATH.

#### Example data

For a test run the freely available 10X Genomics data from PBMC cells can be used. A step by step guideline and example config file are provided in the directory `testdata/`. Note that this test run assumes that cellranger` has been chosen for read mapping.

#### Before running the pipeline

* **internet connection**
Some steps of the scAmpi workflow perform online queries. Please make sure that this is possible on your computing system, e.g. by loading the respective modules to enable the proxy connection. (Most systems will have this enabled per default).


* **config file**
  * input directory
    Before running the pipeline the `config.yaml` file needs to be adapted to contain the **path to input fastq files** for the intended analysis. It is provided in the
    first section (`inputOutput`) of the config file.
  * resource information
    In addition to the input path, further resource information must be provided in the section `resources`. This information is primarily specifying
    input required for the cell type classification and the genomic reference used for the cellranger mapping. An example `config.yaml` file ready for adaptation, as
    well as a brief description of the relevant config blocks, is provided in the directory `config/`.
* **sample map**
Provide a "sample_map", i.e. a tab delimited text file listing all samples that should be analysed (one row per sample).
The sample map must contain a column with the header `sample` (see example below). This ID will be used to name files and identify the sample throughout the pipeline.
An example file ready for adaptation is provided in the directory `config/`.

Sample map example:
```
sample
SAMPLE-1_scR
SAMPLE-2_scR
```



#### Running scAmpi

Different use cases of scAmpi are covered by several snakefiles to choose from.

- to run the basic part only: `workflow/snakefile_basic.smk`
- to run the basic and the clinical part together: `workflow/snakefile_clinical.smk`
- to run the clinical part only: `workflow/snakefile_clinical-only.smk`
- to run scAmpi with STAR instead of Cellranger: `workflow/snakefile_basic_starsolo.smk`

Please find details below.

### scAmpi_basic part

Example call:

```
snakemake -s workflow/snakefile_basic.smk --configfile config/config.yaml -j 1 -p
```

Note that if the pipeline is run on a compute cluster with a job scheduling system (e.g. LSF) the commands need to be adjusted accordingly.


### scAmpi_clinical part

Example call (that includes the basic part as well):

```
snakemake -s workflow/snakefile_clinical.smk --configfile config/config.yaml -j 1 -p
```

### A note on using CIViC

The CIViC query implemented in scAmpi makes use of an offline cache file of the CIViC database. The cache is retrieved with the initial installation of the scAmpi software. Afterwards, users have to manually update the cache file if they want to use a new version.
To update the cache file, load the scAmpi environemnt and open an R session.
Then type:

```
>> civic.update_cache()
```

### Running scAmpi_clinical independently

It is possible to run the scAmpi_clinical part independently of scAmpi_basic, following some restrictions to the file names and formatting.

- Use the master snake file `workflow/snakefile_clinical-only.smk`.
- scAmpi_clinical expects as input the results of a DE analysis on cell cluster level
- The input files must follow the file name convention `SAMPLEID.CLUSTER.txt`
  - SAMPLEID is the sample name specified in the sample map
  - CLUSTER is the cell cluster ID
  - `txt` is the expected suffix

- Provide input files in the subdirectory `results/parse_diff_exp/` that needs to be created.
- The input files must contain five mandatory columns:

```
gene_names  diff    padj      test_statistic  pct_nonzero
ATP1A1      1.679   3.05e-15  14.506          81.42
```
Here, "gene_names" contains the HGNC gene symbols, "diff" contains the fold change or a similar value, "padj" contains the adjusted p-value, "test_statistic" contains the value of the test statistics, and "pct_nonzero" contains the percentage of cells in this cluster with non-zero expression in the respective gene.
Results of this clinical pipeline run are the *in-silico* drug prediction and clinical annotations.
Other side results, e.g. the minimum set cover computation, the plotting of drug predictions on the UMAP, and the gene set enrichment analysis, cannot be created in an independent clinical run as they rely on additional input files generated by the scAmpi_basic part.

#### Adapting/Integrating rules in Snakemake
Snakemake is a Python-based workflow management system for building and executing pipelines. A pipeline is made up of ["rules"](snake/scAmpi_basic_rules.py) that represent single steps of the analysis. In a [yaml config file](config/config_scAmpi.yaml) parameters and rule-specific input can be adjusted to a new analysis without changing the rules. In a ["master" snake file](snake/snake_scAmpi_basic_master.snake) the desired end points of the analysis are specified. With the input and the desired output defined, Snakemake is able infer all steps that have to be performed in-between.

To change one of the steps, e.g. to a different software tool, one can create a new rule, insert a new code block into the config file, and include the input/output directory of this step in the master snake file. It is important to make sure that the format of the input and output of each rule is compatible with the previous and the subsequent rule. For more detailed information please have a look at the excellent [online documentation](https://snakemake.readthedocs.io/en/stable/index.html) of Snakemake.

#### Quick start using test data
To quickly start a scAmpi_basic run with PBMC test data you can follow the following steps:
- clone the scAmpi repository
- make sure you have `snakemake` in your PATH (see Installation instructions)
- prepare Cellranger software and reference directory
- update the path to the cellranger reference directory in `testdata/config.yaml`
- download example data from the 10xGenomics website (for more detailed instructions see `testdata/README_testdata.md`)
- *optional*: to circumvent the time-consuming mapping step create the directory `results/counts_raw/` in your scAmpi repository, copy the raw matrix `testdata/5k_pbmc_v3.h5.tar` into the directory, gunzip the file (e.g. `tar -xvf 5k_pbmc_v3.h5.tar`) and start the test run from this step.
- perform Snakemake dryrun to see a list of steps that will be performed
  `snakemake -s workflow/snakefile_basic.smk --configfile testdata/config.yaml -n -p`
- start analysis run
  `snakemake -s workflow/snakefile_basic.smk --configfile testdata/config.yaml -p`
