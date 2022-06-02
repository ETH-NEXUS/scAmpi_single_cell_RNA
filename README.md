# scAmpi - Single Cell Analysis mRNA pipeline

#### General overview

This scAmpi workflow is organized into two main parts: the `scAmpi_basic` part and the `scAmpi_clinical` part, which can be run independently. scAmpi_basic includes general scRNA processing steps, such as mapping, QC, normalisation, unsupervised clustering, cell type classification, and DE analysis.

scAmpi_clinial includes the search for disease relevant drug targets for differentially expressed genes. Note that the clinical part is only applied if at least one cluster identified in your sample is indicated as a diseased ("malignant") cell type.


![scAmpi_both_rulegraphs](https://user-images.githubusercontent.com/38692323/140029020-6292b989-722d-4c93-909d-1d65c8aacddd.png)



#### Installation instructions

scAmpi follows the best practices of the Snakemake workflow manager in providing the software needed to run the pipeline in per-rule conda environments. Those environmnents are specified in the `envs/` directory in yaml files that are named `{rule_name}.yaml`. The easiest way to install and use the software is by running Snakemake with the `--use-conda` parameter. Snakemake will try to find the environments of the yaml files the rules point to, and install them if they are not already available. The directory for installing the conda environments can be specified with the `--conda-prefix` parameter.

```bash
# prerequisite for the following steps is the snakemake software in your path

# before running scAmpi the software can be installed with the `--conda-create-envs-only`
snakemake --use-conda --conda-prefix /directory/for/conda/envs/ --conda-create-envs-only

#without the `--conda-create-envs-only` parameter the software will be installed and the pipeline steps subsequently.
```



#### Installations of tools for initial read mapping and counting

For the read mapping and UMI counting step scAmpi offers pre-defined rules for using either Cellranger or STARsolo. Both tools are not available for installation via conda and need to be installed separately. Only one of the tools needs to be installed, depending on the method of choice.

- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger): Follow the instructions on the 10xGenomics installation support page to install cellranger and to include the cellranger binary to your path.
  Webpage: [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)
- [STAR](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) as open source alternative to Cellranger. For installation, follow the instructions in the excellent STAR documentation and include STAR in your path upon running the pipeline.

#### Example data

For a test run the freely available 10X Genomics data from PBMC cells can be used. A step by step guideline and example config file are provided in the directory `testdata/`. Note that this test run assumes that the method cellranger has been chosen for read mapping.

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


#### scAmpi_basic part

Example call:

```
snakemake -s workflow/snakefile_basic.smk --configfile config/config.yaml -j 1 -p
```

Note that if the pipeline is run on a compute cluster with a job scheduling system (e.g. LSF) the commands need to be adjusted accordingly (e.g. for lsf: `snakemake --cluster`, indicating memory and time resources for a job (e.g. `--cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}'` )) or a [cluster profile](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) can be used.


#### scAmpi_clinical part

Example call:

```
snakemake -s workflow/snakefile_clinical.smk --configfile config/config.yaml -j 1 -p
```

#### A note on using CIViC

The CIViC query implemented in scAmpi makes use of an offline cache file of the CIViC database. The cache is retrieved with the initial installation of the scAmpi software. Afterwards, users have to manually update the cache file if they want to use a new version.
To update the cache file, load the scAmpi environemnt and open an R session.
Then type:

```
>> civic.update_cache()
```

#### Running scAmpi_clinical independently

It is possible to run the scAmpi_clinical part independently of scAmpi_basic, following some restrictions to the file names and formatting. For this use case please use the master snake file `workflow/snakefile_clinical-only.smk`. Generally, as input scAmpi_clinical expects the results of a DE analysis on a cluster level, and the files should follow the file name convention `SAMPLEID.CLUSTER.txt`, where SAMPLEID is the sample name specified in the sample map, CLUSTER is the cell cluster ID, and `txt` is the expected suffix.
The input files should be provided in the directory `results/parse_diff_exp/` and contain at least five mandatory columns:

```
gene_names  diff    padj      test_statistic  pct_nonzero
ATP1A1      1.679   3.05e-15  14.506          81.42
```

Here, "gene_names" contains the HGNC gene symbols, "diff" contains the fold change or a similar value, "padj" contains the adjusted p-value, "test_statistic" contains the value of the test statistics, and "pct_nonzero" contains the percentage of cells in this cluster with non-zero expression in the respective gene.
The results of this clinical pipeline run are the in-silico drug prediction and clinical annotations.
Other side results, e.g. the minimum set cover computation, the plotting of drug predictions on the UMAP, and the gene set enrichment analysis, cannot be created in an independent clinical run as they rely on input files generated by the scAmpi_basic part.

#### Adapting/Integrating rules in Snakemake

Snakemake is a Python-based workflow management system for building and executing pipelines. A pipeline is made up of ["rules"](snake/scAmpi_basic_rules.py) that represent single steps of the analysis. In a [yaml config file](config/config_scAmpi.yaml) parameters and rule-specific input can be adjusted to a new analysis without changing the rules. In a ["master" snake file](snake/snake_scAmpi_basic_master.snake) the desired end points of the analysis are specified. With the input and the desired output defined, Snakemake is able infer all steps that have to be performed in-between.

To change one of the steps, e.g. to a different software tool, one can create a new rule, insert a new code block into the config file, and include the input/output directory of this step in the master snake file. It is important to make sure that the format of the input and output of each rule is compatible with the previous and the subsequent rule. For more detailed information please have a look at the excellent [online documentation](https://snakemake.readthedocs.io/en/stable/index.html) of Snakemake.

#### Quick start using test data

To quickly start a scAmpi run with PBMC test data you can follow the following steps:

- clone the scAmpi repository
- make sure you have `snakemake` in your PATH (see Installation instructions)
- prepare Cellranger software and reference directory
- update the path to the cellranger reference directory in `testdata/config.yaml`
- download example data from the 10xGenomics website (for more detailed instructions see `testdata/README_testdata.md`)
- *optional*: to circumvent the time-consuming mapping step create the directory `results/counts_raw/` in your cloned scAmpi repository, copy the raw matrix `5k_pbmc_v3.h5.tar` into this directory, gunzip the file (e.g. `tar -xvf 5k_pbmc_v3.h5.tar`) and start the test run from this step.
- perform Snakemake dryrun to see a list of steps that will be performed
  `snakemake -s workflow/snakefile_basic.smk --configfile testdata/config.yaml -n -p`
- start analysis run
  `snakemake -s workflow/snakefile_basic.smk --configfile testdata/config.yaml -p`
