# scAmpi - Single Cell Analysis mRNA pipeline

#### General overview

This scAmpi workflow is organized into two main parts: the `scAmpi_basic` part and the `scAmpi_clinical` part, which can be run independently. scAmpi_basic includes general scRNA processing steps, such as mapping, QC, normalisation, unsupervised clustering, cell type classification, and DE analysis.

scAmpi_clinical includes the search for disease relevant drug targets for differentially expressed genes. Note that the clinical part is only applied if at least one cluster identified in your sample is indicated as a diseased ("malignant") cell type.


![scAmpi_both_rulegraphs](https://user-images.githubusercontent.com/38692323/140029020-6292b989-722d-4c93-909d-1d65c8aacddd.png)



#### Installation instructions

scAmpi provides a yml file to enable installing most software used in the default workflow as a conda environment. After cloning this repository you find the yaml file `scAmpi_scRNA_conda_env.yml` in the subfolder "envs". Additionally, you find in this directory one concrete example of the conda environment used to run scAmpi with detailed package versions, `scAmpi_scRNA_detailed.yml`.

Example (This will install the conda environment in your home. Note that the conda may take a while to install all software packages):
```
>> conda env create -f scAmpi_scRNA_conda_env.yml --name scAmpi_scRNA
```

To activate the environment, type:
```
>> conda activate scAmpi_scRNA
```

Additionally required installations that are not available via conda:
- [Cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger): Follow the instructions on the 10xGenomics installation support page to install cellranger and to include the cellranger binary to your path.
Webpage: [https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)

#### Example data

For a test run the freely available 10X Genomics data from PBMC cells can be used. A step by step guideline and example config file are provided in the directory `testdata`.  

#### Before running the pipeline

Note: some steps of the scAmpi workflow perform online queries. Please make sure that this is possible on your computing system, e.g. by loading the respective modules to enable the proxy connection. (Most systems will have this enabled per default).  
Also note that currently the filtering step of scAmpi_basic is designed to handle human data only.

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
snakemake -s snake_scAmpi_basic_master.snake --configfile config_scAmpi.yaml -j 1 -p -k
```

Note that if the pipeline is run on a compute cluster with a job scheduling system the parameter `--cluster` can be added to the snakemake command, indicating memory and time resources for a job (e.g. `--cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}'` for lsf).


#### scAmpi_clinical part

Example call:

```
snakemake -s snake_scAmpi_clinical_master.snake --configfile config_scAmpi.yaml -j 1 -p -k
```

#### A note on using CIViC

The CIViC query implemented in scAmpi makes use of an offline cache file of the CIViC database. The cache is retrieved with the initial installation of the scAmpi software. Afterwards, users have to manually update the cache file if they want to use a new version. 
To update the cache file, load the scAmpi environemnt and open an R session.
Then type:
```
>> civic.update_cache()
```

#### Running scAmpi_clinical independently

It is possible to run the scAmpi_clinical part independently of scAmpi_basic, following some restrictions to the file names and formatting. For this use case please use the master snake file `snake_scAmpi_clinical-only_master.snake`. Generally, as input scAmpi_clinical expects the results of a DE analysis on a cluster level, and the files should follow the file name convention `SAMPLEID.CLUSTER.txt`, where SAMPLEID is the sample name specified in the sample map, CLUSTER is the cell cluster ID, and `txt` is the expected suffix.
The input files should be provided in the input directory specified in the config file ("input_fastqs") and contain at least five mandatory columns:
```
gene_names  diff    padj      test_statistic  pct_nonzero
ATP1A1      1.679   3.05e-15  14.506          81.42
```
Here, "gene_names" contains the HGNC gene symbols, "diff" contains the fold change or a similar value, "padj" contains the adjusted p-value, "test_statistic" contains the value of the test statistics, and "pct_nonzero" contains the percentage of cells in this cluster with non-zero expression in the respective gene. 
The results of this clinical pipeline run are the in-silico drug prediction and clinical annotations.  
Other side results, e.g. the minimum set cover computation, the plotting of drug predictions on the UMAP, and the gene set enrichment analysis, cannot be created in an independent clinical run as they rely on additional input files generated by the scAmpi_basic part.

#### Adapting/Integrating rules in Snakemake
Snakemake is a Python-based workflow management system for building and executing pipelines. A pipeline is made up of ["rules"](snake/scAmpi_basic_rules.py) that represent single steps of the analysis. In a [yaml config file](config/config_scAmpi.yaml) parameters and rule-specific input can be adjusted to a new analysis without changing the rules. In a ["master" snake file](snake/snake_scAmpi_basic_master.snake) the desired end points of the analysis are specified. With the input and the desired output defined, Snakemake is able infer all steps that have to be performed in-between.

To change one of the steps, e.g. to a different software tool, one can create a new rule, insert a new code block into the config file, and include the input/output directory of this step in the master snake file. It is important to make sure that the format of the input and output of each rule is compatible with the previous and the subsequent rule. For more detailed information please have a look at the excellent [online documentation](https://snakemake.readthedocs.io/en/stable/index.html) of Snakemake.
