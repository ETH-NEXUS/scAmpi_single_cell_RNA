# single_cell_tumor_profiler

#### General overview

This pipeline consists of two main parts: the `sc-transcriptomics` part and the `clinical` part. The first comprises general scRNA preprocessing steps, filtering, normalisation, unsupervised clustering, cell type classification, and DE analysis.

The latter (`clinical`) part includes the search for disease relevant drug targets of differentially expressed genes when comparing the malignant cells to the non-malignant cells of a sample.

#### Software

The installation of Snakemake is prerequisite for this workflow.

The pipeline consists of R and python scripts. Most of the relevant R and python packages used in this workflow can be installed as a conda environment using the provided yaml file `tp_scRNA_conda_env.yml`.

Example:
```
conda env create -f tp_scRNA_conda_env.yml --name tp_scRNA
```

For the normalisation step the R package `sctransform` is used. As there is a known bug in the latest release either an older package version (e.g. 0.2) or the development version should be installed. The development version can be installed within R using `remotes::install_github("ChristophH/sctransform@develop")`.

Phenograph is used for the unsupervised clustering step. It can be installed using `pip install PhenoGraph` (with the conda environment activated).

#### Example data

For a test run the freely available 10X Genomics data from PBMC cells can be used. Please find an example config file and in the directory `config`


#### Before running the pipeline

Before running the pipeline the `config` file needs to be adapted to contain the input and output paths. Those are provided in the first section (`inputOutput`) of the config file.

Apart from the fastq files, further resources must be provided in the section `resources`. Necessary for the cellranger step is the entry `reference_transcriptome`, other entries that should be adapted depending on the sample include `priority_genes` and `colour_config` for the plotting step, or `celltype_lists` and `celltype_config` for the cell type classification.

Also, a "sample_map" must be provided, a tab delimited text file. Here, all samples are listed and each sample is a line in the sample map.

The sample map contains four columns, most interesting is the second column with the sample name that will be used to name files and identify the sample throughout the pipeline.

Example:
```
1	SAMPLE-1_scR	N	1
1	SAMPLE-2_scR	N	1
```



#### sc-transcriptomics part

Example call:

```
snakemake --notemp --latency-wait 60 -s snake_scTranscriptomics_master.snake --configfile config_scTranscriptomics.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k
```



#### Clinical part

Example call:

```
snakemake --notemp --latency-wait 60 -s snake_scTranscriptomics_master.snake --configfile config_scTranscriptomics.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k
```

