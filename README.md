# single-cell-tumor-profiler

* System requirements:

   *  - all r packages have been tested with r 3.5.1
   * Python packages:
   `pip install -r python-requirements.txt`
   * R packages:
   `Rscript r-requirements.R`
   * cellranger needs to be installed (tested with cellranger-3.1.0 and earlier versions)

### Comment on sample integration

The rule `sample_integration` performs an integration of the current sample and a cohort of previous samples using the Seurat 3 framework.
This is done to visualise the cells of all samples together in one UMAP.
Because it is error prone and would use a great amount of resources, the cohort of previous samples does not contain all previously processed samples of the respective indication.
Instead, we selected samples that were of high quality from a technical point of  view and consist of tumour as well as immune cells. The first aspect makes sure that enough common features (genes) are found in the samples.
The latter is necessary for the integration software to find 'similar' cells of the same cell type in different samples that are used as 'anchors'.
Even though the cohort should not exceed a for the purpose reasonable number of samples it will be updated as we get more samples that meet the criteria mentioned above.


## Preprocessing with NovaSeq

When several samples (not necessarily from the same indication) were sequenced together using HiSeq 4000 (ie. NovaSeq), a preprocessing pipeline needs to be applied before the actual analysis can be performed. This pipeline makes use of `snake_scTranscriptomics_novaSeq_master.snake` and `config_scTranscriptomics_novaSeq.json`, and it expects a different directory structure compared to the regular analysis run: a single "root" fastq directory is used which contains one subdirectory per sequenced sample.

Furthermore, in this pipeline the sample map is used to specify sample-specific parameters: the cellranger sample name (ie. short version containing no dots) corresponding to each sample, and the root directory for each sample's corresponding indication (as results from the preprocessing steps will be directly copied into each sample-specific directory for the analysis). Hence, the sample map should have the following example structure:
```
1   melanomaSampleID1   sampleCellRangerName    N   1   /path/to/melanoma/samples/
1   melanomaSampleID2   sampleCellRangerName    N   1   /path/to/melanoma/samples/
1   ovarianSampleID1    sampleCellRangerName    N   1   /path/to/ovarian/samples/
...
```
Beware the mandatory '/' at the end of each root path. Also, the pipeline expects each sample's fasq subdirectory to use the cellranger sample name instead of the "full" sample name.

Example call:
`snakemake --notemp --latency-wait 60 -s /path/to/git/snake/snake_scTranscriptomics_novaSeq_master.snake --configfile /path/to/preprocessing/folder/snake_analysis_files/config_scTranscriptomics.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k -n`





