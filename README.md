# single-cell-tumor-profiler
[![CircleCI](https://circleci.com/gh/cbg-ethz/single-cell-tumor-profiler.svg?style=svg&circle-token=2f10c1d560456ff5407e27d3b688eccf9fcbca21)](https://circleci.com/gh/cbg-ethz/single-cell-tumor-profiler)
<br>

Analysis procedure:
1. Retrieve sample analysis information from LabKey (specifically, sample prefix). Use sample ID to create analysis folder (aka patient directory). 
2. copy all (empty) necessary subdirectories into the patient directory. the example directory used for copying can be found here:
`/cluster/work/bewi/ngs/projects/tumorProfiler/analysis/exampleFolder`
4. within the patient's directory, go into the folder 'openbis'. Then connect to the openbis server and mirror the fastq files for your sample.
5. with the patient's directory, go into the folder 'analysis/fastqc'. Then connect to the openbis server and mirror the fastqc files for your sample.
-> Alternatively, connect to openBIS and retrieve both, the fastq and fastqc files, and move files to respective target directory
7. load module files (see below)
9. run the script `prepare_files.py` to create and adapt necessary snakemake files:
    `/cluster/work/bewi/ngs/projects/tumorProfiler/code/all/scTumorProfiler/scripts/prepare_files.py`
-> example call:
`/cluster/work/bewi/ngs/projects/tumorProfiler/code/all/scTumorProfiler/scripts/prepare_files.py -s USZ-M-21.MSJLF.T-scR.6000c-r1-v1_0-AHYJ7NBGX5.r1-v1_0.39484 -p path/to/folder/singlecell_rna/ -c path/to/snake_analysis_files/config_scTranscriptomics.json -t melanoma -n MSJLF -f path/to/fastqc/folder/ `
10. go to `patientDirectory/snake_analysis_files` and perform the dry run with `./auto_dryrun_sctranscriptomics.sh`
11. perform the run with `./auto_run.sh`
12. For leomed upload, go to `patientDirectory/analysis/result_upload`
-> adapt content of the derived/SummaryFile and create new md5 for the Summary
-> rsync to respective folder on LeoMed (see below)

Upload:
`rsync -rtvP result_upload/\* leomed:/cluster/path/sample_folder/`
After sanity check: create done-file (touch done.txt).

* System requirements:

   * on Euler:
      - all r packages have been tested with r 3.5.1
      - required modules:
      `module load /cluster/project/nexus/utilities/sharedPrograms/cellranger/cellranger-3.1.0/cellranger.modulefile`
      `module load eth_proxy; module load new r/3.5.1`
      - all python-related files, including the snakemake installation, are installed in:
      `/cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/`


   * Python packages:
   `pip install -r python-requirements.txt`
   * R packages:
   `Rscript r-requirements.R`

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
1   melanomaSampleID1   sampleCellRangerName    N   1   /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/
1   melanomaSampleID2   sampleCellRangerName    N   1   /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_melanoma/
1   ovarianSampleID1    sampleCellRangerName    N   1   /cluster/work/bewi/ngs/projects/tumorProfiler/analysis/trial_ovarian/
...
```
Beware the mandatory '/' at the end of each root path. Also, the pipeline expects each sample's fasq subdirectory to use the cellranger sample name instead of the "full" sample name.

Example call:
`/cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/snakemake --notemp --latency-wait 60 -s /path/to/git/snake/snake_scTranscriptomics_novaSeq_master.snake --configfile /path/to/preprocessing/folder/snake_analysis_files/config_scTranscriptomics.json --cluster 'bsub -M {params.mem} -n {threads} -W {params.time} -R "rusage[mem={params.mem},scratch={params.scratch}]" -eo {params.lsferrfile} -oo {params.lsfoutfile}' -j 100 -p -k -n`





