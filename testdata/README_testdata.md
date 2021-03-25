# Guidelines to run testdata with scAmpi

The following instructions provide a short guideline to download example data and execute the scAmpi baic scRNA-seq analysis workflow. 

In the following, we refer to your test directory as `testdir` and to the directory of the scAmpi git repository as `path_git_scAmpi`. Note that the testdata are human PBMCs, so you need to have a version of the cellranger human reference data files available. (For download instructions from cellranger, refer to the [10xGenomics Webpage](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation).
In the following we assume that the reference files are located in `reference_dir`.

(i) Activate environment and export cellranger path

For instructions on how to install the necessary software to run scAmpi, please refer to the README.md file in the scAmpi git repository.

In order to do a test run, you need to activate the conda environment and you need to export the path to the cellranger installation to you path. The following example calls assume that the conda environemt was installed in you HOME idrectory and that cellranger (version 6.0.0) was installed in "cellranger_dir". Note: if you do not want to add cellranger to your PATH, you can also provide the full path to the cellranger installation in the config file at the cellranger_run config block.

```
> conda avtivate scAmpi_scRNA
> export PATH=/cellranger_dir/cellranger-6.0.0:$PATH
```

Note: some steps of the scAmpi workflow perform online queries. Please make sure that this is possible on you computing system, e.g. by loading the respective modules to enable the proxy connection. (Most systems will have this enabled per default).

(ii) Run script to prepare test directory
The helper script `prepare_testrun.sh` will create the fastqs and analysis directories, as well as downloads the testdata, and prepares the config file for the test run:

```
> cd testdir
> sh [path_git_scAmpi]/testdata/prepare_testrun.sh [path_git_scAmpi] [reference_dir]
```

(iii) Run scAmpi

The helper script created a sub folder called "snake_files" in your test directory and prepared the commands for a dry run and full run of the workflow.

```
> cd snake_files/
```

Call the dry run to check if all preparations were correct:

```
> ./dryrun_scAmpi.sh
```

Call the actual run:

```
> ./run_scAmpi.sh
```


