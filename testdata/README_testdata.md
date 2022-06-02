# Guidelines to run testdata with scAmpi

The following instructions provide a short guideline to download example data and execute the scAmpi basic scRNA-seq analysis workflow. 

Note that the testdata are human PBMCs, available for download from the 10xGenomics website with `wget https://cg.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar`, and you need to have a version of the cellranger human reference data files available. (For download instructions of cellranger, refer to the [10xGenomics Webpage](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation)).
In the following we assume that the reference files are located in `reference_dir`.

#### (i) Make sure `snakemake` and `cellranger` are available in PATH

For instructions on how to install the necessary software to run scAmpi, please refer to the README.md file in the scAmpi git repository.

In order to do a test run you need to have `snakemake` and `cellranger` in your PATH. Note: if you do not want to add cellranger to your PATH, you can also provide the full path to the cellranger installation in the config file at the cellranger_count config block.

```
>> export PATH=/cellranger_dir/cellranger-6.0.0:$PATH
```

NOTE: some steps of the scAmpi workflow perform online queries. Please make sure that this is possible on your computing system, e.g. by loading the respective modules to enable the proxy connection. (Most systems will have this enabled per default).

#### (ii) Download example data from 10xGenomics
Make sure you are in the scampi directory.
```
mkdir fastqs
cd fastqs
wget https://cg.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_v3/5k_pbmc_v3_fastqs.tar
tar -xvf 5k_pbmc_v3_fastqs.tar
rm 5k_pbmc_v3_fastqs.tar
mv 5k_pbmc_v3_fastqs/* .
rm -r 5k_pbmc_v3_fastqs/
```
#### (iii) Adapt path to cellranger resources
In the testdata/config.yaml file replace the `resources:` -> `reference_transcriptome` entry with your path to the cellranger reference directory.


#### (iiii) Run scAmpi

Perform a dry run to check if all preparations were correct:
```
>> snakemake -s workflow/snakefile_basic.smk --configfile testdata/config.yaml -n -p
```

Call the actual run:

```
>> snakemake -s workflow/snakefile_basic.smk --configfile testdata/config.yaml -p
```

If you are working on a cluster with a job scheduling system (e.g. LSF, Slurm) you need to adjust the commands accordingly.

```
>> ./run_scAmpi.sh
```

#### Avoid cellranger step
If the cellranger step should be avoided in the test run, the file `5k_pbmc_v3.h5.tar` can be used as a starting point. With `tar -xvf 5k_pbmc_v3.h5.tar` the file can be unpacked. Then, it should be copied into the direcory `results/counts_raw/`. This way, scAmpi_basic starts after the cellranger step, saving time and other resources.
