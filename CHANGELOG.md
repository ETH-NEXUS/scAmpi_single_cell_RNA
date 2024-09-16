# Changelog

## [2.2.0] - 2024-08-12

### Changed
- Added automated fastq file link creation for cell ranger
  If an additional column `file_stem` is provided in the sample map, 
  the fastq files starting with this file stem are assumed to correspond to one sample and links to Cell Ranger compatible folder structure are created automatically. 
  Otherwise, if there is just the sample column, `fastq_dir` is expected to point to Cell Ranger compatible folder structure. 

- Added option to specify `samples` in the config, e.g. to perform a testrun on a one or two samples, or to rerun specific samples. 

- Using containers in the Cell Ranger count step. 
  Cell Ranger version is specified in the config. Currently supported versions 
  include "7.2.0" and "8.0.1".

- added SeaCells und Metacells2 metacell approaches, using the dedicated snakefiles 
  snakefile_seacells.smk or snakefile_metacells2.smk, respectively. These produce
  the same results as on single cells, but in addition all results also on metacells.
  Further, celltyping is compared between metacells and single cells, and report files are generated.

- added files to generate basic snakemake reporting, including rule specific reporting on results. 
  At the current state, this feature is prove of concept only. 
  
- made cell cycle regression optional in sctransform, default is no correction. 
  Can be switched on with --cell_cycle_correction command line argument. 

- added error handling to sctransform::vst. This prevents the script from crashing, 
  as well as reports genes where the algorithm did not converge, which previously
  gave an unspecific warning.

- added some logging. 
## [2.1.1] - 2024-08-15

### Changed

- adapt example cell type classification genes for ovarian cancer and melanoma samples.


## [2.1.0] - 2024-07-30

### Changed

- update cellranger rules  
  have new `cellranger_count_8` rule that includes syntax changes of cellranger v8. The new rule is the default; if an older version of cellranger should be used with the rule `cellranger_count` the  
  `ruleorder: cellranger_count > cellranger_count_8`  
  in the Snakefile must be adapted.  
  Also, add new rule `gunzip_and_link_cellranger` and separate these steps from the cellranger rule.

- update and describe in main README how cellranger expects to find the raw FASTQ files.

- update conda environments

  - `celltyping.yaml`
  - `identify_doublets.yaml`
  - `sctransform_preprocessing.yaml`

- update `identify_doublets` output

  - Have results of `identify_doublets` rule in own subdirectory instead of the `counts_filtered` directory.

- update `generate_qc_plots_*`

  - have resulting QC plots written in own subdirectory instead of same directory as count files.
  - Improve memory usage.
  - Clean up script.

- update `filter_genes_and_cells.R`

  - implement iterative filtering to make sure the selected thresholds for genes and cells apply to all genes/cells of the downstream analyses
  - Clean up script.

- update `plotting.R`

  - add more colours for clusters. Make sure even with a high number of clusters, enough colours are provided.
  - make sure all cell types that are not found in a sample are still shown in the legend (with `show.legend = T`, adapt to new ggplot2 default settings)
  - Clean up script.

- update `create_hdf5.py`
  - make sure the script can work with Human and also Mouse data. Mouse Ensembl gene IDs are longer than 16 characters, and cannot be of type `dtype='S16'`.

### Fixed

- fix `sctransform_preprocessing.R`

  - Filtering of raw input files is not applied to row and column names. This issue should have had no effect as long as filtered input data was provided (with minimum of QC on genes and cells).  
  - Changed to a check that stops the script if unfiltered input is detected.
  - Script linting.


## [2.0.7] - 2023-03-20

### Changed

- specify which library should be used for the function `ggsave` to avoid conflict between the R packages `ggplot2` and `cowplot`

## [2.0.6] - 2023-02-14

### Fixed

- adapt script `query_civic_expr.py` to changed syntax in python package `civicpy` version 3.0. The script no longer works as is with previous versions of the package.
- adapt installation instructions for `civicpy` to require version 3.0

## [2.0.5] - 2023-01-11

### Fixed

- adapt script `query_civic_expr.py` to changed syntax in python package `civicpy` version 2.0. The script no longer works as is with previous versions of the package.
- adapt installation instructions for `civicpy` to require version 2.0

## [2.0.4] - 2022-12-01

### Fixed

- adapt R scripts using ggplot2 for plotting (`filter_genes_and_cells.R`) to changes introduced with ggplot v3.4.0. With new defaults in `scale_*_manual` only factor levels found in the data are shown in the legend.
- in rule `cellranger_count` use consistently full path given in `config[inputOutput][input_fastqs]`. Update respective README sections.
- fix relative path to conda env yaml file in rule inheritance of `create_hdf5_starsolo`.

## [2.0.3] - 2022-11-21

### Fixed

- fix bug in script `filter_genes_and_cells.R` that resulted in colour discrepancy between legend and plot in `{sample}.visualize_filtered_cells.png` in rare cases.

## [2.0.3] - 2023-03-20

### Changed

- specify which library should be used for the function ggsave to avoid conflict between ggplot2 and cowplot

## [2.0.2] - 2022-08-31

### Changed

- change run time keyword in rules and config from "time_min" or "time" to "runtime" as is recommended.
- have template memory values per job and not per thread
- use "mem_mb" consistently in rule files and config files

## [2.0.1] - 2022-07-19

### Changed

- Update example Ovarian cancer cell types to include Mesothelial cells and pDCs

### Fixed

- Fix empty output plot from rule `plot_upsetr`

## [2.0.0] - 2022-07-05

The scAmpi pipeline framework was fully revised to follow the current Snakemake best practices.
No changes were made regarding the analysis steps.

### Changed

- the Snakemake pipeline framework was updated to follow current best practices
- using Snakemake checkpoints the scAmpi basic and clinical parts were joined together to be one workflow
- several samples can be run in parallel (previously that was only possible for scAmpi_basic)

## [1.0.0] - 2022-04-05

First full, publicly available version of the scAmpi single-cell RNA analysis pipeline.
