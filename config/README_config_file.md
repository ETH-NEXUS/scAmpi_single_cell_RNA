# Help on config file content

## General

Before running the pipeline the `config` file needs to be adapted to contain the input and output paths for the intended analysis. Those are provided in the first section (`inputOutput`) of the config file.

    inputOutput:
        input_fastqs: "/dir/fastqs/"
        input_fastqc: "" # fastqc files, if already available
        analysis_output_dir: "/dir/output/"
        analysis_temp_dir: "/dir/snake_temp/"
        sample_map: "/dir/sample_map.tsv" # sample map with information on the samples of this analysis
        malignant_cell_type: "melanoma" # cell type name of malignant/diseased cell types - relevant to select clusters that will be analyzed in the clinical part of scAmpi. Note that cell types are selected via non-strict pattern matching, so also parts of diseased cell type labels can be used.

Apart from the fastq files, further resources must be provided in the section `resources`. Necessary for the cellranger step is the entry `reference_transcriptome`, other entries that should be adapted depending on the sample include `priority_genes` and `colour_config` for the plotting step, or `celltype_lists` and `celltype_config` for the cell type classification.
The scAmpi repository already includes resources for several tissue types, those resources can be used as examples to set up the resources for new tissue types.

    resources:
        pathwayDB: "../required_files/hallmark_pathways_converted_example.gmt"
        drugList: "../required_files/melanoma/melanoma_drug_list.txt" # only required for clinical part
        drugCombinations: "../required_files/melanoma/drug_combinations_melanoma.txt"     # only required for clinical part, can also be empty
        civicDict: "../required_files/drug_synonyms_civic.txt"    # only required for clinical part
        reference_transcriptome: "/resources/refdata-cellranger-GRCh38-3.0.0"  # e.g. downloaded from 10x Genomics resource 
        transcriptome_code: "GRCh38"
        celltype_lists: "../required_files/melanoma/celltype_list_melanoma.gmx"  # tab separated file with marker genes for each cell type we can expect in the respective tissue
        celltype_config: "../required_files/melanoma/celltype_config_melanoma.tsv"  # tab separated file that indicates which of the cell types in "celltype_lists" are major and which are subtypes
        colour_config: "../required_files/melanoma/colour_config_melanoma.txt"  # specify the colour for each cell type mentioned in "celltype_config"
        genesets: "../required_files/hallmark_pathways_example.gmt"  # required for GSVA analysis and clinical part
        priority_genes: "../required_files/melanoma/selected_genes_melanoma.txt"  # tab separated file that lists genes that shall be visualized on a UMAP. Genes are categorized into to user-defined gene-categories to allow better interpretation.


## Disease-specific parameters for CIViC

The scAmpi workflow leverages expert curated clinical data from the [CIViC](https://civicdb.org) database ("Clinical Interpretation of Variants in Cancer"). Variant specific evidence data retrieved from CIVIC can be further filtered for cancer-type specificity using parameters `--cancerTypeList`, `--blackList` and `--highLevelList`. Search terms should be provided in a comma-separated list (no spaces) and multiple words per term are allowed, eg. `ovarian`, `solid tumor,melanoma` and `sex cord-stromal,granulosa cell,glandular pattern` are all valid inputs.

Relevant and non-allowed disease names or terms can be provided in `--cancerTypeList` and `--blackList`, respectively. In both cases, partial matches are sought, eg. `small` will match `non-small cell lung cancer` and `lung small cell carcinoma`, while `non-small` will only match `non-small cell lung cancer`. In the same manner, be aware that `uveal melanoma` will only match `uveal melanoma` and not `melanoma`. As CIVICdb contains some higher-level disease names which are database specific, eg. `cancer` or `solid tumor`, terms provided in the `--highLevelList` only allow exact matches, eg. `cancer` will only match `cancer` and not `lung cancer`.

To select the evidences that will be reported, the following logic is applied based on their associated disease name:
* If any black-listed disease names are provided, partial matches will be sought and evidences associated to the matched diseases excluded from the query.
* Then partial matches to the white-listed terms will be sought, and if any are matched, only associated evidences will be reported using cancer-specific tag `ct`.
* If white-listed terms are not provided or not matched, then exact matches to the high-level terms will be considered as a fall-back case. Only evidences associated to the matched diseases will be reported using general tag `gt`.
* If high-level terms are not provided or not matched, then all available evidences will be reported regardless of the associated disease (except those black-listed, if any) using non-specific tag `nct`.

The above logic is applied separately for each evidence type (one of `Predictive`, `Diagnostic`, `Prognostic` or `Predisposing`). This means that evidences classified as eg. `Predictive` might be associated to a different set of disease names compared to those classified as eg. `Diagnostic`.

To ease the selection of the appropriate parameter terms for a particular disease, we provide a helper file `civic_available_diseases_[DATE].txt` listing all disease names available in CIVICdb as of `[DATE]`. To update this file, run the following script as follows, replacing `[DATE]` withe new date:

    > python /path_to_git_scAmpi/scripts/get_available_diseases_in_civic.py --outFile /path_to_git_scAmpi/required_files/civic_available_diseases_[DATE].txt


