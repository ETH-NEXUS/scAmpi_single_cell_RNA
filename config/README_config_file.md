# Help on config file content

Before running the pipeline the `config` file needs to be adapted to contain the input and output paths for the intended analysis. Those are provided in the first section (`inputOutput`) of the config file.

    "inputOutput":{
        "input_fastqs":"/dir/fastqs/",
        "input_fastqc":" ", # fastqc files, if already available
        "analysis_output_dir":"/dir/output/",
        "analysis_temp_dir":"/dir/snake_temp/",
        "sample_map":"/dir/sample_map.tsv", # sample map with information on the samples of this analysis
        "malignant_cell_type":"melanoma", # cell type name of malignant/diseased cell types - relevant to select clusters that will be analyzed in the clinical part of scAmpi. Note that cell types are selected via non-strict pattern matching, so also parts of diseased cell type labels can be used.
    },

Apart from the fastq files, further resources must be provided in the section `resources`. Necessary for the cellranger step is the entry `reference_transcriptome`, other entries that should be adapted depending on the sample include `priority_genes` and `colour_config` for the plotting step, or `celltype_lists` and `celltype_config` for the cell type classification.
The scAmpi repository already includes resources for several tissue types, those resources can be used as examples to set up the resources for new tissue types.

    "resources":{
        "pathwayDB":"../required_files/hallmark_human_MAPK_and_filtered_converted.gmt",
        "drugList":"../required_files/melanoma_drug_list_converted_corrected.txt", # only required for clinical part
        "drugCombinations":"../required_files/drug_combinations_melanoma.txt",     # only required for clinical part, can also be empty
        "civicDict":"../required_files/drug_synonyms_civic.txt",    # only required for clinical part
        "reference_transcriptome":"/resources/refdata-cellranger-GRCh38-3.0.0",  # e.g. downloaded from 10x Genomics resource 
        "transcriptome_code":"GRCh38",
        "celltype_lists":"../required_files/celltype_list_melanoma.gmx",  # tab separated file with marker genes for each cell type we can expect in the respective tissue
        "celltype_config":"./celltype_config_melanoma.tsv",  # tab separated file that indicates which of the cell types in "celltype_lists" are major and which are sub types
        "colour_config":"./colour_config_melanoma.txt",  # specify the colour for each cell type mentioned in "celltype_config"
        "genesets":"../required_files/hallmark_human_MAPK_and_filtered.gmt",  # only required for clinical part
        "priority_genes":"../required_files/selected_genes_melanoma.txt"  # tab separated file that lists genes that shall be visualized on a UMAP. Genes are categorized into to user-defined gene-categories to allow better interpretation.
    },


