#_################################################
## File name: sample_integration.R
## Author: Anne Bertolini
## Date created: September 2019
## Last update: October 2020
## R Version: 3.5.1
##################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(plyr)
  library(dplyr)
  library(SingleCellExperiment)
  library(optparse)
  #library(cowplot)
  library(reticulate)
  library(ggplot2)
})

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# parse command line arguments
option_list = list(
make_option("--cohort_list", type = "character", help = "List with paths to the input RDS files that contain SCE objects."),
make_option("--sample_data", type = "character", help = "Path to the RDS file containg the SCE object of the current sample."),
make_option("--sampleName", type = "character", help = "Sample name that will be added to the names of all output files."),
make_option("--sampleName_short", type = "character", help = "Short sample name that is used to check if the current sample is already in the cohort. Sample name that is also used for cellranger."),
#make_option("--cluster_celltypes", type = "character", help = "Path to table with clusters and dominant cell type."),
make_option("--colour_config", type = "character", help = "Path to config file table that sorts colours to cell types"),
make_option("--outdir", type = "character", help = "Full path to output directory.")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# specify the python version to use
# This is the python version used in the TP pipeline
use_python("/cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/python")
print("Using python version: /cluster/work/bewi/ngs/projects/tumorProfiler/code/installations/snakemake_v5.1.4/bin/python when calling from R.")

# convenience function for string concatenation
'%&%' = function(a,b) paste(a,b,sep="")


########################
###   Read in data   ###
########################

# read in the data from all the samples, of the cohort and the current sample
cohort_files <- read.csv(opt$cohort_list, header = FALSE, stringsAsFactors = FALSE)
cohort_files <- cohort_files$V1

# check if current sample is already in cohort list and remove it
print("Current sample:")
print(opt$sampleName_short)
mask_current_sample <- grepl(opt$sampleName_short, cohort_files)
print("Checking if current sample is already in the cohort list. If yes, it will be removed to avoid duplication.")
print(table(mask_current_sample))
cohort_files <- cohort_files[!(mask_current_sample)]
print("cohort_files, current sample is removed if it was already in the cohort list:")
print(cohort_files)

# include current sample in list of samples
all_files <- c(cohort_files, opt$sample_data)
print("All files, including cohort and current sample:")
print(all_files)

# make a list of all SCE objects
all_sce_objects <- lapply(seq_along(all_files), function(x) {
  print(paste0("Read in sample ", x))
  print(all_files[x])
  my_sce <- readRDS(all_files[x])
  # sample name is the basename of the input file up to (but excluding) the second "."
  my_sce@metadata$sample_name <- regmatches(basename(all_files[x]), regexpr("[^\\.]+\\.[^\\.]+", basename(all_files[x])))
  print(my_sce@metadata$sample_name)
  my_sce
})

###################################
###   Generate Seurat objects   ###
###################################

# generate list of Seurat objects including all samples
all_samples <- lapply(seq_along(all_sce_objects), function(x) {
  print(paste0("Build Seurat Object of sample ", x, ". SCE object:"))
  print(all_sce_objects[x])
  seurat_obj <- CreateSeuratObject(counts = assay(x = all_sce_objects[[x]], "counts"), min.cells = 0, min.features = 0)
  seurat_obj[["phenograph_clusters"]] <- colData(all_sce_objects[[x]])$phenograph_clusters
  seurat_obj[["celltype_major"]] <- colData(all_sce_objects[[x]])$celltype_major
  seurat_obj[["celltype_final"]] <- colData(all_sce_objects[[x]])$celltype_final
  seurat_obj[["fractionMT"]] <- colData(all_sce_objects[[x]])$fractionMT
  seurat_obj[["n_umi"]] <- colData(all_sce_objects[[x]])$n_umi
  seurat_obj[["n_gene"]] <- colData(all_sce_objects[[x]])$n_gene
  seurat_obj[["log_umi"]] <- colData(all_sce_objects[[x]])$log_umi
  seurat_obj[["g2m_score"]] <- colData(all_sce_objects[[x]])$g2m_score
  seurat_obj[["s_score"]] <- colData(all_sce_objects[[x]])$s_score
  seurat_obj[["cycle_phase"]] <- colData(all_sce_objects[[x]])$cycle_phase
  seurat_obj[["celltype_major_full_ct_name"]] <- colData(all_sce_objects[[x]])$celltype_major_full_ct_name
  seurat_obj[["celltype_final_full_ct_name"]] <- colData(all_sce_objects[[x]])$celltype_final_full_ct_name
  seurat_obj[["sample_name"]] <- all_sce_objects[[x]]@metadata$sample_name
  seurat_obj[["phenograph_cluster"]] <- colData(all_sce_objects[[x]])$phenograph_clusters
  # mark current sample
  if (grepl(all_sce_objects[[x]]@metadata$sample_name, opt$sampleName)){
    seurat_obj[["current_sample"]] <- "current_sample"
  }else{
    seurat_obj[["current_sample"]] <- "cohort_sample"
  }
  # give out seurat object
  seurat_obj
})


##################################
###   Integration of samples   ###
##################################

# Perform SCTransform on each of the new Seurat objects, separately:
cat("\n\nPerform SCTransform.\n")
print(Sys.time())
for (i in seq_len(length(all_samples))) {
  # default for latent_var is c("log_umi")
  all_samples[[i]] <- SCTransform(all_samples[[i]], verbose = FALSE, latent_var = c("log_umi"), vars.to.regress = c("g2m_score", "s_score"))
}
# Select integration features and prepare integration
cat("\n\nPrepare integration.\n")
print(Sys.time())

# The number of features might have to be adapted per project, e.g. reduced if many samples are integrated or if no reference is used
# only those integration features are integrated
integration_features <- SelectIntegrationFeatures(object.list = all_samples, nfeatures = 2000)
all_samples <- PrepSCTIntegration(object.list = all_samples, anchor.features = integration_features)
# Find integration anchors
cat("\n\nFind integration anchors.\n")
print(Sys.time())

integration_anchors <- FindIntegrationAnchors(object.list = all_samples, normalization.method = "SCT",
                                  anchor.features = integration_features, dims = 1:20)
                                  # , reference = reference_datasets)
# Integrate samples
cat("\n\nIntegrate samples.\n")
print(Sys.time())
seurat_integrated <- IntegrateData(anchorset = integration_anchors, normalization.method = "SCT", dims = 1:20)

# run normalisation (scale), PCA and UMAP, find neighbours and clusters
print("\n\nRun normalisation (scale), PCA and UMAP, find neighbours and clusters.\n")
print(Sys.time())
DefaultAssay(object = seurat_integrated) <- "integrated"
seurat_integrated <- RunPCA(object = seurat_integrated, npcs = 40, verbose = TRUE)
seurat_integrated <- RunUMAP(object = seurat_integrated, reduction = "pca", dims = 1:40)
seurat_integrated <- FindNeighbors(object = seurat_integrated)
seurat_integrated <- FindClusters(seurat_integrated, dims.use = 1:40)
seurat_integrated[["barcodes"]] <- rownames(seurat_integrated[[]])

# Write integrated Seurat object into RDS file
filename_out <- opt$outdir %&% opt$sampleName %&% ".integrated_seurat_object.RDS"
print(paste0("File name of output RDS file with integrated Seurat object: ", filename_out))
saveRDS(seurat_integrated, filename_out)


###############################
###   Sample Visualisation  ###
###############################

# read in colour_config and sort colours to cell types
config <- read.csv(opt$colour_config, sep = "\t", stringsAsFactors = FALSE)
config$cell_type <- gsub(pattern = "([^_]+)_.*", "\\1", config$cell_type)
print("Check if the number of cell types and the number of colours for the cell types match.")
stopifnot(length(config$colour) == length(config$cell_type))
stopifnot(config$colour == unique(config$colour))
# have character vector of colour names from the colour config
ct.color = c(config$colour, "grey50", "black")
# give names to colours from other column of colour config
names(ct.color) <- c(config$cell_type, "uncertain", "unknown")
# make sure 'celltype_major' and 'celltype_final' metadata columns contain factors
seurat_integrated[["celltype_major"]] <- as.factor(seurat_integrated[[]]$celltype_major)
seurat_integrated[["celltype_final"]] <- as.factor(seurat_integrated[[]]$celltype_final)
# get integer vector containing the index of the matching colour in ct.color for each
# level of seurat_seurat_integrated[[]]$celltype_major
id.first.ct = match(levels(seurat_integrated[[]]$celltype_major), names(ct.color))
id.first.ct
id.final.ct = match(levels(seurat_integrated[[]]$celltype_final), names(ct.color))
id.final.ct


colstypes <- c("red2", "green4", "blue2", "cyan2", "yellow1", "purple", "brown",
               "chocolate1", "chartreuse2", "darkgoldenrod3", "steelblue1", "slateblue3", "olivedrab4", "gold2",
               "violetred3", "darkcyan","orchid3", "darksalmon", "darkslategrey", "khaki", "indianred2","magenta", "slategray2", 
               "olivedrab1", "mediumaquamarine", "hotpink", "yellow3",
               "bisque4", "darkseagreen1", "dodgerblue3",
               "deeppink4", "sienna4", "mediumorchid4")

p1 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "current_sample", cols = c("grey70", "red2"))
p1
ggsave(filename = paste0(opt$outdir, opt$sampleName, ".sample_integration_highlight_current.png"), width = 30, height = 20, dpi = 600, units = "cm")
print("Written plot 1")
#p2 <- DimPlot(seurat_integrated, reduction = "umap", label = TRUE)
#p2
#ggsave(filename = paste0(opt$outdir, "umap_plot_seurat_integrated.seurat_umap_clusters.png"), width = 30, height = 20, dpi = 600, units = "cm")
p3 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "sample_name")
p3
ggsave(filename = paste0(opt$outdir, opt$sampleName, ".sample_integration_all_samples.png"), width = 30, height = 20, dpi = 600, units = "cm")
print("Written plot 2")
p4 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "celltype_final", cols = ct.color[id.final.ct])
p4
ggsave(filename = paste0(opt$outdir, opt$sampleName, ".sample_integration_celltypes.png"), width = 35, height = 20, dpi = 600, units = "cm")
print("Written plot 3")
p5 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "celltype_major", cols = ct.color[id.first.ct])
p5
ggsave(filename = paste0(opt$outdir, opt$sampleName, ".sample_integration_major_celltypes.png"), width = 30, height = 20, dpi = 600, units = "cm")
print("Written plot 4")
p6 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "seurat_clusters", cols = colstypes)
p6
ggsave(filename = paste0(opt$outdir, opt$sampleName, ".sample_integration_seurat_umap_clusters.png"), width = 30, height = 20, dpi = 600, units = "cm")
print("Written plot 5")

