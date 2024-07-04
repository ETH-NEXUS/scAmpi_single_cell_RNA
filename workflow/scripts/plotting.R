##################################################
#### File name: plotting.R
#### Author: Anne Bertolini
#### Date created: June 2019
#### R Version: 4.0
##################################################

suppressPackageStartupMessages({
  library(optparse)
  library(plyr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(egg)
  library(ggridges)
  library(RColorBrewer)
  library(limma)
  library(cowplot)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

# parse command line arguments
option_list <- list(
  make_option("--sce_in", type = "character", help = "Path to the input RDS file that contains an SCE object."),
  make_option("--genelist", type = "character", help = "Path to a tab separated file with clinically relevant genes in the first column. The expression of those genes will be specifically plotted."),
  make_option("--colour_config", type = "character", help = "Text file that specifies a fixed colour for each cell type."),
  make_option("--sampleName", type = "character", help = "Sample name that will be added to the names of all output files."),
  make_option("--outDir", type = "character", help = "Full path to output directory."),
  make_option("--toggle_label", type = "logical", action = "store_true", default = TRUE, help = "Set gene plot title labels to include user-defined gene aliases. (On by default.)")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# give out session Info
cat("\n\n")
print(Sys.time())
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n")
cat("\nInput files:\n\n")
print(opt)
cat("\n\n")



# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")

# define function
f.round.preserve.sum <- function(x, digits = 0) {
  up <- 10^digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

# set path and sample prefix for all output files
path <- opt$outDir %&% opt$sampleName
# create subfolder for gene expression plots
dir_gene_expr <- opt$outDir %&% "gene_expression/"
path_gene_expr <- opt$outDir %&% "gene_expression/" %&% opt$sampleName
dir.create(dir_gene_expr, showWarnings = FALSE)
# create subfolder for qc plots
dir_qc <- opt$outDir %&% "qc_visualization/"
path_qc <- opt$outDir %&% "qc_visualization/" %&% opt$sampleName
dir.create(dir_qc, showWarnings = FALSE)
print(path)
print(path_gene_expr)
print(path_qc)


# set general options/parameters
options(stringsAsFactors = FALSE)
fontsize <- theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11))
theme_set(theme_bw(12) + fontsize)
col.pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)

# read in sce object from RDS file
my_sce <- readRDS(opt$sce_in)
print("my_sce:")
print(my_sce)



# Match colours to cell types. all_celltypes are all the cell types that are offered to the cells via the config file,
# major and subtypes
cell.types <- metadata(my_sce)$all_celltypes
all.cell.types <- c(sort(cell.types), "uncertain", "unknown")
cols33 <- c(
  "red2", "green4", "blue2", "cyan2", "yellow1", "purple", "brown",
  "chocolate1", "chartreuse2", "darkgoldenrod3", "steelblue1", "slateblue3", "olivedrab4", "gold2",
  "violetred3", "darkcyan", "orchid3", "darksalmon", "darkslategrey", "khaki", "indianred2", "magenta", "slategray2",
  "olivedrab1", "mediumaquamarine", "hotpink", "yellow3",
  "bisque4", "darkseagreen1", "dodgerblue3",
  "deeppink4", "sienna4", "mediumorchid4"
)

# get UMAP coordinates into data frame
umap_coord <- as.data.frame(reducedDims(my_sce)$umap_hvg)
names(umap_coord) <- c("umap1", "umap2")
umap_coord$barcodes <- rownames(umap_coord)

# make sure all cell dimensions are the same
stopifnot(my_sce$barcodes == colData(my_sce)$barcodes)
stopifnot(my_sce$barcodes == umap_coord$barcodes)
# add umap coordinates to the meta data table of cells
cell_attributes <- plyr::join(as.data.frame(colData(my_sce)), umap_coord)

# have cell types ordered alphabetically
levels_alpha <- sort(levels(cell_attributes$celltype_final), decreasing = FALSE)
print("levels_alpha:")
print(levels_alpha)
cell_attributes$celltype_final <- factor(cell_attributes$celltype_final,
  levels = levels_alpha
)

cell_attributes$phenograph_clusters <- factor(cell_attributes$phenograph_clusters)
print("head(cell_attributes):")
print(head(cell_attributes))
print("str(cell_attributes):")
print(str(cell_attributes))

# read in cell type colours
config <- read.csv(opt$colour_config, sep = "\t", stringsAsFactors = FALSE)
config$cell_type <- gsub(pattern = "([^_]+)_.*", "\\1", config$cell_type)
stopifnot(length(config$colour) == length(config$cell_type))
stopifnot(config$colour == unique(config$colour))
ct.color <- c(config$colour, "grey50", "black")
names(ct.color) <- c(config$cell_type, "uncertain", "unknown")
id.first.ct <- match(levels(cell_attributes$celltype_major), names(ct.color))
id.first.ct
id.final.ct <- match(levels(cell_attributes$celltype_final), names(ct.color))
id.final.ct


# Make width of plot dependent of the number of columns in the legend
# second cell type
print("Number of second cell types:")
print(length(all.cell.types))
if (length(all.cell.types) > 18) {
  width_second <- 25
} else if (length(all.cell.types) > 36) {
  width_second <- 30
} else {
  width_second <- 20
}
# first cell type
print("Number of first cell types:")
print(length(levels(cell_attributes$celltype_major)))
if (length(levels(cell_attributes$celltype_major)) > 16) {
  width_first <- 25
} else if (length(levels(cell_attributes$celltype_major)) > 32) {
  width_first <- 30
} else {
  width_first <- 20
}
# Number phenograph clusters
print("Number of phenograph clusters found:")
print(length(unique(cell_attributes$phenograph_clusters)))
if (length(unique(cell_attributes$phenograph_clusters)) > 16) {
  width_pheno <- 25
} else if (length(unique(cell_attributes$phenograph_clusters)) > 32) {
  width_pheno <- 30
} else {
  width_pheno <- 20
}

# plot UMAP (based on highly variable genes), colours = major cell type
p_tirosh_ct <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = celltype_major)) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 15), nrow = 16)) +
  scale_color_manual(name = "Cell type", values = ct.color[id.first.ct], drop = F) +
  coord_fixed(ratio = 1)
p_tirosh_ct_fix <- set_panel_size(p_tirosh_ct, width = unit(12, "cm"), height = unit(12, "cm"))
# p_tirosh_ct
ggsave(path %&% ".first_celltype.png", p_tirosh_ct_fix,
  width = width_first, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = final cell type
p_final_ct <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = celltype_final)) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 2, shape = 15), nrow = 18)) +
  scale_color_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  coord_fixed(ratio = 1)
p_final_ct_fix <- set_panel_size(p_final_ct, width = unit(12, "cm"), height = unit(12, "cm"))
# p_final_ct
ggsave(path %&% ".second_celltype.png", p_final_ct_fix,
  width = width_second, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = phenograph clusters
p_umap_pc <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = as.factor(phenograph_clusters))) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 2, shape = 15), nrow = 16)) +
  scale_color_manual(name = "Phenograph", values = cols33) +
  ggtitle(paste("Phenograph clustering. Modularity score: ", metadata(my_sce)$modularity_score)) +
  theme(plot.title = element_text(size = 12)) +
  coord_fixed(ratio = 1)
p_umap_pc_fix <- set_panel_size(p_umap_pc, width = unit(12, "cm"), height = unit(12, "cm"))
# p_umap_pc
ggsave(path %&% ".phenograph.png", p_umap_pc_fix,
  width = width_pheno, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = number of detected genes
p_umap_ng <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = n_gene)) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_colour_distiller(name = "Nr genes", palette = "RdBu") +
  ggtitle("Number of genes per cell") +
  theme(plot.title = element_text(size = 12)) +
  coord_fixed(ratio = 1)
p_umap_ng_fix <- set_panel_size(p_umap_ng, width = unit(12, "cm"), height = unit(12, "cm"))
# p_umap_ng_fix
ggsave(path_qc %&% ".NrGenes.png", p_umap_ng_fix,
  width = 20, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = library size (number UMI counts per cell)
cell_attributes$log10_umi <- log10(cell_attributes$n_umi)
p_umap_ls <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = log10_umi)) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_colour_distiller(name = "log10 of library size", palette = "RdBu") +
  ggtitle("Log10 of library size (number of UMI) per cell") +
  theme(plot.title = element_text(size = 12)) +
  coord_fixed(ratio = 1)
p_umap_ls_fix <- set_panel_size(p_umap_ls, width = unit(12, "cm"), height = unit(12, "cm"))
# p_umap_ls
ggsave(path_qc %&% ".LibSize_log.png", p_umap_ls_fix,
  width = 20, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = log of library size (number UMI counts per cell)
p_umap_lsn <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = n_umi)) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_colour_distiller(name = "Library size", palette = "RdBu") +
  ggtitle("Library size (number of UMI) per cell") +
  theme(plot.title = element_text(size = 12)) +
  coord_fixed(ratio = 1)
p_umap_lsn_fix <- set_panel_size(p_umap_lsn, width = unit(12, "cm"), height = unit(12, "cm"))
# p_umap_lsn
ggsave(path_qc %&% ".LibSize.png", p_umap_lsn_fix,
  width = 20, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = cell cycle phase
p_umap_ccphase <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = cycle_phase)) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_color_manual(name = "Cell cycle phase", values = cols33) +
  ggtitle("Cell cycle phase") +
  theme(plot.title = element_text(size = 12)) +
  coord_fixed(ratio = 1)
p_umap_ccphase_fix <- set_panel_size(p_umap_ccphase, width = unit(12, "cm"), height = unit(12, "cm"))
# p_umap_ccphase
ggsave(path_qc %&% ".ccphase.png", p_umap_ccphase_fix,
  width = 20, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = s_score)
p_umap_sscore <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = s_score)) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_colour_distiller(name = "s_score", palette = "RdBu") +
  ggtitle("s_score value per cell") +
  theme(plot.title = element_text(size = 12)) +
  coord_fixed(ratio = 1)
p_umap_sscore_fix <- set_panel_size(p_umap_sscore, width = unit(12, "cm"), height = unit(12, "cm"))
# p_umap_sscore
ggsave(path_qc %&% ".s_score.png", p_umap_sscore_fix,
  width = 20, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = g2m_score)
p_umap_g2m <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = g2m_score)) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_colour_distiller(name = "g2m_score", palette = "RdBu") +
  ggtitle("g2m_score value per cell") +
  theme(plot.title = element_text(size = 12)) +
  coord_fixed(ratio = 1)
p_umap_g2m_fix <- set_panel_size(p_umap_g2m, width = unit(12, "cm"), height = unit(12, "cm"))
# p_umap_g2m
ggsave(path_qc %&% ".g2m_score.png", p_umap_g2m_fix,
  width = 20, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = fraction MT genes)
p_umap_MT <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = fractionMT)) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  scale_colour_distiller(name = "fractionMT", palette = "RdBu") +
  ggtitle("Fraction of reads mapping to MT genes per cell") +
  theme(plot.title = element_text(size = 12)) +
  coord_fixed(ratio = 1)
p_umap_MT_fix <- set_panel_size(p_umap_MT, width = unit(12, "cm"), height = unit(12, "cm"))
# p_umap_MT
ggsave(path_qc %&% ".fractionMT.png", p_umap_MT_fix,
  width = 20, height = 15, units = "cm", dpi = 300
)

# plot UMAP (based on highly variable genes), colours = umap clusters
p_umap_uc <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = as.factor(umap_cl))) +
  geom_point(alpha = 0.7, size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 2, shape = 15), nrow = 16)) +
  scale_color_manual(name = "Umap clusters", values = cols33) +
  ggtitle(paste("Umap clustering. Modularity score: ", round(metadata(my_sce)$umap_modularity, digits = 4))) +
  theme(plot.title = element_text(size = 12)) +
  coord_fixed(ratio = 1)
p_umap_uc_fix <- set_panel_size(p_umap_uc, width = unit(12, "cm"), height = unit(12, "cm"))
# p_umap_uc
ggsave(path %&% ".umap_clustering.png", p_umap_uc_fix,
  width = 20, height = 15, units = "cm", dpi = 300
)


# have matrix with normcounts
# remove all genes that have 0 counts in all cells
normcounts_all.zero.removed <- normcounts(my_sce)
mask_all_zero <- apply(normcounts_all.zero.removed, 1, sum) > 0
normcounts_all.zero.removed <- normcounts_all.zero.removed[mask_all_zero, ]
stopifnot(length(normcounts_all.zero.removed[, 1]) == sum(apply(normcounts_all.zero.removed, 1, sum) > 0))

# plot ridge plot fraction MT vs celltypes
ridge1 <- ggplot(cell_attributes, aes(x = fractionMT, y = celltype_final, fill = celltype_final)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4) +
  scale_fill_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  guides(fill = guide_legend(nrow = 18))
ggsave(path_qc %&% ".ridge_fractionMT_per_final_celltype.png", ridge1,
  width = width_second, height = 15, units = "cm", dpi = 300
)

# plot ridge plot fraction MT vs clusters
ridge2 <- ggplot(cell_attributes, aes(x = fractionMT, y = phenograph_clusters, fill = phenograph_clusters)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4)
ggsave(path_qc %&% ".ridge_fractionMT_per_cluster.png", ridge2,
  dpi = 300
)

# plot ridge plot log_umi vs celltype
ridge3 <- ggplot(cell_attributes, aes(x = log_umi, y = celltype_final, fill = celltype_final)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4) +
  scale_fill_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  guides(fill = guide_legend(nrow = 18))
ggsave(path_qc %&% ".ridge_logUMI_per_final_celltype.png", ridge3,
  width = width_second, height = 15, units = "cm", dpi = 300
)

# plot ridge plot log_umi vs cluster
ridge4 <- ggplot(cell_attributes, aes(x = log_umi, y = phenograph_clusters, fill = phenograph_clusters)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4)
ggsave(path_qc %&% ".ridge_logUMI_per_cluster.png", ridge4,
  dpi = 300
)

# plot ridge plot n_gene vs celltype
ridge5 <- ggplot(cell_attributes, aes(x = n_gene, y = celltype_final, fill = celltype_final)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4) +
  scale_fill_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  guides(fill = guide_legend(nrow = 18))
ggsave(path_qc %&% ".ridge_n.gene_per_final_celltype.png", ridge5,
  width = width_second, height = 15, units = "cm", dpi = 300
)

# plot ridge plot n_gene vs cluster
ridge6 <- ggplot(cell_attributes, aes(x = n_gene, y = phenograph_clusters, fill = phenograph_clusters)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4)
ggsave(path_qc %&% ".ridge_n.gene_per_cluster.png", ridge6,
  dpi = 300
)

# plot ridge plot s_score vs celltype
ridge7 <- ggplot(cell_attributes, aes(x = s_score, y = celltype_final, fill = celltype_final)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4) +
  scale_fill_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  guides(fill = guide_legend(nrow = 18))
ggsave(path_qc %&% ".ridge_s.score_per_final_celltype.png", ridge7,
  width = width_second, height = 15, units = "cm", dpi = 300
)

# plot ridge plot s_score vs cluster
ridge8 <- ggplot(cell_attributes, aes(x = s_score, y = phenograph_clusters, fill = phenograph_clusters)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4)
ggsave(path_qc %&% ".ridge_s.score_per_cluster.png", ridge8,
  dpi = 300
)

# plot ridge plot g2m_score vs celltype
ridge9 <- ggplot(cell_attributes, aes(x = g2m_score, y = celltype_final, fill = celltype_final)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4) +
  scale_fill_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  guides(fill = guide_legend(nrow = 18))
ggsave(path_qc %&% ".ridge_g2m.score_per_final_celltype.png", ridge9,
  width = width_second, height = 15, units = "cm", dpi = 300
)

# plot ridge plot s_score vs cluster
ridge10 <- ggplot(cell_attributes, aes(x = g2m_score, y = phenograph_clusters, fill = phenograph_clusters)) +
  geom_density_ridges(alpha = 0.82, jittered_points = TRUE, point_alpha = 0.4, point_size = 0.4)
ggsave(path_qc %&% ".ridge_g2m.score_per_cluster.png", ridge10,
  dpi = 300
)

# plot expression of marker genes in UMAP
# read in gene list
gene.list.all <- read.table(opt$genelist, head = T, sep = "\t", stringsAsFactors = F)
# read protein/alias names of genes for alternate plot titles
if (opt$toggle_label) {
  # use regex to parse gene aliases in "notes" col of list file w/ format: ";ALIAS:alias_name;"
  gene.list.protein <- tapply(gene.list.all$notes, gene.list.all$SYMBOL, function(x) regmatches(x, regexec(";ALIAS:(.*?);", x))[[1]][2])
  # function for converting 'gene name' -> 'gene name/protein name'
  convert_gene_name <- function(gene_list) {
    gene_list <- unique(gene_list)
    out_list <- c()
    for (gene_symbol in gene_list) {
      gene_alias <- gene.list.protein[[gene_symbol]][1]
      if (!is.na(gene_alias) && gene_alias != gene_symbol) {
        final_name <- paste(gene_symbol, "/", gene_alias, sep = "")
      } else {
        final_name <- gene_symbol
      }
      out_list <- c(out_list, final_name)
    }
    return(out_list)
  }
}
# keep only genes that are in matrix
gene.list <- gene.list.all[gene.list.all$SYMBOL %in% rownames(normcounts_all.zero.removed), ]
# have list per group
gene.list <- tapply(gene.list$SYMBOL, gene.list$group, c)
gene.list.all <- tapply(gene.list.all$SYMBOL, gene.list.all$group, c)
# sort gene names alphabetically for better overview in plot
gene.list <- lapply(gene.list, sort)
gene.list.all <- lapply(gene.list.all, sort)
print("head(gene.list):")
print(head(gene.list))
# get indices of genes in matrix
list_groups <- lapply(seq_along(gene.list), function(x) {
  indices <- match(gene.list[[x]], rownames(normcounts_all.zero.removed))
  indices
})
names(list_groups) <- names(gene.list)
print("str(list_groups):")
print(str(list_groups))
# sort by position in matrix, do sanity check with old list
list_groups_sorted <- lapply(list_groups, sort)
list_groups_orig <- limma::ids2indices(gene.list, rownames(normcounts_all.zero.removed), remove.empty = F)
print("str(list_groups_orig):")
print(str(list_groups_orig))
# stopifnot(all.equal(list_groups_sorted, list_groups_orig))


# loop over groups of marker genes to plot gene expression on tSNEs with one png file per group

fontsize <- theme(axis.text = element_text(size = 3), axis.title = element_text(size = 3))
theme_set(theme_bw(4) + fontsize)
# plot UMAP (based on highly variable genes), colours = final cell type, as reference for expression plots
p_final_ref <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = celltype_final)) +
  scale_color_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  geom_point(size = 0.0001) +
  theme(aspect.ratio = 1) +
  # coord_fixed(ratio = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(legend.position = "none")

p_legend_ref <- ggplot(cell_attributes, aes(x = umap1, y = umap2, color = celltype_final)) +
  scale_color_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  geom_point(size = 1) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 0.8, shape = 15), nrow = 17)) +
  theme(
    legend.key.width = unit(0.05, "line"),
    legend.key.height = unit(0.001, "line"), # ,
    # legend.text = element_text(size = 3.5),
    legend.spacing.y = unit(0.001, "line"),
    legend.spacing.x = unit(0.25, "line")
  )

legend_ref <- cowplot::get_legend(p_legend_ref)
# ggsave(path_single_plots %&% ".legend_ct_colours_reference.png", legend_ref, dpi = 600, width = 3.5, height = 4, units = "cm")

for (group in seq(length(gene.list.all))) {
  group_name <- names(gene.list.all[group])
  print("group:")
  print(group_name)

  # get all genes that have zero counts in all cells
  mask.not.expressed <- !(gene.list.all[[group]] %in% rownames(normcounts_all.zero.removed[list_groups[[group_name]], , drop = F]))
  genes.not.expressed <- sort(gene.list.all[[group]][mask.not.expressed])
  if (opt$toggle_label) {
    genes.not.expressed <- convert_gene_name(genes.not.expressed)
  }
  # split number of genes per line based on # of characters of total line
  if (nchar(paste(genes.not.expressed[1:5], collapse = "")) > 32) {
    split_by <- 3
  } else {
    split_by <- 5
  }
  if (length(genes.not.expressed) > 0) {
    genes.not.expressed.split <- split(genes.not.expressed, ceiling(seq_along(genes.not.expressed) / split_by))
  } else {
    genes.not.expressed.split <- NULL
  }
  # genes.not.expressed.split <- split(genes.not.expressed, ceiling(seq_along(genes.not.expressed)/split_by))
  genes.not.expressed.split.collapsed <- lapply(seq_along(genes.not.expressed.split), function(x) {
    paste(genes.not.expressed.split[[x]], collapse = ", ")
  })

  # generate plot that consists of text and lists all genes that were not found in this sample
  if (length(genes.not.expressed.split.collapsed) > 0) {
    text1 <- paste(genes.not.expressed.split.collapsed, collapse = ",\n")
  } else {
    text1 <- "-"
  }
  text2 <- paste0("For the following genes \n no expression was detected in this sample:\n", text1)
  plot_missed <- ggplot() +
    annotate("text", x = 1, y = 1, label = text2, size = 1.2) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      aspect.ratio = 1
    )
  # plot_missed
  # ggsave(path_single_plots %&% "." %&% group_name %&% "_genes_not_expressed.png", plot_missed, dpi = 600, width = 3.5, height = 4, units = "cm")
  # check if there are genes expressed and a submatrix can be generated
  if (length(list_groups[[group_name]]) > 0) {
    # get matrix with all cells and all genes of the respective group
    matrix_group <- normcounts_all.zero.removed[list_groups[[group_name]], , drop = F]
  } else {
    matrix_group <- matrix(0, nrow = 1, ncol = ncol(normcounts(my_sce)))
  }
  # get empty list that will be filled with expression plots per gene
  plot_gene <- list()
  # check first if any of the group's genes are expressed in this sample
  if (length(list_groups[[group_name]]) > 0) {
    # loop over the expressed genes of the current group
    for (gene in seq(length(list_groups[[group_name]]))) {
      # empty list that will be filled with 1) expression plot and 2) violin plot of one gene
      merged_plots <- list()
      # combine cell_attributes with expression of current gene
      cell_attributes_gene <- cell_attributes
      cell_attributes_gene$normcounts <- matrix_group[gene, ]
      # print("table(cell_attributes_gene$normcounts < 0):")
      # print(table(cell_attributes_gene$normcounts < 0))
      # set all values below 0 to 0 for the plotting:
      cell_attributes_gene$normcounts[cell_attributes_gene$normcounts < 0] <- 0
      print("table(cell_attributes_gene$normcounts < 0):")
      print(table(cell_attributes_gene$normcounts < 0))
      # maximum count found for current gene
      max_count <- max(cell_attributes_gene$normcounts)
      # upper limit of gene expression colour scale is either maximum count or 3
      upper_limit <- max(3, max_count)
      # match gene symbol w/ conventional/protein name for plot labels
      gene_symbol <- rownames(matrix_group)[gene]
      if (opt$toggle_label) {
        legend <- convert_gene_name(gene_symbol)
      } else {
        legend <- gene_symbol
      }
      # plot gene expression
      merged_plots[["expr"]] <- ggplot(cell_attributes_gene, aes(x = umap1, y = umap2)) +
        geom_point(aes(color = normcounts), size = rel(0.001)) +
        xlab("UMAP 1") +
        ylab("UMAP 2") +
        #        scale_color_distiller(name="", palette = "Spectral", direction = -1,
        #        scale_color_viridis(option = "plasma", name = "",
        scale_color_gradientn(
          name = "", colours = c(
            "slateblue4", "royalblue1",
            "aquamarine3", "khaki", 383,
            "sienna1", "orangered4"
          ),
          limits = c(0, max(3, upper_limit)),
          breaks = c(floor(upper_limit / 3), round(2 * (upper_limit / 3)), upper_limit)
        ) +
        coord_fixed(ratio = 1) +
        ggtitle(legend) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
          plot.margin = margin(1, 1, 1, 1, "pt"),
          plot.background = element_rect(fill = "white", colour = "black", size = 0.2),
          # legend.text = element_text(size = 8),
          legend.key.size = unit(0.3, "line"),
          legend.text = element_text(size = rel(0.7)),
          legend.margin = margin(0.2, 0.2, 0.2, 0.2),
          axis.text = element_text(size = rel(0.6)),
          aspect.ratio = 1
        )
      merged_plots[["expr"]]
      # ggsave(path_single_plots %&% "." %&% group_name %&% "." %&% legend %&% "_expression_umap.png", merged_plots[["expr"]], dpi = 600, width = 3.5, height = 2.8, units = "cm")
      # scale_color_distiller(palette="RdBu")
      # plot violin plot
      merged_plots[["violin"]] <- ggplot(cell_attributes_gene, aes(x = factor(phenograph_clusters), y = normcounts)) +
        xlab("Phenograph cluster ID") +
        scale_fill_brewer(name = "", palette = "Set2") +
        scale_color_brewer(name = "", palette = "Set2") +
        geom_violin(scale = "count", alpha = 0.2, size = 0.2) +
        ggbeeswarm::geom_quasirandom(alpha = 0.2, groupOnX = T, size = 0.05) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 3),
          plot.margin = margin(3, 5, 6, 1),
          plot.background = element_rect(
            fill = "white",
            colour = "white",
            size = 1
          ),
          axis.text = element_text(size = 3),
          axis.title = element_text(size = 3.5)
        )
      # ggsave(path_single_plots %&% "." %&% group_name %&% "." %&% legend %&% "_expression_violin.png", merged_plots[["violin"]], dpi = 600, width = 3.5, height = 1.2, units = "cm")
      # combine expression plot and violin plot of current gene to one plot
      plot_gene[[gene]] <- cowplot::plot_grid(plotlist = merged_plots, ncol = 1, rel_heights = c(2.1, .9))
      # add combined plot to list that aggregates plots of all genes of the current group
      # plot_gene[[gene]]
    }
  }
  names(plot_gene) <- rownames(matrix_group)
  plot_gene <- c(list(plot_missed, p_final_ref, legend_ref), plot_gene)
  nr.cols <- min(6, 2 * (length(list_groups[[group_name]]) + 3))
  nr.rows <- ceiling((length(list_groups[[group_name]]) + 3) / nr.cols)
  plots_group <- cowplot::plot_grid(plotlist = plot_gene, ncol = nr.cols)
  # save png for each group
  ggsave(path_gene_expr %&% "." %&% group_name %&% "_gene_expression.png", plots_group, dpi = 600, width = 20, height = 4 * nr.rows, units = "cm")
  print(paste0("Gene expression plot of group ", group_name, " plotted."))
}


# set initial general options/parameters
options(stringsAsFactors = FALSE)
fontsize <- theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11))
theme_set(theme_bw(12) + fontsize)
col.pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)

# prepare data for plotting the bar plot
res <- colData(my_sce)
res$new.cell.type <- res$celltype_final
tab <- table(res$phenograph_clusters, res$new.cell.type)
top <- ncol(tab)
left <- nrow(tab)
tab <- addmargins(tab)
tab <- rbind(tab, f.round.preserve.sum(tab[nrow(tab), ] / tab[nrow(tab), ncol(tab)] * 100))

# cell type composition bar plot
# assign colours
t.plot <- tab[nrow(tab), seq(ncol(tab) - 1)]
t.plot <- data.frame(x = "Total " %&% tab[nrow(tab) - 1, ncol(tab)] %&% " cells", "Celltypes" = names(t.plot), value = t.plot)
t.plot$Celltypes <- as.factor(t.plot$Celltypes)

all.cts.bar <- levels(t.plot$Celltypes)[!levels(t.plot$Celltypes) %in% c("uncertain", "unknown")]
id.bar.ct <- match(levels(t.plot$Celltypes), names(ct.color))
id.bar.ct

# plot bar plot of cell types in this sample
p_bar <- ggplot(t.plot, aes(x = x, y = value, fill = Celltypes)) +
  xlab("") +
  ylab("Percent") +
  guides(fill = guide_legend(nrow = 22)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Cell type", values = ct.color[id.bar.ct], drop = F) +
  ggtitle(label = "Celltype composition")
# p_bar
ggsave(path %&% ".celltype_barplot.png", p_bar,
  width = 16, height = 18, units = "cm", dpi = 600
)

