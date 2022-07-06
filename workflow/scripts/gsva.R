################################################################################
## GSVA analysis
## Author: Michael Prummer
## Adapted: Anne Bertolini
################################################################################

lby <- c("optparse", "scran", "reshape2", "uwot", "GSVA", "aroma.light",
        "ggplot2", "pheatmap", "RColorBrewer", "limma", "viridis", "dplyr")
resp <- lapply(lby, require, character.only = T, warn.conflicts = F, quietly = T)
if (!all(unlist(resp))) stop("Could not load one or more packages")
rm(resp, lby)

# general options
options(stringsAsFactors = FALSE)
fontsize <- theme(axis.text = element_text(size = 9), axis.title = element_text(size = 11))
theme_set(theme_bw(12) + fontsize)
col.pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")

# command line arguments are parsed
option_list <- list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_celltypes_noatypical.RDS)."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--geneset", type = "character", help = "Hallmark geneset library gmt file."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

path <- opt$outputDirec %&% opt$sampleName


################################################################################
## main code starts here
################################################################################

## load input data
print(opt$SCE)
sce_data <- readRDS(opt$SCE)

# define color scheme
all.cell.types <- metadata(sce_data)$all_celltypes
cols33 <- c("red2", "green4", "blue2", "cyan2", "yellow1", "purple", "brown",
            "chocolate1", "chartreuse2", "darkgoldenrod3", "steelblue1", "slateblue3", "olivedrab4", "gold2",
            "violetred3", "darkcyan", "orchid3", "darksalmon", "darkslategrey", "khaki", "indianred2", "magenta", "slategray2",
            "olivedrab1", "mediumaquamarine", "hotpink", "yellow3",
            "bisque4", "darkseagreen1", "dodgerblue3",
            "deeppink4", "sienna4", "mediumorchid4")
number_celltypes <- length(all.cell.types)
ct.color <- c(cols33[seq(number_celltypes)], "grey50", "black")

################################################################################

# load ref.gene.list
# make sure all columns are captured when reading in the gene set lists in gmt format
tmp <- readLines(opt$geneset)
tmp <- lapply(tmp, function(x) strsplit(x, "\\\t")[[1]])
names(tmp) <- sapply(tmp, function(x) x[1])
# remove gene set names from gene list
gset <- sapply(tmp, function(x) x[-1])
# remove gene set description from gene list
gset <- sapply(gset, function(x) x[-1])
names(gset) <- gsub("HALLMARK_", "", names(gset))

# select subset gene expression matrix
idxs <- ids2indices(gset, rowData(sce_data)$SYMBOL)
all.genes <- unique(as.numeric(unlist(idxs)))
t.m <- assay(sce_data, "pearson_resid")[all.genes, ]

# estimate geneset-sample matrix from gene-sample matrix
rgsa <- gsva(t.m, gset, method = "gsva")
# prepare gsva results for plotting
gsva_results <- melt(rgsa)
names(gsva_results) <- c("gene.set", "barcodes", "value")

# get UMAP coordinates
umap_coord <- reducedDim(sce_data, "umap_hvg")
umap_coord <- as.data.frame(umap_coord)
umap_coord$barcodes <- rownames(umap_coord)

# add UMAP coordinates to cell meta data
tomerge <- dplyr::left_join(umap_coord, as.data.frame(colData(sce_data)), by = "barcodes")
stopifnot(tomerge$barcodes == as.data.frame(colData(sce_data))$barcodes)

# merge gsva results with cell meta data
# preserve order of cells in `gsva_results`
t.plot <- dplyr::left_join(gsva_results, tomerge, by = "barcodes")
stopifnot(t.plot$barcodes == gsva_results$barcodes)

# Trim outlier values
t.plot$value_limited <- t.plot$value
perc_1 <- quantile(t.plot$value, prob = 0.01)
perc_99 <- quantile(t.plot$value, prob = 0.99)
t.plot$value_limited[t.plot$value_limited <= perc_1] <- perc_1
cat("\n\nNumber of data points set to perc_1:\n")
print(table(t.plot$value_limited == perc_1))

t.plot$value_limited[t.plot$value_limited >= perc_99] <- perc_99
cat("\n\nNumber of data points set to perc_99:\n")
print(table(t.plot$value_limited == perc_99))

cat("\nsummary(t.plot$value_limited) after trimmming the outliers:")
print(summary(t.plot$value_limited))

t.plot$value_limited <- round(t.plot$value_limited, digits = 2)
break_low <- min(t.plot$value_limited)
cat("\nbreak_low:", break_low)
break_high <- max(t.plot$value_limited)
cat("\nbreak_high:", break_high)
label_low <- paste("<=", break_low)
cat("\nlabel_low:", label_low)
label_high <- paste(">=", break_high)
cat("\nlabel_high:", label_high)

# plot UMAP with gsva results
pp <- ggplot(t.plot, aes(x = V1, y = V2)) +
  geom_point(aes(color = value_limited), size = 1) +
  xlab("umap-1") +
  ylab("umap-2") +
  theme(panel.background = element_rect(fill = "grey80")) +
  scale_color_distiller(name = "", palette = "RdBu",
                        breaks = c(break_low, 0, break_high),
                        labels = c(label_low, "0.00", label_high)) +
  facet_wrap(~gene.set, ncol = 6)
ggsave(path %&% ".gsetscore_hvg_umap.png", pp,
       dpi = 600, width = 50, height = 50, units = "cm")


# save sce object for gsva-scores instead of gene counts for later use
sce_gsva <- SingleCellExperiment(assays = list(gsva = rgsa),
                                colData = colData(sce_data),
                                metadata = gset)
saveRDS(sce_gsva, path %&% ".sce_gsva.RDS")


################################################################################
## heatmaps of geneset-cell matrix, with phcl-annotation.

# prepare annotation column for heatmap
annot.col <- data.frame(Phenograph = factor(tomerge$phenograph_clusters),
                       Cell.type = tomerge$celltype_final)
rownames(annot.col) <- tomerge$barcodes
stopifnot(rownames(annot.col) == colnames(rgsa))

# have clusters sorted to colours
all_clusters <- levels(annot.col$Phenograph)
number_clusters <- length(all_clusters)
cluster_colors <- cols33[seq(number_clusters)]
names(cluster_colors) <- all_clusters

# have cell types sorted to colours
celltype_colors <- ct.color
names(celltype_colors) <- levels(annot.col$Cell.type)

# have colours for clusters and cell types together
plotting_colors <- list(Phenograph = cluster_colors, Cell.type = celltype_colors)

col.pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)

# have gsva results ordered by clusters
rgsa_plotting <- rgsa[, order(annot.col$Phenograph)]

# plot heatmap
hm1 <- pheatmap(rgsa_plotting,
               scale = "none",
               clustering_method = "ward.D2",
               show_colnames = F,
               color = col.pal,
               cluster_cols = F,
               annotation_col = annot.col,
               annotation_colors = plotting_colors,
               annotation_names_col = T,
               fontsize_row = 8,
               silent = TRUE)

ggsave(path %&% ".gsetscore_hm.png", hm1$gtable,
       width = 30, height = 24, units = "cm", dpi = 600)
