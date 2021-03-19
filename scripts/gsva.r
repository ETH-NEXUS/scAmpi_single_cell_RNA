
################################################################################
## cell type classification - GSVA
## NOTE: this is an independent piece of analysis. 
##       all output files, including plots, will already be handled in this script
################################################################################

lby = c("optparse", "scran", "reshape2", "uwot", "GSVA", "aroma.light",
        "ggplot2", "pheatmap","RColorBrewer", "limma", "viridis")
resp = lapply(lby, require, character.only=T, warn.conflicts=F, quietly=T)
if(!all(unlist(resp))) stop("Could not load one or more packages")
rm(resp, lby)

# options(error=function()traceback(2))
options(stringsAsFactors = FALSE)
fontsize = theme(axis.text=element_text(size=9), axis.title=element_text(size=11))
theme_set(theme_bw(12) + fontsize)
col.pal = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# convenience function for string concatenation
'%&%' = function(a,b) paste(a,b,sep="")

# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_celltypes_noatypical.RDS)."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--geneset", type = "character", help = "Hallmark geneset library gmt file."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

path = opt$outputDirec %&% opt$sampleName

print(opt$SCE)

path = opt$outputDirec %&% opt$sampleName


################################################################################
## main code starts here
################################################################################
## load input data
sce_data = readRDS(opt$SCE)

# define color scheme
all.cell.types = metadata(sce_data)$all_celltypes
cols33 <- c("red2", "green4", "blue2", "cyan2", "yellow1", "purple", "brown",
            "chocolate1", "chartreuse2", "darkgoldenrod3", "steelblue1", "slateblue3", "olivedrab4", "gold2",
            "violetred3", "darkcyan","orchid3", "darksalmon", "darkslategrey", "khaki", "indianred2","magenta", "slategray2", 
            "olivedrab1", "mediumaquamarine", "hotpink", "yellow3",
            "bisque4", "darkseagreen1", "dodgerblue3",
            "deeppink4", "sienna4", "mediumorchid4")
ct.color = c(cols33[seq(length(all.cell.types))], "grey50", "black")

################################################################################
#
# signature expression painting on tSNE plots (phenograph)
#
# load ref.gene.list
# make sure all columns are captured when reading in the gene set lists in gmt format
tmp = readLines(opt$geneset)
tmp = lapply(tmp, function(x) strsplit(x, "\\\t")[[1]])
names(tmp) = sapply(tmp, function(x) x[1])
# remove gene set names from gene list
gset = sapply(tmp, function(x) x[-1])
# remove gene set description from gene list
gset = sapply(gset, function(x) x[-1])
names(gset) = gsub("HALLMARK_", "", names(gset))

# select subset gene expression matrix
idxs = ids2indices(gset, rowData(sce_data)$SYMBOL)
all.genes = unique(as.numeric(unlist(idxs)))
             
t.m = assay(sce_data, "pearson_resid")[all.genes,]
             
# estimate geneset-sample matrix from gene-sample matrix
rgsa = gsva(t.m, gset, method="gsva")
# plot
t.plot = melt(rgsa)
names(t.plot) = c("gene.set", "barcodes", "value")
tomerge = reducedDim(sce_data, "umap_hvg")
tomerge = as.data.frame(tomerge)
tomerge$barcodes = rownames(tomerge)
tomerge = merge(tomerge, as.data.frame(colData(sce_data)))
t.plot = merge(t.plot, tomerge)
print("str(t.plot):")
print(str(t.plot))
print("summary(t.plot):")
print(summary(t.plot))

# Trim outlier values
t.plot$value_limited <- t.plot$value
perc_1 <- quantile(t.plot$value, prob = 0.01)
perc_99 <- quantile(t.plot$value, prob = 0.99)
t.plot$value_limited[t.plot$value_limited <= perc_1] <- perc_1
print("Number of data points set to perc_1:")
print(table(t.plot$value_limited == perc_1))
t.plot$value_limited[t.plot$value_limited >= perc_99] <- perc_99
print("Number of data points set to perc_99:")
print(table(t.plot$value_limited == perc_99))
print("str(t.plot$value_limited):")
print(str(t.plot$value_limited))
print("summary(t.plot$value_limited) after trimmming the outliers:")
print(summary(t.plot$value_limited))
t.plot$value_limited <- round(t.plot$value_limited, digits = 2)
break_low <- min(t.plot$value_limited)
print(break_low)
break_high <- max(t.plot$value_limited)
print(break_high)
label_low <- paste("<=", break_low)
print(label_low)
label_high <- paste(">=", break_high)
print(label_high)
#t.plot$value_log <- log10(t.plot$value + 1)
print("summary(t.plot$value_limited) after rounding the values to two digits:")
print(summary(t.plot$value_limited))

pp = ggplot(t.plot, aes(x=V1, y=V2) ) +
  geom_point(aes(color=value_limited), size=1) + xlab("umap-1") + ylab("umap-2") +
  theme(panel.background = element_rect(fill = "grey80")) +
  scale_color_distiller(name="", palette = "RdBu",
                        breaks = c(break_low, 0, break_high),
                        labels = c(label_low, "0.00", label_high)) +
  facet_wrap(~gene.set, ncol=6)
ggsave(path %&% ".gsetscore_hvg_umap.png", pp,
       dpi = 600, width = 50, height = 50, units = "cm")

# pp = ggplot(t.plot, aes(x=V1, y=V2, color=celltype_final) ) + 
#   geom_point(size=1) + xlab("umap-1") + ylab("umap-2") + 
#   guides(colour = guide_legend(override.aes = list(size=2, shape=15), ncol=1))+
#   scale_color_manual(name="", values = ct.color, drop=F)
# ggsave(path %&% ".ctfinal_hvg_umap.png", pp,
#        dpi = 600, width = 20, height = 15, units = "cm")
# 
# pp = ggplot(t.plot, aes(x=V1, y=V2, color=factor(umap_cl)) ) + 
#   geom_point(size=1) + xlab("umap-1") + ylab("umap-2") + 
#   guides(colour = guide_legend(override.aes = list(size=2, shape=15), ncol=1))+
#   scale_color_manual(name="", values = ct.color, drop=F)
# ggsave(path %&% ".hvgcl_hvg_umap.png", pp,
#        dpi = 600, width = 20, height = 15, units = "cm")
# 
# pp = ggplot(t.plot, aes(x=V1, y=V2, color=factor(phenograph_clusters)) ) + 
#   geom_point(size=1) + xlab("umap-1") + ylab("umap-2") + 
#   guides(colour = guide_legend(override.aes = list(size=2, shape=15), ncol=1))+
#   scale_color_manual(name="", values = ct.color, drop=F)
# ggsave(path %&% ".phcl_hvg_umap.png", pp,
#        dpi = 600, width = 20, height = 15, units = "cm")

# save sce object for gsva-scores instead of gene counts for later use 
sce_gsva = SingleCellExperiment(assays=list(gsva = rgsa), 
                                colData = colData(sce_data),
                                metadata = gset)
saveRDS(sce_gsva, path %&% ".sce_gsva.RDS")
# optional: save geneset-cell score matrix.
# write.table(rgsa, path %&% ".geneset_scores.tsv", sep="\t", col.names = NA)


################################################################################
## heatmaps of geneset-cell matrix, with phcl-annotation.
annot.col = data.frame(Phenograph=factor(tomerge$phenograph_clusters),
                       Cell.type=tomerge$celltype_final)
# levels(annot.col$Phenograph) = levels(annot.col$Phenograph) %&% " (" %&% 
#   clu.pu[-nrow(clu.pu)+c(0,1), ncol(clu.pu)-1] %&% ")"
rownames(annot.col) = colnames(rgsa)
ancol.col = cols33[seq(length(levels(annot.col$Phenograph)))]
names(ancol.col) = levels(annot.col$Phenograph)
ctan.col = ct.color
names(ctan.col) = levels(annot.col$Cell.type)
ancolor = list(Phenograph=ancol.col, Cell.type=ctan.col)
col.pal = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255)
rgsa = rgsa[ , order(annot.col$Phenograph)]
hm1 = pheatmap(rgsa, scale="none", clustering_method = "ward.D2", 
               show_colnames = F, color = col.pal, cluster_cols = F,
               annotation_col = annot.col, annotation_colors = ancolor,
               annotation_names_col = T, fontsize_row = 8)
ggsave(path %&% ".gsetscore_hm.png", hm1$gtable,
       width = 30, height = 24, units = "cm", dpi = 600)

