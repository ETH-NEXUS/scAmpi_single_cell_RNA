##################################################
#### File name: sctransform_preprocessing.R
#### Author: Michael Prummer
#### Date created: June 2019
#### R Version: 4.0
##################################################


suppressPackageStartupMessages({
  library(optparse)
  library(rhdf5)
  library(scran)
  library(sctransform)
  library(plyr)
  library(ggplot2)
  library(uwot)
  library(igraph)
})
cat('using patched version of vst, with exception handling\n')
# Get command line arguments
args <- commandArgs(trailingOnly = FALSE)
# Find the argument containing the script path
script_arg <- grep("--file=", args, value = TRUE)
# Extract the script path from the argument
script_path <- sub("--file=", "", script_arg)
# Normalize the script path
script_path <- normalizePath(script_path)
# Get the directory of the script
script_dir <- dirname(script_path)
source(file.path(script_dir, 'vst_check.R'))

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")

option_list <- list(
  make_option("--inHDF5", type = "character", help = "Path to hdf5 input file. It includes raw expression matrix, gene & cell attributes."),
  make_option("--sample", type = "character", help = "Sample name."),
  make_option("--number_genes", type = "character", help = "Number of genes with the highest variance in the residuals that will be used for the calculation of umap coordinates and given out into an hdf5 file for the phenograph clustering."),
  make_option("--min_var", type = "character", help = "Minimum variance of the residuals for a gene to be used for the calculation of umap coordinates and given out into an hdf5 file for the phenograph clustering."),
  make_option("--n_nn", type = "character", help = "Number of nearest neighbours for the UMAP calculation."),
  make_option("--outdir", type = "character", help = "Path to output directory.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


dat <- h5read(opt$inHDF5, "raw_counts")

# For each gene, sum of counts
A <- rowSums(dat)
# For each cell, sum of counts
B <- colSums(dat)
# Check if there are all-zero columns or rows. If yes, do not proceed.
if (any(!(A > 0)) || any(!(B > 0))) {
  stop("The data seems to be unfiltered! Make sure the input of this script is a filtered count matrix.")
}

dat = dat[!is.na(A) & A , B]
# cell description table
cell_desc <- as.data.frame(h5read(opt$inHDF5, "cell_attrs"))
colnames(cell_desc)[colnames(cell_desc) == "cell_names"] <- "barcodes"
rownames(cell_desc) <- cell_desc$barcodes
str(cell_desc)
# gene description table
gene_desc <- as.data.frame(h5read(opt$inHDF5, "gene_attrs"))
gene_desc$SYMBOL <- gene_desc$gene_names
str(gene_desc)

# cap at max_count
max_count=as.numeric(strsplit(opt$max_count, 'x', fixed=T)[[1]])
do_capping=TRUE
if (length(max_count) > 1){
  cat(paste0('\nCapping count matrix at ',max_count[1],' quantile x ',max_count[2]))
  rpm=t(t(dat)/colSums(dat))*1e6
  rpm[rpm==0]=NA
  max_rpm=apply(rpm,MAR=1, FUN = quantile, na.rm=TRUE, names=FALSE, probs=max_count[1] )
  max_count=ceiling(max_rpm %*% t(colSums(dat))/1e6)*max_count[2]
  cat(paste0(' (',min(max_count),' to ',max(max_count),' UMIs)\n/n'))
} else if( max_count>0){
  cat(paste0('\nCapping count matrix at ',max_count,' UMIs\n\n'))
} else {
  cat('\n--max_count <= 0 - disable capping of count matrix\n\n')
  do_capping=FALSE
}
if (do_capping){
  n_cap=sum(rowSums(dat>max_count)>0)
  if (n_cap>0){
    cat(paste0('Capping counts for ',n_cap,' genes:\n'))
    cap_genes=gene_desc$SYMBOL[rowSums(dat>max_count)>0]
    cat(cap_genes,sep=', ')
    cat('\n')
    dat = pmin(dat, max_count)
  }
}

# exclude duplicated gene symbols
if (length(which(duplicated(gene_desc$SYMBOL))) > 0) warning("Some gene symbols are duplicated. Only the first is kept.")
keep <- setdiff(seq(nrow(gene_desc)), which(duplicated(gene_desc$SYMBOL)))
gene_desc <- gene_desc[keep, ]
dat <- dat[keep, ]
rownames(gene_desc) <- gene_desc$gene_ids
rownames(dat) <- rownames(gene_desc)
colnames(dat) <- rownames(cell_desc)
# calculate more gene_desc features
gene_desc$mean <- rowMeans(dat)
gene_desc$detection_rate <- rowMeans(dat > 0)
gene_desc$var <- apply(dat, 1, var)
# calculate more cell_desc features
cell_desc$n_umi <- colSums(dat)
cell_desc$n_gene <- colSums(dat > 0)
cell_desc$log_umi <- log(cell_desc$n_umi)
## cell cycle correction scores
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package = "scran"))
print("Start performing scran::cyclone: ")
cecy <- scran::cyclone(dat, pairs = hs.pairs)
cell_desc$g2m_score <- cecy$normalized.scores$G2M
cell_desc$s_score <- cecy$normalized.scores$S
cell_desc$cycle_phase <- cecy$phases
## variance stabilizing transformation:
set.seed(44)
print("Start performing patched  version of vst: ")
vst_out = vst(dat, cell_attr = cell_desc, method="nb_fast",
                           latent_var = c('log_umi'),
                           latent_var_nonreg = c("g2m_score", "s_score"),
                           return_gene_attr = T, return_cell_attr = T)
# potentially concerning warnings about iteration limit reached/ glm.fit algorithm did not converge
print("Start performing sctransform::smooth_via_pca: ")
y_smooth = sctransform::smooth_via_pca(vst_out$y, do_plot = FALSE, max_pc=opt$max_pc_smooth)
print("Start performing sctransform::correct: ")
dat_cor <- sctransform::correct(vst_out,
  data = y_smooth,
  do_round = TRUE, do_pos = FALSE
)

# Print string of the model:
print("Character representation of the model formula:")
print(vst_out$model_str)
print("Character representation of model for non-regularized variables:")
print(vst_out$model_str_nonreg)

# Print information about residual mean
print("str(vst_out$gene_attr$residual_mean:")
print(str(vst_out$gene_attr$residual_mean))
print("mean(vst_out$gene_attr$residual_mean):")
print(mean(vst_out$gene_attr$residual_mean))
print("median(vst_out$gene_attr$residual_mean):")
print(median(vst_out$gene_attr$residual_mean))

# Print information about residual variance
print("str(vst_out$gene_attr$residual_variance:")
print(str(vst_out$gene_attr$residual_variance))
print("mean(vst_out$gene_attr$residual_variance):")
print(mean(vst_out$gene_attr$residual_variance))
print("median(vst_out$gene_attr$residual_variance):")
print(median(vst_out$gene_attr$residual_variance))

# Plot mean variance plot
mvp <- ggplot(vst_out$gene_attr, aes(log10(gmean), log10(residual_variance))) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth()
ggsave(opt$outdir %&% opt$sample %&% ".mean_vs_variance_plot.png", mvp, dpi = 300)
mvp2 <- ggplot(vst_out$gene_attr, aes(log10(gmean), residual_variance)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_density_2d(size = 0.3) +
  geom_smooth()
ggsave(opt$outdir %&% opt$sample %&% ".mean_vs_variance_plot_2.png", mvp2, dpi = 300)


## make SingleCellExperiment object and save to disk
vst_out$gene_attr$gene_ids <- rownames(vst_out$gene_attr)
gene_desc <- merge(vst_out$gene_attr, gene_desc, sort = F, by = "gene_ids")
idx <- match(rownames(dat_cor), rownames(dat))
idy <- match(colnames(dat_cor), colnames(dat))
sce <- SingleCellExperiment(
  assays = list(
    counts = dat[idx, idy],
    normcounts = dat_cor,
    pearson_resid = vst_out$y
  ),
  colData = cell_desc,
  rowData = gene_desc
)


# find the genes with the highest variance in the residuals
# They will be used to calculate the UMAP coordinates
# Later: write them into an extra hdf5 file (for phenograph)
number_genes <- as.integer(opt$number_genes)
min_var <- as.integer(opt$min_var)
variable_genes <- gene_desc[order(gene_desc$residual_variance, decreasing = TRUE), ][seq(number_genes), ]

print("Number of genes that pass residual variance threshold:")
table(gene_desc$residual_variance > min_var)
print("Set value for maximum number of variable genes to consider:")
print(number_genes)

variable_genes <- variable_genes[variable_genes$residual_variance > min_var, ]
print("Number of variable genes that are considered for phenograph and UMAP:")
print(length(variable_genes$residual_variance))
# get corrected count matrix with only the highly variable genes
cor_counts_variable <- dat_cor[rownames(dat_cor) %in% variable_genes$gene_ids, ]
# make sure the gene names and ids are in the same order as in the matrix
cor_counts_variable_order <- as.data.frame(rownames(cor_counts_variable))
names(cor_counts_variable_order) <- "gene_ids"
variable_genes_ordered <- join(cor_counts_variable_order, variable_genes, by = "gene_ids")
stopifnot(rownames(cor_counts_variable) == variable_genes_ordered$gene_ids)
variable_genes_ordered$gene_ids <- as.character(variable_genes_ordered$gene_ids)

# calculate UMAP coordinates on pearson residuals of the variable genes and store in sce object
n_nn <- as.integer(opt$n_nn)
residuals_matrix <- assay(sce, "pearson_resid")
print("head(rownames(residuals_matrix)):")
print(head(rownames(residuals_matrix)))
residuals_matrix_variable <- residuals_matrix[rownames(residuals_matrix) %in% variable_genes$gene_ids, ]
set.seed(184)
umap.hvg <- umap(t(residuals_matrix_variable), n_neighbors = n_nn, pca = 50, spread = 1, min_dist = 0.3, ret_nn = T)
reducedDim(sce, "umap_hvg") <- umap.hvg$embedding
metadata(sce) <- c(metadata(sce), list(umap_hvg = umap.hvg$nn$euclidean))
## nearest neighbors
uu.nn <- umap.hvg$nn$euclidean$idx
uu.nn[uu.nn == 0] <- 1
## weights as 1-distance
wgt <- 1 - umap.hvg$nn$euclidean$dist / min(max(umap.hvg$nn$euclidean$dist), 1e5)
wgt[wgt < -1] <- 0
## convert to adjacancy matrix
adj <- matrix(0, ncol(residuals_matrix_variable), ncol(residuals_matrix_variable))
rownames(adj) <- colnames(residuals_matrix_variable)
colnames(adj) <- colnames(residuals_matrix_variable)
for (i in seq_len(ncol(residuals_matrix_variable))) {
  adj[i, colnames(residuals_matrix_variable)[uu.nn[i, ]]] <- wgt[i, ]
}
# convert to weighted graph
set.seed(0)
g <- graph.adjacency(adj, mode = "undirected", weighted = T, diag = F)
km <- igraph::cluster_louvain(g)
ph.membership <- km$membership
names(ph.membership) <- km$names
colData(sce)$umap_cl <- km$membership
metadata(sce) <- c(metadata(sce), list(umap_modularity = modularity(km)))


# save sce object in RDS file
saveRDS(sce, opt$outdir %&% opt$sample %&% ".corrected.RDS")
print("str(sce):")
print(str(sce))

# write output hdf5 file
outfile <- opt$outdir %&% opt$sample %&% ".corrected.variable_genes.h5"
h5createFile(outfile)
h5createGroup(outfile, "cell_attrs")
h5createGroup(outfile, "gene_attrs")
h5write(variable_genes_ordered$gene_ids, outfile, "gene_attrs/gene_ids")
h5write(variable_genes_ordered$gene_names, outfile, "gene_attrs/gene_names")
h5write(cell_desc$barcodes, outfile, "cell_attrs/cell_names")

# set chunk size for writing h5 to c(1000,1000) or if the matrix is smaller to its dimensions (default)
if (dim(cor_counts_variable)[1] > 1000 && dim(cor_counts_variable)[2] > 1000) {
  chunks <- c(1000, 1000)
} else {
  chunks <- dim(cor_counts_variable)
}
cat("\n\nhdf5 chunk size:", chunks, "\n\n")
h5createDataset(file = outfile, dataset = "cor_counts", dims = dim(cor_counts_variable), chunk = chunks)
h5write(cor_counts_variable, outfile, "cor_counts")

# Plot the ranking of the variance of the residuals
plot_variance_genes <- gene_desc[order(gene_desc$residual_variance, decreasing = TRUE), ]
plot_variance_genes$rank <- seq(length(plot_variance_genes$gene_ids))
print("head(plot_variance_genes):")
print(head(plot_variance_genes))

plot_1 <- ggplot(data = plot_variance_genes, aes(x = rank, y = residual_variance)) +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 20000, 1000)) +
  geom_vline(xintercept = number_genes) +
  geom_hline(yintercept = min_var)
ggsave(
  filename = paste0(opt$outdir, opt$sample, ".plot_ranking_variance_of_residuals_per_gene.png"),
  plot = plot_1,
  width = 30, height = 20, dpi = 600, units = "cm"
)
