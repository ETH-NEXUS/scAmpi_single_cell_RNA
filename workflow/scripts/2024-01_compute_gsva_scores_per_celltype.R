################################################################################
## Run GSVA
################################################################################
# Run GSVA analysis on an SCE object in RDS file.
# Script does not generate any plots.
# Anne Bertolini
# Jan 2024

library(optparse)
library(reshape2)
library(GSVA)
library(limma)
library(SingleCellExperiment)
library(stringr)

options(stringsAsFactors = FALSE)

# command line arguments are parsed
option_list <- list(
  make_option("--SCE", type = "character", help = "Path to SCE object file with input data."),
  make_option("--subset_celltype", type = "character", help = "Celltype(s) that will be included in the analysis. Several cell types can be specified, comma separated."),
  make_option("--outdir", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--geneset", type = "character", help = "Geneset library gmt file."),
  make_option("--sample", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--min_set_size", type = "integer", help = "Minimum gene set size for it to be considered."),
  make_option("--method", type = "character", help = "Method to calculate GSVA score. See GSVA package documentation.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# give out date, time, session info, and input files
cat("\n\n")
print(Sys.time())
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n")
cat("\nInput files:\n\n")
print(opt)
cat("\n\n")


## load input data
my_sce <- readRDS(opt$SCE)
print(my_sce)

## identify tumour cells only
keep_celltypes <- unlist(stringr::str_split(opt$subset_celltype, pattern = ","))
cat("\n\n\n### Only processing cells of type:\n\n")
print(keep_celltypes)
cat("\n\n\n")

# subset SCE object to contain only tumour cells
mask_celltype <- colData(my_sce)$celltype_major %in% keep_celltypes
cat("### Cells included in analysis:\n")
print(table(mask_celltype))
my_sce <- my_sce[, mask_celltype]

# report number of cells
if (as.integer(dim(my_sce)[2]) == 0) {
  cat("\n### Sample", opt$sample, ": Zero cells\n\n")
} else {
  cat("\n### Sample", opt$sample, ":",  dim(my_sce)[2], "cells\n\n")
}

################################################################################

# read gmt file with gene sets
# make sure all columns are captured when reading in the gene set lists in gmt format
# readLines results in a character vector with one string per gene set

tmp <- readLines(opt$geneset)
# turn tmp into list of character vectors
tmp <- lapply(tmp, function(x) strsplit(x, "\\\t")[[1]])
names(tmp) <- lapply(tmp, function(x) x[1])
# remove gene set names from gene list
gset <- lapply(tmp, function(x) x[-1])
# remove gene set description from gene list
gset <- lapply(gset, function(x) x[-1])
cat("\n### Gene sets analysed:\n\n")
names(gset)
cat("\n\n\n")

# check how many genes of the gene sets are found in the matrix
gset_symbols <- unique(unlist(gset))

# check overlap and differences
setdiff_sce_gset <- setdiff(rownames(my_sce), gset_symbols)
setdiff_gset_sce <- setdiff(gset_symbols, rownames(my_sce))
stopifnot(length(intersect(setdiff_sce_gset, setdiff_gset_sce)) == 0)


cat("\n\n### Number of genes in count matrix but not in gene sets:", length(setdiff(rownames(my_sce), gset_symbols)), "\n\n")
cat("\n\n### Number of genes in gene sets but not in count matrix:", length(setdiff(gset_symbols, rownames(my_sce))), "\n\n")
print(setdiff(gset_symbols, rownames(my_sce)))

# have list of gene sets, HGNC symbol strings replaced with index of respective gene in SCE object
idxs <- limma::ids2indices(gset, rownames(my_sce))
# get unique numeric vector of all indices, occurring in gset and my_sce
all_genes <- unique(as.numeric(unlist(idxs)))
cat("\n### Number of unique genes in all gene sets:", length(gset_symbols), "\n\n")
cat("\n### Number of genes found in matrix:", length(all_genes), "\n\n")

# it was decided to not reduce the matrix to the sum of genes of all gene sets and keep all genes
#t.m <- assay(my_sce, "normcounts")[all_genes, ]

count_matrix <- assay(my_sce, "normcounts")
cat("\n### Number of genes in gene sets with non-zero expression in matrix:", table(rowSums(count_matrix) == 0)[1], "\n\n")
cat("\n### Number of genes in gene sets with zero expression in all cells:", table(rowSums(count_matrix) == 0)[2], "\n\n")

# remove genes that have zero counts in all cells
mask_genes_all_zero <- rowSums(count_matrix) != 0
count_matrix <- count_matrix[mask_genes_all_zero, ]

# estimate geneset-sample matrix from gene-sample matrix
# GSVA
gsva_result <- gsva(count_matrix,
                    gset,
                    method = opt$method,
                    min.sz = opt$min_set_size)

# reformat results
gsva_long <- melt(gsva_result)
names(gsva_long) <- c("gene_set", "barcode", "value")
print("str(gsva_long):")
print(head(gsva_long))

cat("\n### Number of gene sets with results available for sample ", opt$sample, ":", length(unique(gsva_long$gene_set)), "\n\n")

# write results into file
filename <- paste0(opt$outdir,'/', opt$sample, "_celltype_GSVA.tsv")
write.table(gsva_long, file = filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
