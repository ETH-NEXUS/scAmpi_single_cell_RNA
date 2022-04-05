# Compare cell barcodes and gene IDs of two hdf5 files containing scRNA count matrices

library(optparse)
library(rhdf5)

option_list <- list(
  make_option("--file_1", type = "character", help = "Path to hdf5 input file."),
  make_option("--file_2", type = "character", help = "Path to hdf5 input file.")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


h5ls(opt$file_1)
h5ls(opt$file_2)

# compare cells
barcodes_1 <- h5read(opt$file_1, "cell_attrs/cell_names")
barcodes_2 <- h5read(opt$file_2, "cell_attrs/cell_names")

cat("\n\n\nCheck if cells in file 1 and file 2 are identical:\n\n")
all.equal(barcodes_1, barcodes_2)
cat("\n\n\nCells in file 1 but not file 2:\n\n")
setdiff(barcodes_1, barcodes_2)
cat("\n\n\nCells in file 2 but not file 1:\n\n")
setdiff(barcodes_2, barcodes_1)

# compare genes
genes_1 <- h5read(opt$file_1, "gene_attrs/gene_ids")
genes_2 <- h5read(opt$file_2, "gene_attrs/gene_ids")

cat("\n\n\nCheck if genes in file 1 and file 2 are identical:\n\n")
all.equal(genes_1, genes_2)
cat("\n\n\nGenes in file 1 but not file 2:\n\n")
setdiff(genes_1, genes_2)
cat("\n\n\nGenes in file 2 but not file 1:\n\n")
setdiff(genes_2, genes_1)
