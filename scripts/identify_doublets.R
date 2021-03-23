##################################################
#### File name: remove_doublets.R
#### Date created: December 2020
#### R Version: 4.0.2
##################################################

library(scDblFinder)
library(optparse)
library(rhdf5)
library(SingleCellExperiment)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# parse command line arguments
option_list = list(
make_option("--hdf5File", type = "character", help = "Path to hdf5 input file containing raw counts after the Cellranger run"),
make_option("--sample", type = "character", help = "Sample name"),
make_option("--outdir", type = "character", help = "Path to output directory where resulting files are created")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# read in raw umi counts
cat("\n\n\nContent of input h5 file:\n\n")
print(h5ls(opt$hdf5File))
umi_counts <- h5read(opt$hdf5File, "raw_counts")
cat("\n\n\nstr(umi_counts)\n\n")
print(str(umi_counts))

# generate SCE object
cat("\n\n\n### Generate SCE object\n\n")
my_sce <- SingleCellExperiment(list(counts = umi_counts),
                               colData = DataFrame(barcodes = h5read(opt$hdf5File, "cell_attrs/cell_names")),
                               rowData = DataFrame(gene_ids = h5read(opt$hdf5File, "gene_attrs/gene_ids")))
rownames(my_sce) <- h5read(opt$hdf5File, "gene_attrs/gene_ids")
colnames(my_sce) <- h5read(opt$hdf5File, "cell_attrs/cell_names")
my_sce

# identify doublets
cat("\n\n\n### Identify doublets\n\n")
set.seed(123)
my_sce <- scDblFinder::scDblFinder(my_sce, verbose = TRUE)
doublet_cells <- colnames(my_sce)[my_sce$scDblFinder.class == "doublet"]
print(table(my_sce$scDblFinder.class))
rm(my_sce)

# save doublet barcodes
doublet_cells_file <- paste0(opt$outdir, opt$sample, ".doublet_barcodes.txt")
cat("\n\n\nOutput file with doublet barcodes:\n", doublet_cells_file, "\n")
write.table(doublet_cells, file = doublet_cells_file,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
