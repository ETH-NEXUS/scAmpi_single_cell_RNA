####    File name: filter_genes_and_cells.R
####    Author: Anne Richter
####    October 2018, adapted March 2022, adapted June 2024
####    R Version: 4.0

# This script takes an hdf5 file of raw expression counts (right after cellranger)

# Filtering Steps:
# 1. Cells are filtered if number of detected genes (NODG) is too low.
# 2. Cells are filtered if the fraction of reads mapped to MT- genes is too high.
# 3. Genes are filtered so that only protein-coding genes are kept (optional).
#    Note: the filtering for protein-coding genes is currently only enabled for human data.
# 4. Genes are filtered so that MT- genes are removed
# 5. Genes are filtered so that genes encoding for ribosomal proteins are removed

###   R libraries   ####

library(rhdf5)
library(optparse)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
suppressMessages(library(scater))

# convenience function for NOT %in%
"%!in%" <- function(x, y) !("%in%"(x, y))


###   Parse command line arguments   ####

option_list <- list(
  make_option("--hdf5File", type = "character", help = "Path to hdf5 input file. It includes raw expression matrix, gene & cell attributes."),
  make_option("--sample", type = "character", help = "Sample name, prefix of all output files"),
  make_option("--doublet_barcodes", type = "character", help = "Path to text file with doublet barcodes"),
  make_option("--remove_doublets", type = "logical", help = "Set TRUE or FALSE if doublets should be removed"),
  make_option("--nmads_NODG", type = "character", help = "Number of median-absolute-deviations away from median required for a value to be called an outlier, e.g. 5"),
  make_option("--nmads_fractionMT", type = "character", help = "Number of median-absolute-deviations away from median required for a value to be called an outlier, e.g. 5"),
  make_option("--outDir", type = "character", help = "Full path to output directory"),
  make_option("--genomeVersion", type = "character", help = "Specify the genome annotation version, either hg19 or GRCh38 are supported. The default is GRCh38."),
  make_option("--threshold_NODG", type = "character", help = "Hard threshold that gives the minimum NODG a cell must have to be further processed. E.g. 250"),
  make_option("--threshold_fractionMT", type = "character", help = "Hard threshold that gives the maximum fraction of MT reads a cell can have and be further processed. E.g. 0.5"),
  make_option("--minNumberCells", type = "character", help = "Minimum number of cells that should express a gene for it to be kept. E.g. 20"),
  make_option("--protein_coding_only", type = "logical", help = "TRUE or FALSE. Keep in matrix protein-coding genes only.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


###   Give out (session) info   ####

cat("\n\n")
print(Sys.time())
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n")
cat("\n###   Input:\n\n")
print(opt)
cat("\n\n")


###   Read in data   ####

# get file name of input hdf5 and sample name
file_name <- basename(opt$hdf5File)
sample <- opt$sample
outdir <- opt$outDir

# print information about input hdf5 file
cat("\n\n###   Content of input h5 file:\n\n")
print(h5ls(opt$hdf5File))

# command line parameters
nmads_fractionMT <- as.integer(opt$nmads_fractionMT)
threshold_MT <- as.numeric(opt$threshold_fractionMT)
nmads_NODG <- as.integer(opt$nmads_NODG)
threshold_NODG <- as.numeric(opt$threshold_NODG)
minNumberCells <- as.integer(opt$minNumberCells)

# get count matrix
umi_matrix <- h5read(opt$hdf5File, "raw_counts")
# link cell and gene information to matrix
rownames(umi_matrix) <- h5read(opt$hdf5File, "gene_attrs/gene_ids")
colnames(umi_matrix) <- h5read(opt$hdf5File, "cell_attrs/cell_names")

# have data frame with cell infos (barcodes)
cell_info <- data.frame(barcodes = h5read(opt$hdf5File, "cell_attrs/cell_names"))

# have data frame with gene info
gene_info <- data.frame(ensembl_gene_id = h5read(opt$hdf5File, "gene_attrs/gene_ids"))
gene_info$hgnc_symbol <- h5read(opt$hdf5File, "gene_attrs/gene_names")
# have ensembl IDs as rownames
rownames(gene_info) <- gene_info$ensembl_gene_id


###   FILTER CELLS   ####

cat("\n\n\n\n#####          FILTER CELLS          #####\n")


###   Filter doublets   ####

# check if doublets should be removed or not (command line parameter)
stopifnot(is.logical(opt$remove_doublets))
if (opt$remove_doublets) {
  cat("\n\n###   Filter doublets\n")
  doublet_barcodes <- read.csv(opt$doublet_barcodes,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  doublet_barcodes <- doublet_barcodes$V1
  stopifnot(colnames(umi_matrix) == cell_info$barcodes)
  cell_info$is_doublet <- cell_info$barcodes %in% doublet_barcodes
} else if (!opt$remove_doublets) {
  cat("\n\n###   No doublet filtering performed\n")
  doublet_barcodes <- character()
  cell_info$is_doublet <- cell_info$barcodes %in% doublet_barcodes
  stopifnot(sum(cell_info$is_doublet) == 0)
}


###   Filter cells with fraction of MT-reads   ####

cat("\n\n###   Filter cells with fraction of reads mapped to MT-genes\n")

# get all mitochondrial genes
gene_info$is_mt <- grepl("^MT-", gene_info$hgnc_symbol)
mt <- gene_info[gene_info$is_mt, "ensembl_gene_id"]
cat("\n### Number of MT- genes: ", length(mt), "\n")
cat("\n### Gene IDs of all MT- genes:\n")
print(mt)
cat("\n\n### Gene HGNC symbols of all MT- genes:\n")
print(gene_info[gene_info$is_mt, "hgnc_symbol"])

# Calculate fraction of MT reads per cell:
cell_info$fractionMTreads <- colSums(umi_matrix[mt, ]) / colSums(umi_matrix)
# find outliers of MT-fraction
cell_info$is_mt_outlier <- scater::isOutlier(
  metric = cell_info$fractionMTreads,
  log = FALSE,
  nmads = nmads_fractionMT,
  type = "higher"
)

# Filter cells with hard threshold of fraction of MT-reads
cell_info$mt_higher_threshold <- cell_info$fractionMTreads > threshold_MT


###   Filter genes with Number Of Detected Genes (NODG) ####

cat("\n\n###   Filter cells with number of detected genes\n")
# Number of detected genes:
cell_info$NODG <- colSums(umi_matrix > 0)
cell_info$rank_nodg <- rank(-cell_info$NODG)

# identify cells that are outliers
cell_info$is_nodg_outlier <- scater::isOutlier(
  metric = cell_info$NODG,
  log = FALSE,
  type = "lower",
  nmads = nmads_NODG
)

# Filter cells with hard threshold of NODG
cell_info$nodg_lower_threshold <- cell_info$NODG < threshold_NODG


###   Summarize CELL filters in temporary "remove_cell" column   ####

cell_info$remove_cell <- rowSums(cell_info[, c(
  "is_doublet",
  "is_mt_outlier",
  "mt_higher_threshold",
  "is_nodg_outlier",
  "nodg_lower_threshold"
)]) > 0


###   APPLY cell filter to count matrix before gene filtering   ####

matrix_cells_removed <- umi_matrix[, !cell_info$remove_cell, drop = FALSE]

# check if any cells are left after filtering
cat("\n\n### NOTE: Only continue if more than 0 cells are left after filtering.\n")
if (dim(matrix_cells_removed)[2] == 0) {
  stop("### ERROR: No cells are left after filtering! Please check sample quality.")
}


###   FILTER GENES   ####

cat("\n\n\n\n#####          FILTER GENES          #####\n\n")

###   Filter non-protein-coding genes   ####

# check if non-protein-coding genes are filtered out
stopifnot(is.logical(opt$protein_coding_only))

if (opt$protein_coding_only) {
  ### Filtering so that only protein-coding genes are kept
  cat("\n\n###   Filter non-protein-coding genes\n")

  # Download ensembl data
  genomeVersion <- opt$genomeVersion
  mart_obj <- biomaRt::useEnsembl(biomart = "ensembl", host = "www.ensembl.org", dataset = "hsapiens_gene_ensembl")

  if (genomeVersion == "hg19") {
    print("Use genome version hg19.")
    mart_obj <- biomaRt::useEnsembl(biomart = "ensembl", host = "www.ensembl.org", dataset = "hsapiens_gene_ensembl", GRCh = 37)
  }
  print(mart_obj)

  # get table with biomart info of only protein-coding genes
  biomart_protein_coding <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
    filters = c("ensembl_gene_id", "biotype"),
    values = list(rownames(umi_matrix), "protein_coding"),
    mart = mart_obj,
    uniqueRows = TRUE
  )

  # get mask if non-protein coding filter was applied for each gene, for exporting the information on filtered genes
  gene_info$non_protein_coding <- !(gene_info$ensembl_gene_id %in% biomart_protein_coding$ensembl_gene_id)

  # some ensembl_gene_ids are not unique in the mart object
  cat("\n\n###   All gene IDs:", length(biomart_protein_coding$ensembl_gene_id), "\n")
  cat("\n###   Unique gene IDs:", length(unique(biomart_protein_coding$ensembl_gene_id)))

  # filter out genes where the mart table contains NA
  biomart_no_NA <- na.omit(biomart_protein_coding)
  # genes found in biomart table but contain NA
  gene_info$biomart_NA <- gene_info$ensembl_gene_id %!in% biomart_no_NA$ensembl_gene_id & gene_info$ensembl_gene_id %in% biomart_protein_coding$ensembl_gene_id

  number_genes_filtered_biomart <- sum(rowSums(gene_info[, c("is_mt", "non_protein_coding", "biomart_NA")]) > 0)
  number_genes_accepted_biomart <- sum(rowSums(gene_info[, c("is_mt", "non_protein_coding", "biomart_NA")]) == 0)
  cat("\n\n###   Genes filtered with biomart:", number_genes_filtered_biomart)
  cat("\n###   Genes accepted with biomart:", number_genes_accepted_biomart)
} else if (!opt$protein_coding_only) {
  cat("\n\n###   Non protein-coding genes are not filtered\n")
}


###   Filter genes that are expressed in less than minimum number of cells   ####

cat("\n\n###   Filter genes not expressed in enough cells\n")

# Identify genes expressed in less than minNumberCells
gene_info$too_few_cells_express <- rowSums(matrix_cells_removed > 0) < minNumberCells


###   Filter genes that encode for ribosomal proteins   ####

cat("\n\n###   Filter genes encoding for ribosomal proteins\n")
gene_info$encodes_ribo_protein <- grepl(
  x = gene_info$hgnc_symbol,
  pattern = "^(RPL|MRPL|RPS|MRPS)"
)


###   Summarize GENE filters in "remove_gene" column   ####

if (opt$protein_coding_only) {
  gene_info$remove_gene <- rowSums(gene_info[, c(
    "is_mt",
    "too_few_cells_express",
    "encodes_ribo_protein",
    "non_protein_coding",
    "biomart_NA"
  )]) > 0
} else if (!opt$protein_coding_only) {
  gene_info$remove_gene <- rowSums(gene_info[, c(
    "is_mt",
    "too_few_cells_express",
    "encodes_ribo_protein"
  )]) > 0
}


###   APPLY cell and gene filter to umi_matrix before iterative gene and cell filtering   ####

matrix_filtered_temp <- umi_matrix[!gene_info$remove_gene, !cell_info$remove_cell, drop = FALSE]


###   Loop to make sure filters are applied   ####
# for NODGs and number of cells expressing a gene

cat("\n\n###   Filter genes and cells iteratively\n")

# collect barcodes and gene IDs for genes and cells to be filtered
loop_NODG_filtered_barcodes <- character()
loop_filtered_genes <- character()

# starting conditions for while loop
dim_matrix_start <- as.integer(c(1, 1))
dim_matrix_end <- as.integer(c(2, 2))

# counter for loop to limit iterations
counter <- 0
max_iterations <- 1000

# loop to enforce NODGs filter for cells and
# number of cells expressing a gene for genes
while (all(dim_matrix_start != dim_matrix_end) && counter < max_iterations) {
  counter <- counter + 1

  dim_matrix_start <- dim(matrix_filtered_temp)

  NODG <- colSums(matrix_filtered_temp > 0)
  too_few_NODG <- NODG < 400
  table(too_few_NODG)
  stopifnot(names(too_few_NODG) == colnames(matrix_filtered_temp))
  loop_NODG_filtered_barcodes <- c(loop_NODG_filtered_barcodes, colnames(matrix_filtered_temp)[too_few_NODG])
  matrix_filtered_temp <- matrix_filtered_temp[, !too_few_NODG]

  number_of_cells_express <- rowSums(matrix_filtered_temp > 0)
  too_few_cells_express <- number_of_cells_express < minNumberCells
  table(too_few_cells_express)
  stopifnot(names(too_few_cells_express) == rownames(matrix_filtered_temp))
  loop_filtered_genes <- c(loop_filtered_genes, rownames(matrix_filtered_temp)[too_few_cells_express])
  matrix_filtered_temp <- matrix_filtered_temp[!too_few_cells_express, ]

  dim_matrix_end <- dim(matrix_filtered_temp)
}
cat("\n\n###   Number of iterations: ", counter, "\t")

# Include iterative filters in gene_info and cell_info
cell_info$nodg_lower_threshold_iterative <- cell_info$barcodes %in% loop_NODG_filtered_barcodes
gene_info$too_few_cells_express_iterative <- gene_info$ensembl_gene_id %in% loop_filtered_genes


###   Summarize GENE filters in "remove_gene" column   ####

gene_info$remove_gene <- NULL
if (opt$protein_coding_only) {
  gene_info$remove_gene <- rowSums(gene_info[, c(
    "is_mt",
    "too_few_cells_express",
    "encodes_ribo_protein",
    "non_protein_coding",
    "biomart_NA",
    "too_few_cells_express_iterative"
  )]) > 0
} else if (!opt$protein_coding_only) {
  gene_info$remove_gene <- rowSums(gene_info[, c(
    "is_mt",
    "too_few_cells_express",
    "encodes_ribo_protein",
    "too_few_cells_express_iterative"
  )]) > 0
}


###   Summarize CELL filters in "remove_cell" column   ####

cell_info$remove_cell <- NULL
cell_info$remove_cell <- rowSums(cell_info[, c(
  "is_doublet",
  "is_mt_outlier",
  "mt_higher_threshold",
  "is_nodg_outlier",
  "nodg_lower_threshold",
  "nodg_lower_threshold_iterative"
)]) > 0


###   Export information about filtered genes and cells   ####

# give out info about removed cells
cat("\n\n\n###   Total number of cells: ", length(cell_info$barcodes))
cat("\n###   Number cells removed: ", sum(cell_info$remove_cell))
cat("\n###   Percentage cells removed: ", signif((sum(cell_info$remove_cell) / length(cell_info$barcodes)) * 100, digits = 4))

cat("\n\n###   Doublets: ", sum(cell_info$is_doublet))
cat("\n###   MT outlier: ", sum(cell_info$is_mt_outlier))
cat("\n###   MT higher than threshold: ", sum(cell_info$mt_higher_threshold))
cat("\n###   NODG outlier: ", sum(cell_info$is_nodg_outlier))
cat("\n###   NODG lower than threshold: ", sum(cell_info$nodg_lower_threshold))
cat("\n###   NODG lower than threshold in iterative filter: ", sum(cell_info$nodg_lower_threshold_iterative))

# Write file with full cell info
txtname <- paste0(outdir, file_name, ".cell_info.txt")
write.table(cell_info, txtname, sep = "\t", row.names = F, col.names = T, quote = F)

# give out info about removed genes
cat("\n\n\n\n###   Total number of genes: ", length(gene_info$ensembl_gene_id))
cat("\n###   Number genes removed: ", sum(gene_info$remove_gene))

cat("\n\n###   MT genes: ", sum(gene_info$is_mt))
cat("\n###   Genes expressed in too few cells: ", sum(gene_info$too_few_cells_express))
cat("\n###   Genes encoding ribosomal proteins: ", sum(gene_info$encodes_ribo_protein))
cat("\n###   Genes expressed in too few cells in iterative filter: ", sum(gene_info$too_few_cells_express_iterative))
# only relevant if protein-coding filter was applied
if (opt$protein_coding_only) {
  cat("\n###   Non-protein-coding: ", sum(gene_info$non_protein_coding))
  cat("\n###   NA in biomart table: ", sum(gene_info$biomart_NA), "\n\n")
}

# Export table with info about filtered genes
txtname <- paste0(outdir, file_name, ".gene_info.txt")
write.table(gene_info, txtname, sep = "\t", row.names = F, col.names = T, quote = F)


###   APPLY cell and gene filter to umi_matrix   ####

matrix_filtered <- umi_matrix[!gene_info$remove_gene, !cell_info$remove_cell, drop = FALSE]


###   Filter gene and cell info   ####

gene_info_filtered <- gene_info[gene_info$ensembl_gene_id %in% rownames(matrix_filtered), ]

# cell barcodes
output_cell_names <- colnames(matrix_filtered)
# fraction read mapping to MT genes
fractionMTreads_out <- cell_info$fractionMTreads[cell_info$barcodes %in% output_cell_names]

# remove colnames and rownames from matrix_filtered matrix
rownames(matrix_filtered) <- NULL
colnames(matrix_filtered) <- NULL


###   Write output hdf5 file   ####

outfile <- paste0(outdir, sample, ".genes_cells_filtered.h5")
h5createFile(outfile)
h5createGroup(outfile, "cell_attrs")
h5createGroup(outfile, "gene_attrs")
h5write(gene_info_filtered$ensembl_gene_id, outfile, "gene_attrs/gene_ids")
h5write(gene_info_filtered$hgnc_symbol, outfile, "gene_attrs/gene_names")
h5write(output_cell_names, outfile, "cell_attrs/cell_names")
h5write(fractionMTreads_out, outfile, "cell_attrs/fractionMT")

# set chunk size for writing h5 to c(1000,1000) or if the matrix is smaller to its dimensions (default)
if (dim(matrix_filtered)[1] > 1000 && dim(matrix_filtered)[2] > 1000) {
  chunks <- c(1000, 1000)
} else {
  chunks <- dim(matrix_filtered)
}
cat("\n\n### hdf5 chunk size:", chunks, "\n\n")
h5createDataset(file = outfile, dataset = "raw_counts", dims = dim(matrix_filtered), chunk = chunks)
h5write(matrix_filtered, outfile, "raw_counts")


###   Plot filtered cells   ####

cell_info$log_library_size <- log2(colSums(umi_matrix))
cell_info$log2_nodg <- log2(cell_info$NODG)

cell_info$col <- rep("black", nrow(cell_info))
cell_info$col[cell_info$is_mt_outlier] <- "red"
cell_info$col[cell_info$mt_higher_threshold] <- "red"
cell_info$col[cell_info$is_nodg_outlier] <- "cyan"
cell_info$col[cell_info$nodg_lower_threshold] <- "cyan"

both_filters <- rowSums(cell_info[, c("is_mt_outlier", "mt_higher_threshold")]) > 0 & rowSums(cell_info[, c("is_nodg_outlier", "nodg_lower_threshold")]) > 0
cell_info$col[both_filters] <- "orange"

cell_info$col[cell_info$nodg_lower_threshold_iterative] <- "magenta"
cell_info$col[cell_info$is_doublet] <- "green"
cell_info$col <- factor(cell_info$col, levels = c("black", "cyan", "green", "orange", "red", "magenta"))

# for manual colour scale
my_cols <- c(
  "black" = "black",
  "cyan" = "cyan",
  "green" = "green",
  "orange" = "orange",
  "red" = "red",
  "magenta" = "magenta",
  "test" = "tomato3",
  "test2" = "turquoise3"
)

# have size of point that represents cell be dependent on library size of the cell
point_size <- 0.4 + (cell_info$log_library_size - min(cell_info$log_library_size)) / diff(range(cell_info$log_library_size))

plot_cells_filtered <- ggplot(cell_info, aes(x = log2_nodg, y = fractionMTreads, alpha = 0.5)) +
  geom_point(aes(colour = col), size = point_size) +
  scale_colour_manual(
    values = my_cols,
    limits = my_cols,
    labels = c(
      "dot size corresponds to library size",
      "NODGs too low",
      "doublets",
      "filtered both criteria",
      "fraction MT too high",
      "NODGs too low after gene filtering",
      "absolute threshold fraction MT",
      "absolute threshold NODG"
    )
  ) +
  xlab("Log2(Number of Detected Genes)") +
  ylab("Fraction of MT reads") +
  ggtitle("Cell filtering based on number of detected genes and fraction of reads mapping to MT- genes") +
  geom_hline(aes(yintercept = threshold_MT), color = "tomato3", show.legend = T) +
  geom_vline(aes(xintercept = log2(threshold_NODG)), color = "turquoise3", show.legend = F) +
  guides(
    colour = guide_legend(
      title = NULL,
      override.aes = list(
        colour = my_cols,
        shape = c(rep(16, 6), NA, NA),
        linetype = c(rep(0, 6), 1, 1),
        size = c(rep(3, 6), 1.3, 1.3)
      )
    ),
    size = "none",
    alpha = "none"
  ) +
  theme_bw(base_size = 11.5) +
  theme(
    legend.position = c(0.8, 0.86),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.key.size = unit(.4, "cm"),
    legend.text = element_text(size = 9),
    title = element_text(size = 9),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  )

plotname <- paste0(outdir, file_name, ".visualize_filtered_cells.png")
ggplot2::ggsave(filename = plotname, plot = plot_cells_filtered, width = 19, height = 14, units = "cm", dpi = 300)


###   Plot ranking of cells for NODG, with threshold   ####

p_rank_nodg <- ggplot(cell_info, aes(x = rank_nodg, y = NODG)) +
  geom_point() +
  xlab("Cell rank") +
  ggtitle("Cell ranking according to number of detected genes (NODG)") +
  geom_hline(aes(yintercept = threshold_NODG, color = "red"), show.legend = T) +
  scale_colour_manual(
    name = NULL,
    labels = "threshold for NODG filtering",
    values = "red"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = c(0.8, 0.9))

filename <- paste0(outdir, file_name, ".cell_ranking_nodgs.png")
ggplot2::ggsave(filename = filename, plot = p_rank_nodg, width = 23, height = 18, units = "cm", dpi = 300)
