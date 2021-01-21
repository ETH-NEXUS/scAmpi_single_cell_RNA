########################################

####    File name: filter_genes_and_cells.R
####    Author: Anne Richter
####    Created in October 2018
####    R Version: 3.5.1

########################################

# Filtering of non-coding genes was included into this script
# Description:
## File name: select_protein_coding_genes.R
## Author: Dr. Franziska Singer and Mustafa Anil Tuncel
## Date created: 22.03.2018
## R Version: 3.4.1

#######################################


# This script takes an hdf5 file of raw expression counts (right after cellranger)

# Filtering Steps:
# 1. Cells are filtered if number of detected genes (NODG) is too low.
# 2. Cells are filtered if the fraction of reads mapped to MT- genes is too high.
# 3. Genes are filtered so that only protein-coding genes are kept.
# 4. Genes are filtered so that MT- genes are removed
# 5. Genes are filtered so that genes encoding for ribosomal proteins are removed

# Filtering is based on the scRNA tutorial by Atul Sethi and Panagiotis Papasaikas on ECCB in Athens

# R libraries
library(rhdf5)
library(optparse)
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
suppressMessages(library(scater))

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
# convenience function for NOT %in%
"%!in%" <- function(x, y) !("%in%"(x, y))

# parse command line arguments
option_list <- list(
  make_option("--hdf5File", type = "character", help = "Path to hdf5 input file. It includes raw expression matrix, gene & cell attributes."),
  make_option("--sample", type = "character", help = "Sample name, prefix of all output files"),
  make_option("--doublet_barcodes", type = "character", help = "Path to text file with doublet barcodes"),
  make_option("--nmads_NODG", type = "character", help = "Number of median-absolute-deviations away from median required for a value to be called an outlier, e.g. 5"),
  make_option("--nmads_fractionMT", type = "character", help = "Number of median-absolute-deviations away from median required for a value to be called an outlier, e.g. 5"),
  make_option("--outDir", type = "character", help = "Full path to output directory"),
  make_option("--genomeVersion", type = "character", help = "Specify the genome annotation version, either hg19 or GRCh38 are supported. The default is GRCh38."),
  make_option("--threshold_NODG", type = "character", help = "Hard threshold that gives the minimum NODG a cell must have to be further processed. E.g. 250"),
  make_option("--threshold_fractionMT", type = "character", help = "Hard threshold that gives the maximum fraction of MT reads a cell can have and be further processed. E.g. 0.5")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

###   set variable for filtering:
# minimum number of cells that should express a gene for it to be kept
minNumberCells <- 20

################################
#####     READ IN DATA     #####
################################

# get file name of input hdf5 and sample name
file_name <- basename(opt$hdf5File)
cat("\n\nInput file:", file_name)
sample <- opt$sample
outdir <- opt$outDir
cat("\nOutput directory:", outdir)

# print information about input hdf5 file
cat("\n\nContent of input h5 file:\n")
print(h5ls(opt$hdf5File))

# get count matrix
umi_counts <- h5read(opt$hdf5File, "raw_counts")
# link cell and gene information to matrix
rownames(umi_counts) <- h5read(opt$hdf5File, "gene_attrs/gene_ids")
colnames(umi_counts) <- h5read(opt$hdf5File, "cell_attrs/cell_names")

# get vector with cell_names (barcodes).
input_cell_names <- h5read(opt$hdf5File, "cell_attrs/cell_names")

# create data frame with mapping of gene_ids to hgnc gene_names
ENS2HGNC <- as.data.frame(h5read(opt$hdf5File, "gene_attrs/gene_ids"))
names(ENS2HGNC) <- "ensembl_gene_id"
ENS2HGNC$hgnc_symbol <- h5read(opt$hdf5File, "gene_attrs/gene_names")
rownames(ENS2HGNC) <- h5read(opt$hdf5File, "gene_attrs/gene_ids")
print("head(ENS2HGNC):")
print(head(ENS2HGNC))



cat("\n\n#####                                  #####\n")
cat("#####          CELL FILTERING          #####\n")
cat("#####                                  #####\n")

### Filter doublets
cat("\n\n\n###   Filtering doublets   ###\n")
doublet_barcodes <- read.csv(opt$doublet_barcodes, header = FALSE, stringsAsFactors = FALSE)
doublet_barcodes <- doublet_barcodes$V1
cat("\n\n\nDoublet barcodes object:\n")
print(str(doublet_barcodes))
doublets_mask <- colnames(umi_counts) %in% doublet_barcodes
indices_doublets <- which(doublets_mask)
perc_doublets <- signif((length(indices_doublets) / length(input_cell_names)) * 100, digits = 2)
# print out numbers of filtered cells:
print("Total number of cells:")
print(length(input_cell_names))
print("Doublets removed:")
print(length(indices_doublets))
print("Percentage of cells removed as doublets:")
print(paste(perc_doublets, "%"))


### Filtering according to fraction of MT-reads
cat("\n\n\n#############     Filtering cells according to fraction of reads mapped to MT-genes     #############\n")
nmads_fractionMT <- as.integer(opt$nmads_fractionMT)

# Grep all mitochondrial genes:
mt <- as.vector(ENS2HGNC[rownames(umi_counts), 1][grep("^MT-", ENS2HGNC[rownames(umi_counts), 2])])
print("Number of MT- genes:")
print(length(mt))
print("Gene IDs of all MT- genes:")
print(mt)

# Calculate fraction of MT reads per cell:
fractionMTreads <- colSums(umi_counts[mt, ]) / colSums(umi_counts)
# find outliers in terms of MT-fraction
MToutliers <- isOutlier(metric = fractionMTreads, log = FALSE, nmads = nmads_fractionMT, type = "higher")
cat("\n\ntable(MToutliers):\n")
print(table(MToutliers))
high_MT <- which(MToutliers)
perc_high_MT <- signif((length(high_MT) / length(input_cell_names)) * 100, digits = 2)
# print out numbers of filtered cells:
print("Total number of cells:")
print(length(input_cell_names))
print("Cells filtered because fraction of reads mapped to MT- genes makes them outliers:")
print(length(high_MT))
print("Percentage of cells filtered because fraction of reads mapped to MT- genes makes them outliers:")
print(paste(perc_high_MT, "%"))

# Filter cells according to hard threshold of fraction of MT- reads
threshold_MT <- as.numeric(opt$threshold_fractionMT)
cat("\n\nthreshold_MT:\n")
print(threshold_MT)
higher_than_threshold <- fractionMTreads > threshold_MT
high_abs_MT <- which(higher_than_threshold)
perc_high_threshold_MT <- signif((length(high_abs_MT) / length(input_cell_names)) * 100, digits = 2)
# print out numbers of filtered cells:
print("Total number of cells:")
print(length(input_cell_names))
print("Cells filtered because fraction of reads mapped to MT- genes is higher than threshold:")
print(length(high_abs_MT))
print("Percentage of cells filtered because fraction of reads mapped to MT- genes is higher than threshold:")
print(paste(perc_high_threshold_MT, "%"))


### Filtering according to Number of detected genes (NODG)
cat("\n\n\n#############     Filtering cells according to number of detected genes     #############\n")
nmads_NODG <- as.integer(opt$nmads_NODG)
# Number of detected genes:
NODG <- colSums(umi_counts > 0)
# identify cells that are outliers
outlierNODG <- isOutlier(metric = NODG, log = FALSE, type = "lower", nmads = nmads_NODG)
cat("\n\ntable(outlierNODG):\n")
print(table(outlierNODG))
low_NODG <- which(outlierNODG)
# print out numbers of filtered cells
perc_low_NODG <- signif((length(low_NODG) / length(input_cell_names)) * 100, digits = 2)
print("Total number of cells:")
print(length(input_cell_names))
print("Cells filtered because NODG makes them outliers:")
print(length(low_NODG))
print("Percentage of cells filtered because NODG makes them outliers:")
print(paste(perc_low_NODG, "%"))

# Filter cells according to hard threshold of NODG
threshold_NODG <- as.numeric(opt$threshold_NODG)
print("threshold_NODG:")
print(threshold_NODG)
lower_than_threshold <- NODG < threshold_NODG
low_abs_NODG <- which(lower_than_threshold)
# print out numbers of filtered cells
perc_low_threshold_NODG <- signif((length(low_abs_NODG) / length(input_cell_names)) * 100, digits = 2)
print("Total number of cells:")
print(length(input_cell_names))
print("Cells filtered because NODG is lower than threshold:")
print(length(low_abs_NODG))
print("Percentage of cells filtered because NODG is lower than threshold:")
print(paste(perc_low_threshold_NODG, "%"))


### Apply filters combined
# Merge all filtered cells:
filtered_cells <- unique(c(low_NODG, high_MT, high_abs_MT, low_abs_NODG))
# get cell barcodes of cell that will be filtered out
names_filtered_cells <- unique(names(c(low_NODG, high_MT, high_abs_MT, low_abs_NODG)))
names_filtered_cells <- unique(c(names_filtered_cells, doublet_barcodes))
print("str(names_filtered_cells):")
print(str(names_filtered_cells))

cat("\n\nTotal number of cells filtered from dataset:\n")
print(length(names_filtered_cells))
both_filters <- intersect(low_NODG, high_MT)
perc_filtered_out <- signif((length(unique(c(low_NODG, high_MT, high_abs_MT, low_abs_NODG))) / length(input_cell_names)) * 100, digits = 2)

# Write file with the cell barcodes of the cells that are filtered out
txtname <- outdir %&% file_name %&% "." %&% nmads_fractionMT %&% "_nmads_fractionMT." %&% nmads_NODG %&% "_nmads_NODG.filtered_cells.txt"
write.table(names_filtered_cells, txtname, sep = "\t", row.names = F, col.names = F, quote = F)

# Remove filtered cells from the dataset:
clean_umi_counts <- umi_counts[, !colnames(umi_counts) %in% names_filtered_cells, drop = FALSE]
print("Total number of cells before filtering:")
print(length(input_cell_names))
print("Number of cells after filtering:")
print(length(colnames(clean_umi_counts)))
print("Percentage of cells filtered out:")
print(paste(perc_filtered_out, "%"))

###   Show in plot ranking of cells for NODG, with threshold
# Number of detected genes:
plotname <- outdir %&% file_name %&% ".cell_ranking_nodgs.png"
png(plotname, width = 2200, height = 1800, res = 300)
# Plot NODGs ordered by rank (rank-size distribution)
plot (rank(-NODG), NODG,  pch = 19, xlab = "Cell rank", main = "Cell ranking according to number of detected genes (NODG)", cex.main = 1)
legend("topright", legend = c("threshold for NODG filtering"),
  lty = 1, col = c("red"), cex = 0.8)
# Threshold cells with low NODG:
abline(threshold_NODG, 0, col = "red")
dev.off()

# check if any cells are left after filtering
cat("\n\nOnly continue if more than 0 cells are left after filtering.\n")
if (dim(clean_umi_counts)[2] == 0) {
  stop("No cells are left after filtering! Please check sample quality.")
}

cat("\n\n\n#####                                  #####\n")
cat("#####          GENE FILTERING          #####\n")
cat("#####                                  #####\n")

print("Total number of genes before filtering:")
print(length(rownames(clean_umi_counts)))

### Filtering so that only protein-coding genes are kept
cat("\n\n#############     Filtering genes out if not protein coding     #############\n")

# Download ensembl data
genomeVersion <- opt$genomeVersion
mart_obj <- useEnsembl(biomart = "ensembl", host = "www.ensembl.org", dataset = "hsapiens_gene_ensembl")

if (genomeVersion == "hg19") {
  print("Use genome version hg19.")
  mart_obj <- useEnsembl(biomart = "ensembl", host = "www.ensembl.org", dataset = "hsapiens_gene_ensembl", GRCh = 37)
}
print(mart_obj)

entrezGeneMapping_proteinCoding <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "description"),
  filters = c("ensembl_gene_id", "biotype"), values = list(rownames(umi_counts), "protein_coding"), mart = mart_obj, uniqueRows = T)
print("str(entrezGeneMapping_proteinCoding):")
print(str(entrezGeneMapping_proteinCoding))
# get mask if non protein coding filter was applied for each gene, for exporting the information on filtered genes
non_prcoding_mask <- !(rownames(umi_counts) %in% entrezGeneMapping_proteinCoding$ensembl_gene_id)
cat("\n\ntable(non_prcoding_mask):\n")
print(table(non_prcoding_mask))
# see that some ensembl_gene_ids are not unique in the mart object
cat("\nAll gene IDs:", length(entrezGeneMapping_proteinCoding$ensembl_gene_id), "\n")
cat("\nUmique gene IDs:", length(unique(entrezGeneMapping_proteinCoding$ensembl_gene_id)), "\n")

# filter NAs and multiple entries
# filter NAs
entrezGeneMapping_proteinCoding_noNa <- na.omit(entrezGeneMapping_proteinCoding)
print("str(entrezGeneMapping_proteinCoding_noNa):")
print(str(entrezGeneMapping_proteinCoding_noNa))
# get mask if non na filter was applied for each gene, for exporting the information on filtered genes
na_mask <- rownames(umi_counts) %!in% entrezGeneMapping_proteinCoding_noNa$ensembl_gene_id & rownames(umi_counts) %in% entrezGeneMapping_proteinCoding$ensembl_gene_id
cat("\n\ntable(na_mask):\n")
print(table(na_mask))

# filter duplicated ensembl gene IDs
entrezGeneMapping_proteinCoding_noNa_unique <- entrezGeneMapping_proteinCoding_noNa[!duplicated(entrezGeneMapping_proteinCoding_noNa$ensembl_gene_id), ]
print("str(entrezGeneMapping_proteinCoding_noNa_unique), after removing duplicate ensembl IDs from mart table:")
print(str(entrezGeneMapping_proteinCoding_noNa_unique))

print(paste("Mapped to entrez and kept only protein coding genes. Remaining #genes:", dim(entrezGeneMapping_proteinCoding_noNa_unique)[1], sep = " "))

# creates a boolean (logical) mask vector
res <- rownames(clean_umi_counts) %in% entrezGeneMapping_proteinCoding_noNa_unique$ensembl_gene_id

cat("\n\n\nNumber of genes that are protein coding, do not contain NAs and are no duplicates in ensembl IDs:\n")
print(sum(res))

# Apply filter for protein coding genes to matrix and gene_id vector
clean_umi_counts <- clean_umi_counts[res, ]

# Print out numbers
cat("\n\nNumber of genes after protein-coding gene filtering:", length(rownames(clean_umi_counts)), "\n")


###   Filter genes that are expressed in less than the minimum number of cells set at the top of the script.
cat("\n\n\n#############     Filtering genes not expressed in enough cells     #############\n")
# Identify genes expressed in less than minNumberCells
absent_genes <- rowSums(clean_umi_counts > 0) < minNumberCells
cat("\n\ntable(absent_genes):")
print(table(absent_genes))
names_genes <- row.names(clean_umi_counts)
names_absent_genes <- names_genes[absent_genes]

# remove absent genes from dataset:
clean_umi_counts <- clean_umi_counts[!absent_genes, ]

# Print out numbers about filtered genes
perc_absent <- signif((sum(absent_genes) / length(ENS2HGNC$ensembl_gene_id)) * 100, digits = 2)
print("Number of genes that are not detected in enough of the remaining cells:")
print(sum(absent_genes))
print("Percentage of total number of genes that are not detected in any of the remaining cells:")
print(paste(perc_absent, "%"))
print("Number of genes after absent gene filtering:")
print(length(rownames(clean_umi_counts)))

# get absent gene mask for exporting the information on filtered genes
absent_genes_mask <- rownames(umi_counts) %in% names_absent_genes
cat("\n\ntable(absent_genes_mask):\n")
print(table(absent_genes_mask))


###   Filter genes that are MT- genes
cat("\n\n#############     Filtering out MT-genes     #############\n")
# Identify MT- genes
mt_genes <- rownames(clean_umi_counts) %in% mt

# remove MT- genes from dataset
clean_umi_counts <- clean_umi_counts[!mt_genes, ]

# Print out numbers about filtered genes
perc_MT <- signif((sum(mt_genes) / length(ENS2HGNC$ensembl_gene_id)) * 100, digits = 2)
print("Number of MT- genes:")
print(sum(mt_genes))
print("Percentage of total number of genes that are MT- genes:")
print(paste(perc_MT, "%"))
print("Number of genes after MT- gene filtering:")
print(length(rownames(clean_umi_counts)))

# get mask of MT genes for exporting the information on filtered genes
MT_mask <- rownames(umi_counts) %in% mt
cat("\n\ntable(MT_mask):\n")
print(table(MT_mask))


### Filter genes that encode for ribosomal proteins
cat("\n\n#############     Filtering out genes encoding for ribosomal proteins     #############\n")
# get gene_ids of genes encoding for ribosomal proteins
riboall <- as.vector(ENS2HGNC[rownames(clean_umi_counts), 1][grep("^(RPL|MRPL|RPS|MRPS)", ENS2HGNC[rownames(clean_umi_counts), 2])])
# get gene names
ribo_names <- as.vector(ENS2HGNC[rownames(clean_umi_counts), 2][grep("^(RPL|MRPL|RPS|MRPS)", ENS2HGNC[rownames(clean_umi_counts), 2])])

# Identify genes encoding for ribosomal proteins:
ribo_genes <- rownames(clean_umi_counts) %in% riboall

# remove ribosomal genes from dataset:
clean_umi_counts <- clean_umi_counts[!ribo_genes, ]

# Print out numbers about filtered genes
perc_ribo <- signif((sum(ribo_genes) / length(ENS2HGNC$ensembl_gene_id)) * 100, digits = 2)
print("Number of ribosomal genes:")
print(sum(ribo_genes))
print("Percentage of total number of genes that are ribosomal genes:")
print(paste(perc_ribo, "%"))
print("Number of genes after ribosomal gene filtering:")
print(length(rownames(clean_umi_counts)))

# get mask of ribosomal proteins for exporting the information on filtered genes
ribo_mask <- rownames(umi_counts) %in% riboall
cat("\n\ntable(ribo_mask):\n")
print(table(ribo_mask))



### Write file with the information on filtered genes
filtered_genes_out <- as.data.frame(rownames(umi_counts))
names(filtered_genes_out) <- "gene_ids"
filtered_genes_out$gene_names <- ENS2HGNC$hgnc_symbol
filtered_genes_out$contains_NAs <- na_mask
# filtered_genes_out$duplicated_hgnc <- duplicated_hgnc_mask
filtered_genes_out$non_protein_coding <- non_prcoding_mask
filtered_genes_out$absent_all_filtered_cells <- absent_genes_mask
filtered_genes_out$MT_genes <- MT_mask
filtered_genes_out$ribosomal_protein_genes <- ribo_mask

txtname <- outdir %&% file_name %&% ".information_on_genes_filtered.txt"
write.table(filtered_genes_out, txtname, sep = "\t", row.names = F, col.names = T, quote = F)


# filter cells again in case they do not have counts for any of the remaining genes
cat("\n\n#############     Filtering out cells that express none of the remaining genes     #############\n")
cells_filtered_second <- colSums(clean_umi_counts) == 0
cat("\n\ntable(cells_filtered_second):\n")
print(table(cells_filtered_second))
indic_filtered_second <- which(cells_filtered_second)
names_filtered_second <- names(indic_filtered_second)

clean_umi_counts <- clean_umi_counts[, !cells_filtered_second]
print("Number of cells before second filtering:")
print(length(cells_filtered_second))
print("Number of cells after second filtering:")
print(length(colnames(clean_umi_counts)))
print("Number of cells filtered out in second round:")
print(length(names_filtered_second))

if (length(names_filtered_second) > 0) {
  # Write file with the cell barcodes of the cells that are filtered out
  txtname <- outdir %&% file_name %&% ".filtered_cells_second_round.txt"
  write.table(names_filtered_second, txtname, sep = "\t", row.names = F, col.names = F, quote = F)
}

### Show in plot which cells were filtered:
# Log transformed umi counts:
Log_library_size <- log2(colSums(umi_counts))
# Point size proportional to library size :
point.size <- 0.25 + (Log_library_size - min(Log_library_size)) / diff(range(Log_library_size))
# Set a different color for the filtered cells:
col <- rep("black", ncol(umi_counts))
col[high_MT] <- "red"
col[higher_than_threshold] <- "red"
col[low_NODG] <- "cyan"
col[lower_than_threshold] <- "cyan"
col[both_filters] <- "orange"
col[indic_filtered_second] <- "magenta"
col[higher_than_threshold & lower_than_threshold] <- "orange"
col[doublets_mask] <- "green"
color_transparent <- adjustcolor(col, alpha.f = 0.5)
# Plot the fraction of MT reads as a function of the number of detected genes
plotname <- outdir %&% file_name %&% ".visualize_filtered_cells.png"
png(plotname, width = 2000, height = 1800, res = 300)
plot(log2(colSums(umi_counts > 0)), colSums(umi_counts[mt, ]) / colSums(umi_counts), pch = 19, cex = point.size, col = color_transparent, lwd = 0,
  xlab = "Log2(Number of Detected Genes)", ylab = "Fraction of MT reads",
  main = "Cell filtering based on Number of detected genes and fraction of reads mapping to MT- genes", cex.main = 0.8)
abline(h = threshold_MT, col = "red")
abline(v = log2(threshold_NODG), col = "cyan")
legend("topright", legend = c("fraction_MT_too_high", "NODG_too_low", "filtered_both_criteria", "no genes expressed after gene filtering", "dot size corresponds to library size", "absolute threshold fractionMT", "absolute threshold NODG", "doublets"),
  pch = c(19, 19, 19, 19, 19, 45, 45), col = c("red", "cyan", "orange", "magenta", "black", "red", "cyan", "green"), cex = 0.8)
dev.off()



######       Export cleaned dataset as hdf5 file        ########

cat("\n\n\nNumber of gene_names and gene_ids before filtering:", length(ENS2HGNC$hgnc_symbol), "\n")

# get filtered gene names that are not linked to matrix as rownames
gene_attrs_filtered <- ENS2HGNC[ENS2HGNC$ensembl_gene_id %in% rownames(clean_umi_counts), ]
cat("\nNumber of gene_names and gene_ids after filtering:", length(gene_attrs_filtered$hgnc_symbol), "\n")

# get remaining cell barcodes
output_cell_names <- colnames(clean_umi_counts)
cat("\nNumber of cell barcodes after filtering:", length(output_cell_names), "\n")

# add MT fraction to h5 output
fractionMTreads_out <- fractionMTreads[names(fractionMTreads) %in% output_cell_names]
# Make sure MT fractions are in order
stopifnot(all.equal(names(fractionMTreads_out), output_cell_names))
names(fractionMTreads_out) <- c()

# remove colnames and rownames from clean_umi_counts matrix
rownames(clean_umi_counts) <- c()
colnames(clean_umi_counts) <- c()

cat("\n\nCheck if rownames and colnames are NULL\n")
print(head(rownames(clean_umi_counts)))
print(head(colnames(clean_umi_counts)))

print("dim(clean_umi_counts) after removing row- and colnames")
print(dim(clean_umi_counts))
print("str(clean_umi_counts)")
print(str(clean_umi_counts))


# write the output hdf5 file
outfile <- outdir %&% sample %&% ".genes_cells_filtered.h5"
h5createFile(outfile)
h5createGroup(outfile, "cell_attrs")
h5createGroup(outfile, "gene_attrs")
h5write(gene_attrs_filtered$ensembl_gene_id, outfile, "gene_attrs/gene_ids")
h5write(gene_attrs_filtered$hgnc_symbol, outfile, "gene_attrs/gene_names")
h5write(output_cell_names, outfile, "cell_attrs/cell_names")
h5write(fractionMTreads_out, outfile, "cell_attrs/fractionMT")

# set chunk size for writing h5 to c(1000,1000) or if the matrix is smaller to its dimensions (default)
if (dim(clean_umi_counts)[1] > 1000 & dim(clean_umi_counts)[2] > 1000) {
  chunks <- c(1000,1000)
} else {
  chunks <- dim(clean_umi_counts)
}
cat("\n\nhdf5 chunk size:", chunks, "\n\n")
h5createDataset(file = outfile, dataset = "raw_counts", dims = dim(clean_umi_counts), chunk = chunks)
h5write(clean_umi_counts, outfile, "raw_counts")

