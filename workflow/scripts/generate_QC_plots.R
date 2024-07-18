########################################
####    File name: generate_QC_plots.R
####    Author: Anne Richter
####    Created in October 2018
####    R Version: 4.0
#########################################

# This script takes an hdf5 file of sc-RNA expression data and returns several QC plots
# Plots are saved in the same directory where the input file is
# The plots are based on the scRNA tutorial on ECCB in Athens

library(rhdf5)
library(optparse)
library(plyr)
library(ggplot2)
library(RColorBrewer)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")

# parse command line arguments
option_list <- list(
  make_option("--hdf5File", type = "character", help = "Path to hdf5 input file. It includes expression matrix, gene & cell attributes"),
  make_option("--sample_name", type = "character", help = "Sample name."),
  make_option("--sample_status", type = "character", help = "Sample status, if sample is filtered or unfiltered."),
  make_option("--outdir", type = "character", help = "Output directory.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


outdir <- opt$outdir
print(outdir)

# get file name of input hdf5 and sample name
file_name <- paste0(opt$outdir, opt$sample_name, ".", opt$sample_status)
print(file_name)


# read in hdf5 file with expression matrix
# link cell and gene information to matrix
umi_counts <- h5read(opt$hdf5File, "raw_counts")
rownames(umi_counts) <- h5read(opt$hdf5File, "gene_attrs/gene_ids")
colnames(umi_counts) <- h5read(opt$hdf5File, "cell_attrs/cell_names")
dim(umi_counts)
print(head(umi_counts[, 1:3]))

# create data frame with mapping of gene_ids to hgnc gene_names
ENS2HGNC <- as.data.frame(h5read(opt$hdf5File, "gene_attrs/gene_ids"))
names(ENS2HGNC) <- "ensembl_gene_id"

ENS2HGNC$hgnc_symbol <- h5read(opt$hdf5File, "gene_attrs/gene_names")

rownames(ENS2HGNC) <- h5read(opt$hdf5File, "gene_attrs/gene_ids")
print(head(ENS2HGNC))
print(dim(ENS2HGNC))


# Plot of gene detection probability ####
# (fraction of cells where a gene is detected) as a function of gene expression (total gene count)
plotname <- file_name %&% ".transcript_capture_efficiency.png"
png(plotname, width = 2200, height = 1800, res = 300)
plot1 <- smoothScatter(log2(rowSums(umi_counts) + 1), rowSums(umi_counts > 0) / ncol(umi_counts),
  nrpoints = 0,
  xlab = expression(Log[2] ~ "Total gene count"), ylab = "Detection probability",
  main = "Transcript Capture efficiency
    Gene detection probability (fraction of cells where a gene is detected)
    as a function of gene expression (total gene count)", cex.main = 1
)
dev.off()


# Histogram of library sizes in the dataset ####
plotname <- file_name %&% ".histogram_library_sizes.png"
png(plotname, width = 2000, height = 1800, res = 300)
hist(colSums(umi_counts) / 1e6, xlab = "Library size (millions)", breaks = 20, col = "grey80", ylab = "Number of cells", main = "Library Sizes")
dev.off()


# Histograms with number of detected genes and number of dropout values ####
plotname <- file_name %&% ".histogram_nodgs_dropout.png"
png(plotname, width = 3000, height = 1800, res = 300)
par(mfrow = c(1, 2))
# Histogram of number of detected genes:
hist(colSums(umi_counts > 0), xlab = "Number of detected genes", breaks = 20, col = "grey80", ylab = "Number of cells", main = "NODGs")
# Histogram of number of dropout values:
hist(colSums(umi_counts == 0), xlab = "Number of dropout values", breaks = 20, col = "grey80", ylab = "Number of cells", main = "Dropouts")
dev.off()


# Plot genes accounting for the majority of reads #####

# Sort genes in decreasing order in terms of their expression
ReadsPerGene <- sort(rowSums(umi_counts), decreasing = TRUE)
# Cumulative fraction of reads for the sorted genes:
CumulFraction <- cumsum(ReadsPerGene) / sum(ReadsPerGene)

# Fraction of reads coming from the top N genes calculated per cell:
for (N in c(25, 50, 100)) {
  topN <- names(ReadsPerGene)[1:N]
  ReadFraction <- apply(umi_counts, 2, function(x) x[topN] / sum(x))
  # Percentage of reads coming from the top N genes:
  # (signif just rounds)
  f <- signif(CumulFraction[N] * 100, digits = 3)
  # Produce a boxplot for the fraction of reads coming from the top N genes:
  plotname <- file_name %&% ".top" %&% N %&% "_genes_majority_of_reads.png"
  png(plotname, width = 2200, height = (1800 * (N / 35)), res = 300)
  title <- paste("Top ", N, " genes account for ", f, "% of reads", sep = "")
  boxplot(ReadFraction[N:1, ],
    use.cols = FALSE, horizontal = TRUE, outline = FALSE,
    boxwex = 0.5, names = rev(ENS2HGNC[topN, 2]), col = "orange", main = title,
    las = 2, par(cex.axis = 0.6, cex.main = 0.9), xlab = "Fraction of Reads in Cell"
  )
  dev.off()
}

# Plot proportion of mitochondrial (MT) reads #####
# Grep all mitochondrial genes:
mt <- as.vector(ENS2HGNC[, 1][grep("^MT-", ENS2HGNC[, 2])])
print(length(mt))
print(mt)

# Calculate fraction of MT reads:
df <- as.data.frame(colSums(umi_counts > 0))
names(df) <- "x"
df$y <- colSums(umi_counts[mt, ]) / colSums(umi_counts)
quant90 <- quantile(df$y, probs = 0.9)
quant95 <- quantile(df$y, probs = 0.95)
quant98 <- quantile(df$y, probs = 0.98)


# Plot fraction of MT reads as function of number of detected genes ####
plotname <- file_name %&% ".fraction_MT_reads.png"
png(plotname, width = 2200, height = 1800, res = 300)
ggplot(data = df, aes(x = x, y = y)) +
  geom_point(shape = 1) +
  scale_x_continuous(name = "Number of detected genes") +
  scale_y_continuous(name = "Fraction MT reads") +
  geom_hline(aes(yintercept = quant90, linetype = "90 percent quantile"), color = "red") +
  geom_hline(aes(yintercept = quant95, linetype = "95 percent quantile"), color = "orange") +
  geom_hline(aes(yintercept = quant98, linetype = "98 percent quantile"), color = "blue") +
  scale_linetype_manual(
    name = "", values = c(1, 1, 1),
    guide = guide_legend(override.aes = list(color = c("red", "orange", "blue")))
  ) +
  ggtitle("Fraction of reads mapping to MT- transcripts
as function of number of detected genes")
dev.off()


# Plot histogram of fraction of MT reads ####
plotname <- file_name %&% ".histogram_fraction_MT_reads.png"
png(plotname, width = 2200, height = 1800, res = 300)
plot_mt_fraction <- ggplot(data = df, aes(df$y)) +
  geom_histogram(stat = "bin", binwidth = 0.002) +
  xlab("Fraction of MT- reads") +
  ylab("Number of cells") +
  ggtitle("Fraction of reads mapping to transcripts of mitochondrial (MT-) genes")
plot_mt_fraction
dev.off()


# Fraction of reads mapping to genes encoding for ribosomal proteins ####

# Calculate fraction of ribosomal reads
riboall <- as.vector(ENS2HGNC[, 1][grep("^(RPL|MRPL|RPS|MRPS)", ENS2HGNC[, 2])])
print(length(riboall))
print(head(riboall))

ribo_names <- as.vector(ENS2HGNC[, 2][grep("^(RPL|MRPL|RPS|MRPS)", ENS2HGNC[, 2])])
print(length(ribo_names))
print(head(ribo_names))

fraction_riboall <- as.data.frame(colSums(umi_counts[riboall, ]) / colSums(umi_counts))
names(fraction_riboall) <- "fraction_ribosomal_genes_reads"
fraction_riboall$number_detected_genes <- colSums(umi_counts > 0)
head(fraction_riboall)

# Plot histogram of fraction of reads mapping to ribosomal protein transcripts ####
plotname <- file_name %&% ".histogram_fraction_ribosomal_reads.png"
png(plotname, width = 2200, height = 1800, res = 300)
ggplot(data = fraction_riboall, aes(fraction_riboall$fraction_ribosomal_genes_reads)) +
  geom_histogram(stat = "bin", binwidth = 0.003) +
  ggtitle("Fraction of reads mapping to transcripts of genes encoding for ribosomal proteins") +
  xlab("Fraction of (RPL|MRPL|RPS|MRPS) reads") +
  ylab("Number of cells")
dev.off()


# Plot the fraction of ribosomal reads as a function of number of detected genes:
plotname <- file_name %&% ".fraction_ribosomal_reads.png"
png(plotname, width = 2200, height = 1800, res = 300)
ggplot(data = fraction_riboall, aes(x = number_detected_genes, y = fraction_ribosomal_genes_reads)) +
  geom_point(shape = 1) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "top") +
  scale_x_continuous(name = "Number of detected genes") +
  scale_y_continuous(name = "Fraction of ribosomal reads") +
  ggtitle("Fraction of reads mapping to transcripts of genes
encoding for ribosomal proteins")
dev.off()


# Genes that are detected in at least one cell:
DetectedGenes <- which(rowSums(umi_counts) > 0)
# calculate the genes' mean expression (with pseudocounts):
mean_GE <- rowMeans(umi_counts[DetectedGenes, ] + 1 / ncol(umi_counts))
# calculate the genes' coefficient of variation for:
gene_cv <- apply(umi_counts[DetectedGenes, ], 1, function(x) sd(x) / mean(x + 1 / length(x)))
# Log transform expression and cv:
X1 <- log2(mean_GE)
Y1 <- log2(gene_cv + 1 / ncol(umi_counts))
# linear fit of log(cv) as a function of log(gene expression):
m <- lm(Y1 ~ X1)
plotname <- file_name %&% ".mean-variance_trend.png"
png(plotname, width = 2200, height = 1800, res = 300)
# scatterplot of log(cv) as a function of log(mean expression)
plot(X1, Y1,
  xlab = "log2(mean gene expression)", ylab = "log2(coefficent of variation)",
  main = "Gene dispersion as a function of their mean expression (Mean-variance trend)
     Counts are not normalized", pch = ".", cex = 2, col = "#00000055", cex.main = 1
)
# Add regression line
abline(coef(m)[1], coef(m)[2], col = "red", lwd = 2, lty = 2)
# Slope in m-v trend according to poisson distribution:
abline(0, -0.5, col = "grey", lwd = 2, lty = 2)
legend("topright", legend = c("Regression line", "Poisson distribution"), lty = c(2, 2), lwd = c(2, 2), col = c("red", "grey"))
dev.off()


# Ranking plots for cell filtering ####
# Number of detected genes:
NODG <- colSums(umi_counts > 0)
plotname <- file_name %&% ".nodgs_cell_ranking.png"
png(plotname, width = 2200, height = 1800, res = 300)
# Plot NODGs ordered by rank  (rank-size distribution)
plot(rank(-NODG), NODG, pch = 19, xlab = "Cell rank", main = "Cell ranking according to number of detected genes (NODG)", cex.main = 1)
# Threshold cells with low NODG:
# abline(1000,0,col="red")
dev.off()


# Fraction of MT reads per cell ####
fractionMTreads <- colSums(umi_counts[mt, ]) / colSums(umi_counts)
plotname <- file_name %&% ".cell_ranking_mt_reads.png"
png(plotname, width = 2200, height = 1800, res = 300)
df_fractionMTreads <- as.data.frame(fractionMTreads)
names(df_fractionMTreads) <- "fractionMTreads"
df_fractionMTreads$rank <- rank(fractionMTreads)
ggplot(data = df_fractionMTreads, aes(x = df_fractionMTreads$rank, y = df_fractionMTreads$fractionMTreads)) +
  geom_point(shape = 1) +
  scale_x_continuous(name = "Cell rank") +
  scale_y_continuous(name = "fractionMTreads", breaks = seq(0, 1, 0.1)) +
  ggtitle("Cell ranking according to fraction of MT- reads") +
  geom_hline(aes(yintercept = quant90, linetype = "90 percent quantile"), color = "red") +
  geom_hline(aes(yintercept = quant95, linetype = "95 percent quantile"), color = "orange") +
  geom_hline(aes(yintercept = quant98, linetype = "98 percent quantile"), color = "blue") +
  scale_linetype_manual(
    name = "", values = c(1, 1, 1),
    guide = guide_legend(override.aes = list(color = c("red", "orange", "blue")))
  )
dev.off()
