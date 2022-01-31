# ########################################

####    File name: get_cluster_gene_expression.R
####    Author: Anne Bertolini
####    Created in March 2019
####    R Version: 3.5.1

# ########################################

# This script takes as input:
# 1) RDS file that contains a SCE object
# 2) a tab separated file with priority genes (e.g. .../required_files/melanoma/drug_targets_extended.tsv)
# 3) the filtering threshold for genes per sample
# 4) the filtering threshold for genes per cluster
# 5) the filtering type for the sample
# 6) the filtering type for the cluster
# 7) and the output directory

# It returns three tab separated tables:
# 1) one with one row per gene and two columns per cluster, one column for average expression values and one column for percentage of cells that are non-zero together with the rank of expression, separated by "/". Two additional columns for gene_ids and gene_names.
# 2) a subset of output file 1) containing only the priority genes
# 3) one file with one row per rank of expression and one column per cluster that contains the gene_names and mean expression values separated by "/" and ordered in descending order. One additional column for the rank of the genes in this row.


suppressPackageStartupMessages({
library(optparse)
library(plyr)
library(dplyr)
library(ggplot2)
library(purrr)
})

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# get division string to have a better overview for printed out results
division_string <- "           -------------------------------------          "

# parse command line arguments
option_list = list(
make_option("--sce_in", type = "character", help = "Path to RDS input file that contains a SCE object. It includes among other data the expression matrix normalised and corrected for cell cycle and gene & cell attributes."),
make_option("--sample_name", type = "character", help = "Sample name."),
make_option("--priority_genes", type = "character", help = "Tab separated file with header, first column contains gene names."),
make_option("--filtering_threshold_sample", type = "character", help = "Filtering threshold to not include genes if either number or percentage of cells with non-zero counts is below threshold. Filter is applied per sample. 
For percentage or total number see argument --filter_type. Example value: 20"),
make_option("--filter_type_sample", type = "character", help = "Argument if genes are filtered with an absolute number of cells or a percentage of cells that have to have non-zero counts for this gene to be included. Concerns filter per sample. Either 'number_cells' or 'percentage_cells'"),
make_option("--outDir", type = "character", help = "Full path to output directory")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# read in priority genes
priority_genes <- read.csv(opt$priority_genes, header=TRUE, stringsAsFactors=FALSE, sep = "\t")
print("str(priority_genes):")
print(str(priority_genes))
print("head(priority_genes):")
print(head(priority_genes))


# read in data from sce object in RDS file
my_sce <- readRDS(opt$sce_in)
print("my_sce:")
print(my_sce)
# count matrix
count_matrix = assay(my_sce, "normcounts")
# cell table
all_cell_attrs <- colData(my_sce)
row.names(all_cell_attrs) <- all_cell_attrs$barcodes
all_cell_attrs$phenograph_clusters <- as.integer(all_cell_attrs$phenograph_clusters)
print("head(all_cell_attrs):")
print(head(all_cell_attrs))
# gene table
all_gene_attrs <- rowData(my_sce)
print("head(all_gene_attrs):")
print(head(all_gene_attrs))

# create a named list to use as a dict
# gene ids as keyL and gene names as value
lookup = as.list(rowData(my_sce)$gene_names)
lookup = setNames(lookup, rowData(my_sce)$gene_ids)
print("head(lookup):")
print(head(lookup))

# assign the rownames to keep track of the gene IDs
row.names(count_matrix) = rowData(my_sce)$gene_ids
colnames(count_matrix) = colData(my_sce)$barcodes
print("str(count_matrix):")
print(str(count_matrix))


# Filter out genes with a too small number of cells with counts in the whole sample
print(paste0(division_string, "Filter genes with filtering threshold ", division_string))
print("dim(count_matrix):")
print(dim(count_matrix))
counts_nonzero <- count_matrix > 0
print("counts_nonzero[1:5]:")
print(counts_nonzero[1:5])
# number of cells with non-zero counts per gene
rowsums_counts_nonzero <- rowSums(counts_nonzero)
print("head(rowsums_counts_nonzero):")
print(head(rowsums_counts_nonzero))

# retrieve filtering threshold for genes
print("class(opt$filter_type_cluster[1]):")
print(class(opt$filter_type_cluster[1]))

threshold_sample <- NULL

if (opt$filter_type_sample[1] == "number_cells"){
        threshold_sample <- as.numeric(opt$filtering_threshold_sample)
} else if (opt$filter_type_sample[1] == "percentage_cells"){
        threshold_sample <- ((as.numeric(opt$filtering_threshold_sample))*length(colnames(count_matrix)))/100
} else {
        stop("The command line argument '--filter_type_sample' is not correctly specified! Please check the help message for details.")
}


# Filter out cells from cluster_0:
print(paste0(division_string, "Filter out cells from cluster_0 ", division_string))

# From cluster data filter out cells with cluster ID 0
all_cell_attrs = all_cell_attrs[!all_cell_attrs$phenograph_clusters == "0",]
print("str(all_cell_attrs) without cells of cluster 0: ")
print(str(all_cell_attrs))

# From matrix filter out cells with cluster ID 0
count_matrix = count_matrix[, colnames(count_matrix) %in% all_cell_attrs$barcodes]
print("str(count_matrix):")
print(str(count_matrix))
print("str(colnames(count_matrix)):")
print(str(colnames(count_matrix)))


# continue filtering out genes
print(paste0(division_string, "Filter genes with filtering threshold ", division_string))
# get mask for genes with counts above the filtering threshold
print("Filtering type sample: ")
print(opt$filter_type_sample[1])
print("threshold_sample:")
print(threshold_sample)
mask_filtering <- rowsums_counts_nonzero >= threshold_sample
print("str(mask_filtering):")
print(str(mask_filtering))

number_genes_considered <- sum(mask_filtering)
print("number_genes_considered:")
print(number_genes_considered)
print("number_genes_total:")
print(length(rownames(count_matrix)))

# Filter count_matrix
count_matrix <- count_matrix[mask_filtering,]
print("dim(count_matrix) after filtering genes sample-wide:")
print(dim(count_matrix))


# get IDs of all cell clusters
all_clusters = sort(unique(all_cell_attrs$phenograph_clusters))
all_clusters <- c(all_clusters, "total_sample")
print("str(all_clusters):")
print(str(all_clusters))
print("all_clusters:")
print(all_clusters)

# initiate lists to fill with data per cluster
list_mean <- list()
list_ranked <- list()

# iterate over each cluster
for (cluster in all_clusters){
    print("#######################################     NEXT CLUSTER    #########################################")
    print(paste0("Cluster ", cluster))

    # Get subset of the count matrix with only cells from current cluster
    # get subset of all_cell_attrs
    if (cluster == "total_sample"){
        # get all cell attrs
        sub_all_cell_attrs <- all_cell_attrs
        # get whole matrix
        sub_matrix <- count_matrix
    }else{
        # get subset of cell attrs
        sub_all_cell_attrs <- all_cell_attrs[all_cell_attrs$phenograph_clusters == cluster,]
        # get sub_matrix
        sub_matrix <- count_matrix[, colnames(count_matrix) %in% sub_all_cell_attrs$barcodes]
    }
    print("head(sub_all_cell_attrs):")
    print(head(sub_all_cell_attrs))
    cells_in_cluster <- as.numeric(length(sub_all_cell_attrs$barcodes))
    print(paste0("Number of cells in cluster ", cluster, ":"))
    print(cells_in_cluster)

    # mask of cells that have non-zero count
    counts_nonzero <- sub_matrix > 0
    print("str(counts_nonzero):")
    print(str(counts_nonzero))

    # count number of cells that have non-zero counts per gene
    rowsums_counts_nonzero <- rowSums(counts_nonzero)
    print("str(rowsums_counts_nonzero):")
    print(str(rowsums_counts_nonzero))

    # calculate mean expression
    cluster_subtable <- as.data.frame(apply(sub_matrix, 1, mean, na.rm=TRUE))
    names(cluster_subtable) <- "mean_expr"
    cluster_subtable$gene_ids <- row.names(sub_matrix)
    print("str(cluster_subtable):")
    print(str(cluster_subtable))
    print("summary(cluster_subtable):")
    print(summary(cluster_subtable))

    # gene names in additional column
    cluster_subtable$gene_names = as.character(lapply(cluster_subtable$gene_ids, function(x) {lookup[toString(x)][[1]]}))
    print("sapply(cluster_subtable, class):")
    print(sapply(cluster_subtable, class))
    print("str(cluster_subtable):")
    print(str(cluster_subtable))

    # sort by mean expression
    # two same values are given the rank of the latter value (e.g. ranks 1,3,3 for three values 1.5, 1.0, 1.0)
    cluster_subtable$rank <- rank(-cluster_subtable$mean_expr, na.last = "keep", ties.method = "max")
    print("head(cluster_subtable):")
    print(head(cluster_subtable))
    # round values
    cluster_subtable$mean_expr <- round(cluster_subtable$mean_expr, digits = 3)
    # Have column with percentage of cells in this cluster with non-zero counts
    cluster_subtable$pct_nonzero <- round((rowsums_counts_nonzero/cells_in_cluster)*100, digits = 1)
    # get column with percentage of cells expressing gene and rank
    cluster_subtable$`pct_nonzero/rank` <- paste(cluster_subtable$pct_nonzero, cluster_subtable$rank, sep = "/")
    # get column with gene name and mean expression
    cluster_subtable$`gene_name/mean_expr` <- paste(cluster_subtable$gene_names, cluster_subtable$mean_expr, sep = "/")
    print("str(cluster_subtable):")
    print(str(cluster_subtable))
    print("head(cluster_subtable):")
    print(head(cluster_subtable))

    # have cluster ID in column names
    colnames(cluster_subtable) <- gsub("^mean_expr", paste0("cl_", cluster, "_mean_expr"), colnames(cluster_subtable), perl = TRUE)
    colnames(cluster_subtable) <- gsub("pct_nonzero", paste0("cl_", cluster, "_pct_nonzero"), colnames(cluster_subtable))
    colnames(cluster_subtable) <- gsub("gene_name/mean_expr", paste0("cl_", cluster, "_gene_name/mean_expr"), colnames(cluster_subtable))
    print("head(cluster_subtable):")
    print(head(cluster_subtable))

    # get data frame that will be put in list per cluster, then be merged for all clusters and given out as tsv
    mean_out <- as.data.frame(select(cluster_subtable, c("gene_ids", paste0("cl_", cluster, "_mean_expr"), paste0("cl_", cluster, "_pct_nonzero/rank"))))
    print("head(mean_out):")
    print(head(mean_out))

    # add data frame to list
    list_mean[[cluster]] <- mean_out

    # get rank data out
    ranked_out <- as.data.frame(select(cluster_subtable, c("rank", paste0("cl_",cluster, "_gene_name/mean_expr"))))
    print("head(ranked_out):")
    print(head(ranked_out))
 
    ranked_out <- ranked_out[order(ranked_out$rank, decreasing = FALSE, na.last = TRUE),]
    ranked_out <- as.data.frame(select(ranked_out, c(paste0("cl_", cluster, "_gene_name/mean_expr"))))
    ranked_out$rank <- seq(1,length(ranked_out[,1]), by = 1)
    print("head(ranked_out):")
    print(head(ranked_out))
    print("tail(ranked_out):")
    print(tail(ranked_out))

    # add data frame to list
    list_ranked[[paste0("cl_", cluster, "_gene_name/mean_expr")]] <- ranked_out

    ranked_out <- ranked_out[order(ranked_out$rank, decreasing = FALSE, na.last = TRUE),]
    ranked_out <- as.data.frame(select(ranked_out, c(paste0("cl_", cluster, "_gene_name/mean_expr"))))
    print("head(ranked_out):")
    print(head(ranked_out))
    ranked_out$rank <- seq(1,length(ranked_out[,1]), by = 1)
    print("head(ranked_out):")
    print(head(ranked_out))
}

print(paste0(division_string, "After looping over clusters merge data sets together ", division_string))

# structure of list of data frames with mean expressions
print("str(list_mean):")
print(str(list_mean))

# check if all data frames in the list_mean have the same number of cols and rows
dims_dfs_mean <- lapply(list_mean, dim)
print("str(dims_dfs_mean):")
print(str(dims_dfs_mean))
number_rows_mean <- sapply(dims_dfs_mean, "[", 1)
number_cols_mean <- sapply(dims_dfs_mean, "[", 2)
uni_rows_mean <- unique(number_rows_mean)
uni_cols_mean <- unique(number_cols_mean)
print("str(uni_rows_mean):")
print(str(uni_rows_mean))
print("str(uni_cols_mean):")
print(str(uni_cols_mean))
stopifnot(length(uni_rows_mean) == 1)
stopifnot(length(uni_cols_mean) == 1)

# structure of list of data frames with ranked expressions
print("str(list_ranked):")
print(str(list_ranked))
# check if all data frames in the list_mean have the same number of cols and rows
dims_dfs_ranked <- lapply(list_ranked, dim)
print("str(dims_dfs_ranked):")
print(str(dims_dfs_ranked))
number_rows_ranked <- sapply(dims_dfs_ranked, "[", 1)
number_cols_ranked <- sapply(dims_dfs_ranked, "[", 2)
uni_rows_ranked <- unique(number_rows_ranked)
uni_cols_ranked <- unique(number_cols_ranked)
print("str(uni_rows_ranked):")
print(str(uni_rows_ranked))
print("str(uni_cols_ranked):")
print(str(uni_cols_ranked))
stopifnot(length(uni_rows_ranked) == 1)
stopifnot(length(uni_cols_ranked) == 1)

# merge data frames in list_mean to one data frame
final_mean_df <- list_mean %>% purrr::reduce(inner_join, by = "gene_ids")
print("table(final_mean_df < 0.0):")
print(table(final_mean_df < 0.0))
final_mean_df[final_mean_df < 0.0] <- 0.0
print("table(final_mean_df < 0.0):")
print(table(final_mean_df < 0.0))

print("length(final_mean_df[,1]):")
print(length(final_mean_df[,1]))
stopifnot(length(final_mean_df[,1]) == uni_rows_mean)
print("str(final_mean_df):")
print(str(final_mean_df))
# get gene names into table
final_mean_df$gene_names = as.character(lapply(final_mean_df$gene_ids, function(x) {lookup[toString(x)][[1]]}))
print("str(final_mean_df):")
print(str(final_mean_df))
# reformat to have gene_names column first
final_mean_df <- final_mean_df %>% select("gene_names", everything())
# have "total_sample" column name without "cluster"
colnames(final_mean_df) <- gsub("cl_total_sample", "total_sample", colnames(final_mean_df))
print("str(final_mean_df):")
print(str(final_mean_df))

# get final_mean_table for priority genes only
priority_final_mean_df <- final_mean_df[final_mean_df$gene_names %in% priority_genes$SYMBOL,]
print("str(priority_final_mean_df):")
print(str(priority_final_mean_df))

# merge data frames in list_ranked to one data frame
final_ranked_df <- list_ranked %>% purrr::reduce(inner_join, by = "rank")
print("length(final_ranked_df[,1]):")
print(length(final_ranked_df[,1]))
stopifnot(length(final_ranked_df[,1]) == uni_rows_ranked)
print("str(final_ranked_df):")
print(str(final_ranked_df))
# reformat to have rank column first
final_ranked_df <- final_ranked_df %>% select("rank", everything())
# have "total_sample" column name without "cluster"
colnames(final_ranked_df) <- gsub("cl_total_sample", "total_sample", colnames(final_ranked_df))
print("tail(final_ranked_df) before replacing values below 0:")
print(tail(final_ranked_df))
print("str(final_ranked_df):")
print(str(final_ranked_df))
# set all mean values below 0 to 0
final_ranked_df <- apply(final_ranked_df, 2, function(x) gsub(pattern = "([a-zA-Z0-9_-]+/)-\\d+.?\\d*",
                                                              replacement = "\\10", x, perl = TRUE))
final_ranked_df <- as.data.frame(final_ranked_df, stringsAsFactors = FALSE)
final_ranked_df[,1] <- as.integer(final_ranked_df[,1])
print("tail(final_ranked_df) after replacing values below 0:")
print(tail(final_ranked_df))
print("str(final_ranked_df):")
print(str(final_ranked_df))


# write data frame with mean expressions into output file
output_file_name_mean <- paste(opt$outDir, opt$sample_name, ".gene_expression_per_cluster.tsv", sep="")
print(output_file_name_mean)
write.table(final_mean_df, file = output_file_name_mean, quote = F, sep = "\t", row.names = F)

# write data frame with mean expressions of priority genes into output file
output_file_name_mean_priority <- paste(opt$outDir, opt$sample_name, ".gene_expression_per_cluster_priority_genes.tsv", sep="")
print(output_file_name_mean_priority)
write.table(priority_final_mean_df, file = output_file_name_mean_priority, quote = F, sep = "\t", row.names = F)

# write data frame with ranked results into output file
output_file_name_ranked <- paste(opt$outDir, opt$sample_name, ".gene_expression_per_cluster_ranked.tsv", sep="")
print(output_file_name_ranked)
write.table(final_ranked_df, file = output_file_name_ranked, quote = F, sep = "\t", row.names = F)
