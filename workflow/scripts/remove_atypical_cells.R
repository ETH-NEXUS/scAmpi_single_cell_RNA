#################################################
####    File name: remove_atypical_cells.R
####    Author: Anne Bertolini
####    Created June 2019
####    Updated December 2020
####    R Version: 4.0
################################################

suppressPackageStartupMessages({
  library(rhdf5)
  library(optparse)
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(scran)
  library(limma)
  library(purrr)
  library(Hmisc)
})

h5disableFileLocking()
print("h5disableFileLocking() called")

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# parse command line arguments
option_list <- list(
  make_option("--sce_in", type = "character", help = "Path to RDS input file that contains sce object."),
  make_option("--cluster_table", type = "character", help = "Table with overview over cluster and the cell types of the cells in each cluster. '*.phenograph_celltype_association.txt'"),
  make_option("--celltype_config", type = "character", help = "Config text file with two columns where the Major cell types and the respective subtypes are specified."),
  make_option("--threshold_filter", type = "character", help = "Threshold number. With a percentage or number of cells (see argument --) below this number cells will be filtered out."),
  make_option("--min_threshold", type = "character", help = "Absolute minimum number of cells for the threshold. Threshold for cell number that has to be exceeded so that the respective cell type is not filtered out cannot be lower."),
  make_option("--threshold_type", type = "character", help = "Set if cells are filtered out based on cell number or percentage. Possible parameters are 'number_cells' and 'percentage_cells'"),
  make_option("--sample_name", type = "character", help = "Sample name that will be the prefix of all output files."),
  make_option("--outDir", type = "character", help = "Full path to output directory")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")

# read in thresholds for filtering
threshold_filter <- as.integer(opt$threshold_filter)
min_threshold <- as.integer(opt$min_threshold)
cat("\n\nthreshold_filter:", threshold_filter, "(", opt$threshold_type, ")\n")
cat("\n\nmin_threshold:", min_threshold, "cells\n\n")

# read in sce object from RDS file
my_sce <- readRDS(opt$sce_in)
cat("\n\nmy_sce:\n\n")
print(my_sce)

# read in cell type config file
celltype_config <- read.csv(opt$celltype_config, sep = "\t", stringsAsFactors = FALSE)
celltypes_overview <- as.data.frame(celltype_config$Major, stringsAsFactors = FALSE)
names(celltypes_overview) <- "Major"
celltypes_overview$Minor <- celltypes_overview$Major

# format table with all major (`Major`) and final (`Minor`) cell types
for (my_celltype in 1:dim(celltype_config)[1]) {
  Minor <- unlist(strsplit(celltype_config[my_celltype, 2], ",", fixed = TRUE))
  add_to_df <- as.data.frame(Minor, stringsAsFactors = FALSE)
  add_to_df$Major <- celltype_config[my_celltype, 1]
  print(add_to_df)
  celltypes_overview <- rbind(celltypes_overview, add_to_df)
}

# remove entries with "none"
celltypes_overview <- celltypes_overview[!celltypes_overview$Minor == "none", ]
# add 'uncertain' and 'unknown'
uncertain_unknown <- data.frame(c("uncertain", "unknown"), stringsAsFactors = FALSE)
names(uncertain_unknown) <- "Major"
uncertain_unknown$Minor <- uncertain_unknown$Major
celltypes_overview <- rbind(celltypes_overview, uncertain_unknown)
# shorten cell type names
celltypes_overview$Major <- gsub(pattern = "([^_]+)_.*", "\\1", celltypes_overview$Major)
celltypes_overview$Minor <- gsub(pattern = "([^_]+)_.*", "\\1", celltypes_overview$Minor)
# check if Minor cell types are unique
stopifnot(length(celltypes_overview$Minor) == length(unique(celltypes_overview$Minor)))
# have minor types as rownames to use the table as dictionary
# dictionary will be used to translate dominant cell type (can be minor) to its major cell type
rownames(celltypes_overview) <- celltypes_overview$Minor


# read in cell type composition (dominant) of phenograph clusters
ct_composition <- read.csv(opt$cluster_table, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

# get the dominant cell types of the clusters into data frame
# Here it is assumed that there is no cluster 0 present in the data!
dominant_overview <- as.data.frame(ct_composition$Dominant.celltype, stringsAsFactors = FALSE)
names(dominant_overview) <- "dominant"
# remove rows with empty string
dominant_overview <- dominant_overview[!dominant_overview$dominant == "", , drop = FALSE]
# replace "majorCT_mixed" with "majorCT"
dominant_overview$dominant_celltype <- gsub(pattern = "([^_]+)_mixed", "\\1", dominant_overview$dominant)
# add cluster ID to dominant_overview
cluster_id <- ct_composition$Cluster
cluster_id <- cluster_id[!cluster_id %in% c("Sum", "Percent (95% ci)")]
dominant_overview$cluster_id <- cluster_id
# have the major type of the dominant cell type in separate column
# use rownames of celltypes_overview to retrieve the respective Major type (dictionary)
dominant_overview$major_of_dominant <- celltypes_overview[dominant_overview$dominant_celltype, "Major"]


# Have table with cell information from the SCE object (colData)
df_colData_cells <- as.data.frame(colData(my_sce), stringsAsFactors = FALSE)


# initiate data frame that will contain the "atypical" cells
colData_atypical <- setNames(data.frame(matrix(ncol = length(names(df_colData_cells)), nrow = 0)),
                             names(df_colData_cells))


# iterate over all clusters and identify cells to be removed
for (cluster in dominant_overview$cluster_id) {
  cluster_index <- which(dominant_overview$cluster_id == cluster)
  cluster_dominant <- dominant_overview$dominant[cluster_index]
  major_of_dominant <- dominant_overview$major_of_dominant[cluster_index]

  # cells are removed only from clusters with clear cell type, not classified "mixed"
  if (cluster_dominant %in% c("mixed", "uncertain", "unknown")) {
    cat("\n\n\n###   Cluster ", cluster, "is '", cluster_dominant, "' and not filtered   ###")
  }
  else if (!(cluster_dominant %in% c("mixed", "uncertain", "unknown"))) {
    cluster <- as.integer(cluster)
    cat("\n\n\n###   Start filtering cluster", cluster, "   ###\n")
    df_colData_cells$celltype_final <- as.character(df_colData_cells$celltype_final)
    df_colData_cells$celltype_major <- as.character(df_colData_cells$celltype_major)
    # colData of only cells of current cluster
    colData_cells.in.cluster <- subset(df_colData_cells, df_colData_cells$phenograph_clusters == cluster)
    cat("\nNumber of cells:", length(colData_cells.in.cluster$barcodes), "\n")
    cat("Dominant cell type:", cluster_dominant, "\n")
    cat("Its major type:", major_of_dominant, "\n")

    # set threshold that is specific to cluster
    if (opt$threshold_type == "percentage_cells") {
      threshold_cluster <- (threshold_filter * length(colData_cells.in.cluster$barcodes)) / 100
      threshold_cluster <- round(threshold_cluster)
    } else if (opt$threshold_type == "number_cells") {
      threshold_cluster <- opt$threshold_filter
    } else {
      stop("The command line argument '--threshold_type' is not correctly specified! Please check the help message for details.")
    }
    # threshold must be at least min_threshold
    if (threshold_cluster < min_threshold) {
      threshold_cluster <- min_threshold
      cat("Threshold was set to min_threshold\n")
    }
    cat("Threshold for cluster", cluster, ":", threshold_cluster, "\n\n")

    # find all cell types present in cluster
    all_cts_in_cluster <- unique(as.character(colData_cells.in.cluster$celltype_major))
    cat("Cell types found in cluster", cluster, ":\n")
    print(all_cts_in_cluster)
    all_unrelated_types <- all_cts_in_cluster[!all_cts_in_cluster == major_of_dominant]
    cat("Unrelated cell types filtered out if cell number is not above threshold:\n")
    print(all_unrelated_types)
    cat("\n")

    # iterate over unrelated cell types and filter out respective cells if cell number is not above threshold.
    for (my_celltype in all_unrelated_types) {
      subset_cells_type <- colData_cells.in.cluster[colData_cells.in.cluster$celltype_major == my_celltype, , drop = FALSE]
      number_of_cells <- length(subset_cells_type$barcode)
      if (number_of_cells < threshold_cluster) {
        cat(my_celltype, "removed:", number_of_cells, "\n")
        # add information about cells that will be removed to data frame
        colData_atypical <- rbind(colData_atypical, subset_cells_type)
      } else {
        cat(my_celltype, ": not removed because", number_of_cells, "is above filtering threshold\n")
      }
    }
  }
}


# remove cells from sce object that are filtered out
# get binary mask for removing the cells
mask_remove_cells <- colData(my_sce)$barcodes %in% colData_atypical$barcodes
number_removed_cells <- sum(mask_remove_cells)
cat("\n\nNumber of cells that are filtered out:", number_removed_cells, "\n")

# filter assays and colData of sce object
my_sce_filtered <- my_sce[, !mask_remove_cells]

# remove filtered cells from the second dimension of the phenograph distance matrix
reducedDims(my_sce_filtered)$phenodist <- reducedDims(my_sce_filtered)$phenodist[, !mask_remove_cells]

# Write filtered sce object into RDS file
filename_out <- opt$outDir %&% opt$sample_name %&% ".atypical_removed.RDS"
cat("\nFile name of output RDS file with filtered sce object:\n", filename_out)
saveRDS(my_sce_filtered, filename_out)

# Write table with information about atypical cells
colData_atypical_file <- opt$outDir %&% opt$sample_name %&% ".atypical_cells_colData.tsv"
cat("\n\nFile name of output table with atypical cells:\n", colData_atypical_file, "\n")
write.table(colData_atypical, colData_atypical_file, sep = "\t", row.names = F, quote = FALSE)


# write output hdf5 file
outfile <- opt$outDir %&% opt$sample_name %&% ".atypical_removed.h5"
# create file
h5createFile(outfile)
# create groups
h5createGroup(outfile, "cell_attrs")
h5createGroup(outfile, "gene_attrs")

# write gene attributes
h5write(rowData(my_sce_filtered)$gene_ids, outfile, "gene_attrs/gene_ids")
h5write(rowData(my_sce_filtered)$gene_names, outfile, "gene_attrs/gene_names")
h5write(rowData(my_sce_filtered)$detection_rate.x, outfile, "gene_attrs/detection_rate.x")
h5write(rowData(my_sce_filtered)$gmean, outfile, "gene_attrs/gmean")
h5write(rowData(my_sce_filtered)$variance, outfile, "gene_attrs/variance")
h5write(rowData(my_sce_filtered)$residual_mean, outfile, "gene_attrs/residual_mean")
h5write(rowData(my_sce_filtered)$residual_variance, outfile, "gene_attrs/residual_variance")
h5write(rowData(my_sce_filtered)$SYMBOL, outfile, "gene_attrs/SYMBOL")
h5write(rowData(my_sce_filtered)$mean, outfile, "gene_attrs/mean")
h5write(rowData(my_sce_filtered)$detection_rate.y, outfile, "gene_attrs/detection_rate.y")
h5write(rowData(my_sce_filtered)$var, outfile, "gene_attrs/var")

# write cell attributes
h5write(colData(my_sce_filtered)$barcodes, outfile, "cell_attrs/cell_names")
h5write(colData(my_sce_filtered)$g2m_score, outfile, "cell_attrs/g2m_score")
h5write(colData(my_sce_filtered)$s_score, outfile, "cell_attrs/s_score")
h5write(colData(my_sce_filtered)$barcodes, outfile, "cell_attrs/barcodes")
h5write(colData(my_sce_filtered)$fractionMT, outfile, "cell_attrs/fractionMT")
h5write(colData(my_sce_filtered)$n_umi, outfile, "cell_attrs/n_umi")
h5write(colData(my_sce_filtered)$n_gene, outfile, "cell_attrs/n_gene")
h5write(colData(my_sce_filtered)$log_umi, outfile, "cell_attrs/log_umi")
h5write(colData(my_sce_filtered)$cycle_phase, outfile, "cell_attrs/cycle_phase")
h5write(colData(my_sce_filtered)$umap_cl, outfile, "cell_attrs/umap_cl")
h5write(colData(my_sce_filtered)$phenograph_clusters, outfile, "cell_attrs/phenograph_clusters")
h5write(as.character(colData(my_sce_filtered)$celltype_major_full_ct_name), outfile, "cell_attrs/celltype_major_full_ct_name")
h5write(as.character(colData(my_sce_filtered)$celltype_major), outfile, "cell_attrs/celltype_major")
h5write(as.character(colData(my_sce_filtered)$celltype_final_full_ct_name), outfile, "cell_attrs/celltype_final_full_ct_name")
h5write(as.character(colData(my_sce_filtered)$celltype_final), outfile, "cell_attrs/celltype_final")

# set chunk size for writing h5 to c(1000,1000) or if the matrix is smaller to its dimensions (default)
if (dim(assay(my_sce_filtered, "counts"))[1] > 1000 & dim(assay(my_sce_filtered, "counts"))[2] > 1000) {
  chunks <- c(1000,1000)
} else {
  chunks <- dim(assay(my_sce_filtered, "counts"))
}
cat("\n\nhdf5 chunk size:", chunks, "\n\n")
# write assays
h5createDataset(file = outfile, dataset = "raw_counts", dims = dim(assay(my_sce_filtered, "counts")), chunk = chunks)
h5write(assay(my_sce_filtered, "counts"), outfile, "raw_counts")
h5createDataset(file = outfile, dataset = "normcounts", dims = dim(assay(my_sce_filtered, "normcounts")), chunk = chunks)
h5write(assay(my_sce_filtered, "normcounts"), outfile, "normcounts")
h5createDataset(file = outfile, dataset = "pearson_resid", dims = dim(assay(my_sce_filtered, "pearson_resid")), chunk = chunks)
h5write(assay(my_sce_filtered, "pearson_resid"), outfile, "pearson_resid")


################################################################################
# Write out overview table over cell types and phenograph clusters
# Code is copied from the celltyping.r script written by Michael Prummer
################################################################################

f.round.preserve.sum <- function(x, digits = 0) {
  up <- 10^digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

################################################################################
## compare cell type with phenograph_clusters
res <- colData(my_sce_filtered)
tab <- table(res$phenograph_clusters, res$celltype_final)
top <- ncol(tab)
left <- nrow(tab)
tab <- addmargins(tab)
tab <- rbind(tab, f.round.preserve.sum(tab[nrow(tab), ] / tab[nrow(tab), ncol(tab)] * 100))
rownames(tab)[nrow(tab)] <- "Percent"

#########################################
# read in celltype config file, subtypes are individual for each major type
type_config <- as.data.frame(read.table(opt$celltype_config, sep = "\t", head = T, stringsAsFactors = F))
major_types <- type_config$Major
print("str(major_types):")
print(str(major_types))
minor_types <- lapply(type_config$Subtype, function(x) strsplit(x, ",")[[1]])
names(minor_types) <- type_config$Major
# remove "none"
minor_types <- lapply(minor_types, function(x) x[x != "none"])
# remove NA
idx <- sapply(minor_types, function(x) length(which(is.na(x))))
minor_types <- minor_types[which(idx == 0)]
# remove empty list items
idx <- sapply(minor_types, length)
minor_types <- minor_types[idx > 0]
# final list
print("str(minor_types):")
print(str(minor_types))

#########################################
## begin new cluster purity calculations
# if at least one minor_type, make major-minor dictionary:
all.major.types <- sort(c(metadata(my_sce)$all_major_types, "uncertain", "unknown"))
if (length(minor_types) > 0) {
  minor_types <- lapply(minor_types, function(x) gsub("([^_]+)_.*", "\\1", x))
  names(minor_types) <- gsub("([^_]+)_.*", "\\1", names(minor_types))
  # add major type also to minor type column (as it can be final cell type as well)
  minor_types <- sapply(names(minor_types), simplify = FALSE, function(x) {
    minor_types[x] <- c(minor_types[[x]], x)
  })
  major_minor_dict <- melt(minor_types)
  names(major_minor_dict) <- c("minor", "major")
  # also have all major types in the dictionary
  major_minor_dict$minor <- as.character(major_minor_dict$minor)
  for (m_type in all.major.types) {
    if (!m_type %in% major_minor_dict$minor) {
      major_minor_dict[nrow(major_minor_dict) + 1, ] <- c(m_type, m_type)
    }
  }
}
# loop over all phenograph clusters. Need not be continuous or numeric but unique
clu.pu <- data.frame("Cluster" = rownames(tab[seq(left), ]), "Dominant.celltype" = NA,
  "Celltype composition" = NA, stringsAsFactors = F, check.names = F)
for (ii in unique(res$phenograph_clusters)) {
  # ii = unique(res$phenograph_clusters)[1]
  rowid <- which(rownames(tab) == ii)
  # calculate percentage of each cell type in cluster ii
  purity <- tab[rowid, seq(top)] / sum(tab[rowid, seq(top)])
  purity <- f.round.preserve.sum(purity * 100)
  # a "dominant.type" is a cell type that contributes more than 20% to the cluster
  dominant.type <- which(purity > 20)
  if (length(dominant.type) == 1) {
    # if only one dominant cell type is found in cluster ii insert information into
    # the column "Dominant.celltype"
    clu.pu$Dominant.celltype[rowid] <- names(dominant.type)
    ct.comp <- names(dominant.type) %&% " (" %&% purity[dominant.type] %&% "%)"
    clu.pu$`Celltype composition`[rowid] <- ct.comp
  } else {
    if (length(dominant.type) > 1) {
      clu.pu$Dominant.celltype[rowid] <- "mixed"
      # order dominant celltype labels by decreasing proportion
      ct.order <- order(purity[dominant.type], decreasing = T)
      ct.comp.1 <- names(dominant.type)[ct.order]
      ct.comp.2 <- "(" %&% purity[dominant.type[ct.order]] %&% "%)"
      clu.pu$`Celltype composition`[rowid] <- paste(ct.comp.1, ct.comp.2, collapse = ", ")
      # if all dominant.types have the same major_celltype, change to "mayor_type_mixed"
      if (exists("major_minor_dict")) {
        idx <- match(names(dominant.type), major_minor_dict$minor)
        idx <- idx[!is.na(idx)]
        id.major <- major_minor_dict$major[idx]
        # only if all dominant types have the same major_type
        if (length(id.major) > 0 & length(unique(id.major)) == 1) {
          dominant_major <- unique(id.major)
          clu.pu$Dominant.celltype[rowid] <- dominant_major %&% "_mixed"
        }
      }
    } else {
      clu.pu$Dominant.celltype[rowid] <- "mixed"
      clu.pu$`Celltype composition`[rowid] <- "all<20"
    }
  }
}
clu.pu <- rbind(clu.pu, c("Sum", "", ""), c("Percent (95% ci)", "", ""))
tmp <- as.data.frame(as.matrix(tab))
tmp$Cluster <- c(sort(unique(res$phenograph_clusters)), "Sum", "Percent (95% ci)")
clu.pu <- merge(tmp, clu.pu, sort = F)
# # confidence intervals for proportions (Wilson)
# round(binconf(0, 4446, method="wilson")*100,1)
tab.ci <- binconf(tab[nrow(tab) - 1, seq(ncol(tab) - 1)],
  rep(tab[nrow(tab) - 1, ncol(tab)], ncol(tab) - 1),
  method = "wilson") * 100
tab.ci <- apply(tab.ci, 2, function(x) as.character(signif(x, 2)))
tab.ci <- apply(tab.ci, 1, function(x) sprintf("%s (%s;%s)", x[1], x[2], x[3]))
clu.pu[nrow(clu.pu), 1 + seq(length(tab.ci))] <- tab.ci
# to remove confidence intervals, comment out the previous 4 lines.
ct_association_filename <- opt$outDir %&% opt$sample_name %&% ".atypical_removed.phenograph_celltype_association.txt"
write.table(clu.pu, file = ct_association_filename, sep = "\t", row.names = F)
## end new cluster purity calculations
#########################################

# print info about cells filtered out
cat("\n\n##########     Total number of cells filtered out:      ##########\n\n")
cat(length(df_colData_cells$barcodes), "total cells in sample\n")
cat(length(colData_atypical$barcodes), "cells filtered out\n")
percentage_filtered_out <- round((length(colData_atypical$barcodes) / length(df_colData_cells$barcodes)) * 100, 2)
cat(percentage_filtered_out, "percent of cells filtered out\n")
cat(length(df_colData_cells$barcodes) - length(colData_atypical$barcodes), "cells kept\n")
atypical_major <- as.data.frame(table(colData_atypical$celltype_major))
atypical_major <- atypical_major[order(atypical_major$Freq, decreasing = TRUE), ]
cat("\n\n Numbers major cell types of atypical cells:\n")
print(atypical_major)

# print out percentages of cell types, major and final
perc_major <- as.data.frame(round(prop.table(table(res$celltype_major)), 3) * 100)
perc_major <- perc_major[order(perc_major$Freq, decreasing = TRUE), ]
cat("\n\n Frequencies major cell types:\n")
print(perc_major)


perc_final <- as.data.frame(round(prop.table(table(res$celltype_final)), 3) * 100)
perc_final <- perc_final[order(perc_final$Freq, decreasing = TRUE), ]
cat("\n\n Frequencies final cell types:\n")
print(perc_final)

cat("\n\n Number of cells: ", length(res$barcodes))

cat("\n\n Number of clusters: ", length(unique(res$phenograph_clusters)), "\n\n")

new_dom_types <- clu.pu$Dominant.celltype[!clu.pu$Dominant.celltype == ""]
dom_types <- as.data.frame(table(new_dom_types))
dom_types <- dom_types[order(dom_types$Freq, decreasing = TRUE), ]
cat("\n\n Frequencies dominant cell types of clusters:\n")
print(dom_types)

