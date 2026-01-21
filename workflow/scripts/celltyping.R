################################################################################
## cell type classification - perform cell typing by comparing gene expression
## profiles with pre-defined lists
################################################################################

################################################################################
## Load packages and parse command line options

# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
# convenience function for loading packages
f_pck_load <- function(x) {
  suppressWarnings(
    suppressMessages(
      require(x, character.only = T, warn.conflicts = F, quietly = T)
    )
  )
}
# Load packages
lby <- c(
  "optparse",
  "reshape2",
  "Hmisc",
  "SingleCellExperiment",
  "scROSHI"
)
resp <- lapply(lby, f_pck_load)
if (!all(unlist(resp))) {
  print(resp)
  stop("Could not load one or more packages")
} else {
  rm(resp, lby)
}

# command line arguments are parsed
option_list <- list(
  make_option("--SCE", type = "character", help = "Path to sce onject file with input data (sce_basic.RDS)."),
  make_option("--celltype_lists", type = "character", help = "Path to the file containing marker genes all cell types."),
  make_option("--celltype_config", type = "character", help = "Path to the file containing a config file that specifies major and sub cell types."),
  make_option("--min_genes", type = "character", help = "Minimum number of celltype specific genes expressed. If less: classify unknown."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# give out session Info
cat("\n\n")
print(Sys.time())
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n")
cat("\nInput files:\n\n")
print(opt)
cat("\n\n")


path <- opt$outputDirec %&% opt$sampleName
min_genes <- as.numeric(opt$min_genes)
if (!is.finite(min_genes)) stop("min_genes was not provided correctly.")

# Set global options
options(stringsAsFactors = FALSE)
# options(error=function()traceback(2))

f.round.preserve.sum <- function(x, digits = 0) {
  up <- 10^digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}

################################################################################
## main code starts here  ####
################################################################################
## load input data ####
sce_data <- readRDS(opt$SCE)

# celltype config file, sub types are individual for each major type
type_config <- as.data.frame(read.table(opt$celltype_config, sep = "\t", head = T, stringsAsFactors = F))
major_types <- type_config$Major
cat("\n\nstr(major_types):\n")
str(major_types)
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
cat("\n\nstr(minor_types):\n")
str(minor_types)

##################################################### scROSHI ####

sce_data <- scROSHI(
  sce_data = sce_data,
  celltype_lists = opt$celltype_lists,
  type_config = type_config,
  count_data = "normcounts",
  gene_symbol = "SYMBOL",
  cell_scores = TRUE,
  min_genes = min_genes,
  min_var = 1.5,
  n_top_genes = 2000,
  n_nn = 5,
  thresh_unknown = 0.05,
  thresh_uncert = 0.1,
  thresh_uncert_second = 0.8
)

##################################################### scROSHI

# make sure celltype_final is factor including all final cell types as levels
minor_types_levels_long <- unique(unlist(minor_types))
minor_types_levels <- gsub("^([^_]*)_.*$",
                           replacement = "\\1",
                           x = minor_types_levels_long,
                           perl = T)
major_types_levels <- gsub("^([^_]*)_.*$",
                           replacement = "\\1",
                           x = major_types,
                           perl = T)

all_celltype_levels <- c(minor_types_levels,
                         major_types_levels,
                         "uncertain",
                         "unknown")

colData(sce_data)$celltype_final <- factor(colData(sce_data)$celltype_final,
                                           levels = all_celltype_levels)


# remove cells that are in phenograph cluster 0
# those are cells that were not similar enough to other cells to be
# assigned to a cluster
mask_keep_cells <- colData(sce_data)$phenograph_clusters != 0
cat("\n\nNumber of cells that remain after filtering out cluster 0:\n\n")
print(table(mask_keep_cells))

## write information about cells in cluster 0 into table
sce_data_cluster0 <- sce_data[, !mask_keep_cells]
info_cluster0 <- as.data.frame(colData(sce_data_cluster0))
write.table(info_cluster0, file = path %&% ".info_cells_cluster_0.tsv",
            sep = "\t",
            row.names = F)
rm(info_cluster0)
rm(sce_data_cluster0)

sce_data <- sce_data[, mask_keep_cells]
reducedDims(sce_data)$phenodist <- reducedDims(sce_data)$phenodist[, mask_keep_cells]


## write final celltype to disk
res <- colData(sce_data)[, c("barcodes", "celltype_final")]
write.table(res, file = path %&% ".cts_final.txt", sep = "\t", row.names = F)


################################################################################
## compare cell type with phenograph_clusters
res <- colData(sce_data)
tab <- table(res$phenograph_clusters, res$celltype_final)
top <- ncol(tab)
left <- nrow(tab)
tab <- addmargins(tab)
tab <- rbind(tab, f.round.preserve.sum(tab[nrow(tab), ] / tab[nrow(tab), ncol(tab)] * 100))
rownames(tab)[nrow(tab)] <- "Percent"


#########################################
## begin new cluster purity calculations ####
# if at least one minor_type, make major-minor dictionary:

all.major.types <- sort(c(metadata(sce_data)$all_major_types, "uncertain", "unknown"))

if (length(minor_types) > 0) {
  minor_types <- lapply(minor_types, function(x) gsub("([^_]+)_.*", "\\1", x))
  names(minor_types) <- gsub("([^_]+)_.*", "\\1", names(minor_types))
  # add major type also to minor type column (as it can be final cell type as well)
  minor_types <- sapply(names(minor_types), simplify = FALSE, function(x) {
    minor_types[x] <- c(minor_types[[x]], x)
  })
  major_minor_dict <- reshape2::melt(minor_types)
  names(major_minor_dict) <- c("minor", "major")
  # also have all major types in the dictionary
  for (m_type in all.major.types) {
    if (!m_type %in% major_minor_dict$minor) {
      major_minor_dict[nrow(major_minor_dict) + 1, ] <- c(m_type, m_type)
    }
  }
}
# loop over all phenograph clusters. Need not be continuous or numeric but unique
clu.pu <- data.frame(
  "Cluster" = rownames(tab[seq(left), ]), "Dominant.celltype" = NA,
  "Celltype composition" = NA, stringsAsFactors = F, check.names = F
)
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
# round(Hmisc::binconf(0, 4446, method="wilson")*100,1)
tab.ci <- Hmisc::binconf(tab[nrow(tab) - 1, seq(ncol(tab) - 1)],
  rep(tab[nrow(tab) - 1, ncol(tab)], ncol(tab) - 1),
  method = "wilson"
) * 100
tab.ci <- apply(tab.ci, 2, function(x) as.character(signif(x, 2)))
tab.ci <- apply(tab.ci, 1, function(x) sprintf("%s (%s;%s)", x[1], x[2], x[3]))
clu.pu[nrow(clu.pu), 1 + seq(length(tab.ci))] <- tab.ci
# to remove confidence intervals, comment out the previous 4 lines.
write.table(clu.pu, file = path %&% ".celltyping.phenograph_celltype_association.txt", sep = "\t", row.names = F)
## end new cluster purity calculations


# print out percentages of cell types, major and final ####
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

################################################################################
## save sce object for later use ####
saveRDS(sce_data, path %&% ".celltyping.RDS")
