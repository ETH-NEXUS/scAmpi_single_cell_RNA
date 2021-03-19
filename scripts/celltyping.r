
################################################################################
## cell type classification - perform cell typing by comparing gene expression
## profiles with pre-defined lists
################################################################################


################################################################################
## INIT: set some tuning parameters
min_var <- 1.5 # minimum variance for highly variable genes
n_top_genes <- 2000 # maximum number of highly variable genes
n_nn <- 5 # number of nearest neighbors for umap for assignment of celltypes
# to uncertain/unknowns (from major celltyping)
thresh_unknown <- 0.05
thresh_uncert <- 0.1
thresh_uncert_second <- 0.8

################################################################################
## Load packages and parse command line options

# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
# convenience function for loading packages
f_pck_load <- function(x) {
  suppressWarnings(
    suppressMessages(
      require(x, character.only = T, warn.conflicts = F, quietly = T))
  )
}
# Load packages
lby <- c("optparse", "reshape2", "scran", "limma", "uwot", "igraph", "Hmisc",
  "pheatmap", "RColorBrewer", "cowplot")
resp <- lapply(lby, f_pck_load)
if (!all(unlist(resp))) {
  print(resp)
  stop("Could not load one or more packages")
} else {
  rm(resp, lby)
}

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

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

path <- opt$outputDirec %&% opt$sampleName
min_genes <- as.numeric(opt$min_genes)
if (!is.finite(min_genes)) stop("min_genes was not provided correctly.")
print(opt$SCE)

# Set global options
options(stringsAsFactors = FALSE)
# options(error=function()traceback(2))


################################################################################
## Functions

# calculate wilcox p-value or return preset values
f.my.wilcox.test <- function(x, y, min_genes = 5) {
  if (length(x) >= min_genes & length(y) >= min_genes) {
    return(wilcox.test(x, y, alternative = "greater")$p.value)
  } else {
    return(1)
  }
}


# calculate celltype score: U-test
f.score.ctgenes.U <- function(sce, gset, min_genes = 5) {
  cat("\nCalculate scores for cell type classification.")
  # gset = cell.type[these.ct]
  # sce = these.cells
  # min_genes = 5
  # select count table with all.ct.genes
  all.genes <- unique(as.character(unlist(gset)))
  cat("\n\nNumber of genes on cell type specific lists:", length(all.genes))
  # count matrix including all cells and only the cell type specific genes
  mat <- assay(sce, "normcounts")[rowData(sce)$SYMBOL %in% all.genes, , drop = F]
  cat("\n\nNumber of genes included in matrix:", dim(mat)[1])
  # keep genes only if they have counts in at least one cell
  keep <- rowSums(mat)
  mat <- mat[keep > 0, , drop = F]
  cat("\n\nNumber of genes with with a sum of counts > 0:", dim(mat)[1], "\n\n\n")
  colnames(mat) <- sce$barcodes
  # get row index of dd for each gene in a list of celltypes
  idxs <- ids2indices(gene.sets = gset, rownames(mat), remove.empty = F)
  # generate gene x celltype matrix
  # all values are 0
  ds <- matrix(0, nrow = nrow(mat), ncol = length(idxs))
  rownames(ds) <- rownames(mat)
  colnames(ds) <- names(idxs)
  # fill matrix with 1 where a gene is specific for a cell type
  for (cell_type in seq(length(idxs))) {
    ds[idxs[[cell_type]], cell_type] <- 1
  }
  # perform wilcox test using gene x sample and gene x celltype matrices.
  # wilcox.test(mat[ds[,1]==1,1],mat[ds[,1]==0,1])$p.value
  m.cts <- apply(mat, 2, function(x) {
    # on each column of mat (each cell)
    # perform wilcox test with each list of cell type genes against all other cell type genes
    apply(ds, 2, function(y) f.my.wilcox.test(x[y == 1], x[y == 0], min_genes))
  })
  m.cts[is.na(m.cts)] <- 1
  return(m.cts)
}


# choose best matching celltype
f.annot.ctgenes <- function(m.cts, unknown, uncertain) {
  ct.annot <- apply(m.cts, 2, function(x) names(which.min(x)))
  ct.annot <- data.frame(barcodes = names(ct.annot),
    cell.type = as.character(ct.annot),
    stringsAsFactors = F)
  best <- apply(m.cts, 2, min)
  ct.annot$celltype.score <- -log10(best)
  # deal with low certainty annotations:
  secondbest <- apply(m.cts, 2, function(x) sort(x, decreasing = F)[2])
  ct.annot$cts.2 <- -log10(secondbest)
  ct.annot$reclass <- NA
  # if best and second best are too similar & secondbest is at all relevant
  # cell gets reclassified as "uncertain"
  ct.annot$reclass[best / secondbest > uncertain & secondbest < unknown] <- "uncertain"
  # deal with low confidence annotations:
  # if total cell type score is too low cell gets reclassified as "unknown"
  ct.annot$reclass[ct.annot$celltype.score < -log10(unknown)] <- "unknown"
  # cell.type col is overwritten with `reclass` col where necessary
  ct.annot$cell.type[!is.na(ct.annot$reclass)] <- ct.annot$reclass[!is.na(ct.annot$reclass)]
  return(ct.annot[, 1:2])
}


f.round.preserve.sum <- function(x, digits = 0) {
  up <- 10^digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}


# score heatmap plot function
f_plot_final_score <- function(final.score, major.type) {
  bad <- apply(final.score, 2, function(x) length(unique(x)))
  notbad <- which(bad != 1)
  if (length(notbad) > 1) {
    hm1 <- pheatmap(-log10(final.score[, notbad, drop = F]), scale = "none", color = col.pal,
      silent = T, show_rownames = T, show_colnames = F,
      clustering_method = "ward.D2", treeheight_row = 0, treeheight_col = 20)
    hm2 <- pheatmap(-log10(final.score[, notbad, drop = F]), scale = "column", color = col.pal,
      silent = T, show_rownames = T, show_colnames = F,
      clustering_method = "ward.D2", treeheight_row = 0, treeheight_col = 20)
    pp <- plot_grid(hm1$gtable, hm2$gtable, nrow = 2)
    ggsave(path %&% "." %&% major.type %&% "_final_celltype_scores.png", pp,
      width = 30, height = 13, units = "cm", dpi = 300)
  }
  write.table(t(final.score), path %&% "." %&% major.type %&% "_final_celltype_scores.tsv",
    sep = "\t", row.names = T, col.names = NA)
}


################################################################################
## main code starts here
################################################################################
## load input data
sce_data <- readRDS(opt$SCE)

################################################################################
## load celltype gene list -> all types in one file, subtype distinction via config file
if (endsWith(opt$celltype_lists, ".gmt")) {
  tmp <- readLines(opt$celltype_lists)
  tmp <- lapply(tmp, function(x) strsplit(x, "\\\t")[[1]])
  names(tmp) <- sapply(tmp, function(x) x[1])
  cell.type <- sapply(tmp, function(x) x[-1])
} else if (endsWith(opt$celltype_lists, ".gmx")) {
  cell.type <- read.table(opt$celltype_lists, sep = "\t", head = T, stringsAsFactors = F)
  cell.type <- apply(cell.type[-1, ], 2, function(x) x[x != ""])
} else {
  print("Error. File type of cell type list must be either .gmt or .gmx.")
  quit(status = 1)
}

all.ct.genes <- as.character(unlist(cell.type))
nr.genes <- sapply(cell.type, length)
cat("\n\nNumber of genes on each cell type list:\n\n")
print(nr.genes)
cat("\n\n")

# celltype config file, sub types are individual for each major type
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

# Keep list of all possible final cell types
all_celltypes_full_ct_name <- apply(type_config, 1, function(x) paste(x, collapse = ","))
all_celltypes_full_ct_name <- paste(all_celltypes_full_ct_name, collapse = ",")
all_celltypes_full_ct_name <- gsub("^NA,|,NA", "", all_celltypes_full_ct_name)
all_celltypes_full_ct_name <- strsplit(all_celltypes_full_ct_name, ",")[[1]]
idx <- which(all_celltypes_full_ct_name != "none")
all_celltypes_full_ct_name <- unique(all_celltypes_full_ct_name[idx])

if (hasName(metadata(sce_data), "all_celltypes_full_ct_name")) {
  metadata(sce_data)$all_celltypes_full_ct_name <- all_celltypes_full_ct_name
} else {
  metadata(sce_data) <- c(metadata(sce_data), list(all_celltypes_full_ct_name = all_celltypes_full_ct_name))
}
xx <- gsub(pattern = "([^_]+)_.*", "\\1", all_celltypes_full_ct_name)
if (hasName(metadata(sce_data), "all_celltypes")) {
  metadata(sce_data)$all_celltypes <- xx
} else {
  metadata(sce_data) <- c(metadata(sce_data), list(all_celltypes = xx))
}
# Keep list of all possible major cell types
if (hasName(metadata(sce_data), "all_major_types_full_ct_name")) {
  metadata(sce_data)$all_major_types_full_ct_name <- major_types
} else {
  metadata(sce_data) <- c(metadata(sce_data), list(all_major_types_full_ct_name = major_types))
}
xx <- gsub(pattern = "([^_]+)_.*", "\\1", major_types)
if (hasName(metadata(sce_data), "all_major_types")) {
  metadata(sce_data)$all_major_types <- xx
} else {
  metadata(sce_data) <- c(metadata(sce_data), list(all_major_types = xx))
}

# filter out non-unique genes as long as more than min_genes are left
genes_dupl <- which(duplicated(unlist(cell.type[major_types])))
genes_dupl <- as.character(unlist(cell.type[major_types])[genes_dupl])
# prelim new celltype-specific gene list
test_major_types <- lapply(cell.type[major_types], function(x) setdiff(x, genes_dupl))
# check if length > min_genes
tmp <- which(sapply(test_major_types, length) >= min_genes)
tmp <- setdiff(names(test_major_types), names(tmp))
if (length(tmp) > 0) {
  # final new celltype specific gene list
  # if there is a cell type in temp that has less than min_genes genes
  # replace it with the cell type list BEFORE filtering for unique genes
  for (ii in seq(length(tmp))) {
    test_major_types[[tmp[ii]]] <- as.character(cell.type[[tmp[[ii]]]])
  }
}


# perform the first celltyping
first.score <- f.score.ctgenes.U(sce_data, cell.type[major_types], min_genes)
first.class <- f.annot.ctgenes(first.score, thresh_unknown, thresh_uncert)
first.class$cell.type <- factor(first.class$cell.type, levels = c(major_types, "unknown", "uncertain"))
# attach major celltype to SCE object
colData(sce_data)$celltype_major_full_ct_name <- first.class$cell.type
colData(sce_data)$celltype_major <- as.factor(gsub(pattern = "([^_]+)_.*", "\\1", first.class$cell.type))
all.major.types <- sort(c(metadata(sce_data)$all_major_types, "uncertain", "unknown"))
cat("\n\nall.major.types in alphabetical order:\n\n")
print(all.major.types)
sce_data$celltype_major <- factor(sce_data$celltype_major, levels = all.major.types)

################################################################################
## plot major score heatmaps
bad <- apply(first.score, 2, function(x) length(unique(x)))
notbad <- which(bad != 1)
col.pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
hm1 <- pheatmap(-log10(first.score[, notbad]), scale = "none", color = col.pal, silent = T,
  show_rownames = T, show_colnames = F, clustering_method = "ward.D2",
  treeheight_row = 0, treeheight_col = 20)
hm2 <- pheatmap(-log10(first.score[, notbad]), scale = "column", color = col.pal, silent = T,
  show_rownames = T, show_colnames = F, clustering_method = "ward.D2",
  treeheight_row = 0, treeheight_col = 20)
pp <- plot_grid(hm1$gtable, hm2$gtable, nrow = 2)
ggsave(path %&% ".first_celltype_scores.png", pp,
  width = 30, height = 13, units = "cm", dpi = 300)
write.table(t(first.score[, notbad]), path %&% ".first_celltype_scores.tsv", sep = "\t", row.names = T, col.names = NA)

#########################################################
## Rphenograph of celltype-specific genes.
tmp <- assay(sce_data, "pearson_resid")
idx <- match(all.ct.genes, rownames(tmp))
idx <- idx[!is.na(idx)]
tmp <- tmp[idx, ]
r.umap <- umap(t(tmp), n_neighbors = n_nn, spread = 1, min_dist = 0.01,
  ret_nn = T)
reducedDim(sce_data, "umap_ct") <- r.umap$embedding
if (hasName(metadata(sce_data), "umap_ct")) {
  metadata(sce_data)$umap_ct <- r.umap$nn$euclidean
} else {
  metadata(sce_data) <- c(metadata(sce_data), list(umap_ct = r.umap$nn$euclidean))
}

#########################################################
## assign unknown or uncertain cells to the majority celltype of their nearest neighbors
prelim.celltype <- sce_data$celltype_major_full_ct_name
idx <- which(prelim.celltype %in% c("unknown", "uncertain"))
uu.nn <- r.umap$nn$euclidean$idx[idx, , drop = F]
# get preliminary cell types of the nearest neighbors of each uncertain/unknown cell
# disregard the cell itself
uu.nn.ct <- prelim.celltype[as.numeric(uu.nn[, -1])]
# have preliminary cell types shaped as matrix
# for each uncertain/unknown cell there is a row with (n_nn - 1) columns.
# the columns contain the cell types of the nearest neighbours of the cell.
uu.nn.ct <- matrix(uu.nn.ct, ncol = ncol(uu.nn) - 1, byrow = F)
# deal with definitely unknown/uncertain cases:
id_clear <- apply(uu.nn.ct, 1, function(x) all(x %in% c("unknown", "uncertain")))
id_unclear <- which(!id_clear)
# deal with rest:
f.vote <- function(x) names(sort(table(x[!x %in% c("unknown", "uncertain")]), decreasing = T)[1])
# f.vote(uu.nn.ct[1,])
vote.ct <- apply(uu.nn.ct[id_unclear, , drop = F], 1, f.vote)
# replace "unknown" and "uncertain" with the majority celltype of their nearest neighbors
prelim.celltype[idx[id_unclear]] <- vote.ct
cat("\nResults of first cell type classification:\n\n")
as.matrix(table(sce_data$celltype_major_full_ct_name))


# perform second celltyping
# prefill final ct column with results of major cell type classification
colData(sce_data)$celltype_final_full_ct_name <- as.character(sce_data$celltype_major_full_ct_name)
# iterate over all major cell types (ii is index of major cell type in list `minor_types`)
for (ii in seq_len(length(minor_types))) {
  # for each major cell type with minor defined, get all cells with this major cell type as prelim.celltype
  idy <- which(prelim.celltype %in% names(minor_types)[ii])
  cat("\n\nPerform second cell typing for", length(idy), "cells of major type:   ", names(minor_types)[ii], "\n")
  if (length(idy) > 0) {
    these.cells <- sce_data[, idy]
    these_major <- these.cells$celltype_major_full_ct_name
    these.ct <- match(minor_types[[ii]], names(cell.type))
    these.ct <- these.ct[!is.na(these.ct)]
    second.score <- f.score.ctgenes.U(these.cells, cell.type[these.ct], min_genes)
    second.class <- f.annot.ctgenes(second.score, thresh_unknown, thresh_uncert_second)
    table(second.class$cell.type)
    # cells with (second) label "unknown", "uncertain" keep the original first.class as final cell type
    idx <- which(second.class$cell.type %in% c("unknown", "uncertain"))
    if (length(idx) > 0) second.class$cell.type[idx] <- as.character(these_major[idx])
    sce_data$celltype_final_full_ct_name[idy] <- second.class$cell.type
    # for all successful second classifications, replace major celltype with parent of minor type
    idx <- which(!second.class$cell.type %in% c("unknown", "uncertain"))
    if (length(idx) > 0) sce_data$celltype_major_full_ct_name[idy][idx] <- names(minor_types)[ii]
    f_plot_final_score(second.score, names(minor_types)[ii])
  }
}

colData(sce_data)$celltype_final <- gsub(pattern = "([^_]+)_.*", "\\1", colData(sce_data)$celltype_final_full_ct_name)
all.cell.types <- sort(c(metadata(sce_data)$all_celltypes, "uncertain", "unknown"))
print("all.cell.types in alphabetical order:")
print(all.cell.types)
sce_data$celltype_final <- factor(sce_data$celltype_final, levels = all.cell.types)
cat("\n\nsce_data$celltype_final after second cell typing:\n\n")
as.matrix(table(sce_data$celltype_final))
cat("\nAdapted first cell type classification:\n\n")
as.matrix(table(sce_data$celltype_major_full_ct_name))
# have colData(sce_data)$celltype_final_full_ct_name as factors
all.cell.types.full <- sort(c(metadata(sce_data)$all_celltypes_full_ct_name, "uncertain", "unknown"))
# print("all.cell.types.full in alphabetical order:")
# print(all.cell.types.full)
colData(sce_data)$celltype_final_full_ct_name <- factor(colData(sce_data)$celltype_final_full_ct_name, levels = all.cell.types.full)


# keep sce object including cluster 0
sce_data_incl_cluster0 <- sce_data

# remove cells that are in cluster 0
mask_keep_cells <- colData(sce_data)$phenograph_clusters != 0
cat("\n\nNumber of cells that remain after filtering out cluster 0:\n\n")
print(sum(mask_keep_cells))
cat("\n\nsce object before filtering out cells that are in cluster 0:\n\n")
print(sce_data)

sce_data <- sce_data[, mask_keep_cells]
reducedDims(sce_data)$phenodist <- reducedDims(sce_data)$phenodist[, mask_keep_cells]
cat("\n\nsce object after filtering out cells that are in cluster 0:\n\n")
print(sce_data)

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
## begin new cluster purity calculations
# if at least one minor_type, make major-minor dictionary:
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
write.table(clu.pu, file = path %&% ".phenograph_celltype_association.txt", sep = "\t", row.names = F)
## end new cluster purity calculations


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

################################################################################
## save sce object for later use
saveRDS(sce_data, path %&% ".RDS")
# also save sce_data object that including the cells that are in cluster 0
saveRDS(sce_data_incl_cluster0, path %&% ".all_cells.RDS")
