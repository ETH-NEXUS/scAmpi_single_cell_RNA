################################################################################
## Merge sce object from scRNA preprocessing with the phenograph dist. matrix, mod.-score and clusters
## File name: prepare_celltyping.R
## Authors: Michael Prummer, Anne Bertolini, Franziska Singer
## Date created: June 2019
## R Version: 4.0
################################################################################

# This script merges the results of phenograph (the cluster IDs of the cells, the modularity score of the clustering and the distance matrix) into a given sce object.


lby = c("optparse", "rhdf5", "reshape2", "scran", "plyr")
resp = lapply(lby, function (x) suppressWarnings(suppressMessages(require(x, character.only=T,
                                                 warn.conflicts=F, quietly=T))))
if(!all(unlist(resp))) stop("Could not load one or more packages")
rm(resp, lby)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

#options(error=function()traceback(2))

# convenience function for string concatenation
'%&%' = function(a,b) paste(a,b,sep="")


option_list = list(
  make_option("--in_sce", type = "character", help = "Path to RDS file with sce object stored inside."),
  make_option("--phenograph_cluster", type = "character", help = "Path to the file clusters_phenograph.csv."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--distanceMatrix", type = "character", help = "Distance matrix used by phenograph."),
  make_option("--modularity_score", type = "character", help = "txt file containing modularity score of the phenograph clustering.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

path = opt$outputDirec %&% opt$sampleName

# read in sce object
my_sce = readRDS(opt$in_sce)

# read in phenograph cluster table
phcl = read.table(opt$phenograph_cluster, sep=",", head=T, stringsAsFactors = F)
names(phcl) = c("barcodes", "phenoclust")
print("Read in phenograph cluster table: ")
print(str(phcl))
print(head(phcl))

# integrate phenograph ID into sce object
# make sure the phenograph clusters IDs are in the right order, same as the cells in the sce object
stopifnot(length(which(colData(my_sce)$barcodes %in% phcl$barcodes)) == length(colData(my_sce)$barcodes))
cells_ordered <- as.data.frame(as.character(colData(my_sce)$barcodes))
names(cells_ordered) <- "barcodes"
cells_ordered <- join(cells_ordered, phcl, by = "barcodes")
colData(my_sce)$phenograph_clusters <- cells_ordered$phenoclust

# read in phenograph distance matrix
phenodist = read.table(opt$distanceMatrix, sep="\t", head=F)
phenodist = as.matrix(phenodist)
rownames(phenodist) = phcl$barcodes
colnames(phenodist) = phcl$barcodes
print("Read in phenograph distance matrix: ")
print(str(phenodist))
stopifnot(rownames(phenodist) == colData(my_sce)$barcodes)
all.equal(rownames(phenodist), colData(my_sce)$barcodes)

# read in modularity score of phenograph clustering
modularity_score = readLines(opt$modularity_score)
modularity_score <- as.numeric(modularity_score)
modularity_score <- signif(modularity_score, digits = 3)
print("Read in modularity_score of phenograph clustering: ")
print(modularity_score)

# Integrate modularity score into the sce object
metadata(my_sce)$modularity_score <- modularity_score
print("metadata(my_sce)$modularity_score:")
metadata(my_sce)$modularity_score

# Integrate the phenograph distance matrix into the sce object
reducedDim(my_sce, "phenodist") <- phenodist
print("reducedDim(my_sce, 'phenodist')[1:5,1:5]")
reducedDim(my_sce, "phenodist")[1:5,1:5]

# change the rownames of the sce object from Ensembl IDs to HGNC gene names
rownames(my_sce) <- rowData(my_sce)$SYMBOL
print("my_sce:")
print(my_sce)

# Save sce object
print("Save sce object in the following path: ")
print(paste0(opt$outputDirec, opt$sampleName, ".prep_celltyping.RDS"))
saveRDS(my_sce, opt$outputDirec %&% opt$sampleName %&% ".prep_celltyping.RDS")
