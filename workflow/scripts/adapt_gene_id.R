# Load required packages
library(Matrix)
library(EnsDb.Mmusculus.v79)
library(optparse)

# Define command line options
option_list <- list(
  make_option(c("-m", "--matrix"), type = "character", default = NULL, 
              help = "Path to the input feature barcode matrix file (.mtx)", metavar = "MATRIX"),
  make_option(c("-f", "--features"), type = "character", default = NULL, 
              help = "Path to the input feature file", metavar = "FEATURES")

)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if file option is provided
if (is.null(opt$matrix) || is.null(opt$features)) {
  stop("Please provide the paths to the input feature barcode matrix file (-m), the input barcode matrix file (-b), and the input feature file (-f).")
}

# Read in matrix
my_matrix <- readMM(opt$matrix)

# Read in the feature file
my_features <- read.table(opt$features, sep="\t", header = FALSE)

# Extract the SYMBOLS from the first column
symbol_ids <- my_features[, 1]
print(head(symbol_ids))

# Map SYMBOL IDs to Ensembl IDs
ensembl_ids <- mapIds(EnsDb.Mmusculus.v79, keys = symbol_ids, column = "GENEID", keytype = "SYMBOL", multiVals = "first")
print(ensembl_ids)

# Replace IDs with Ensembl IDs in the first column
my_features[, 1] <- ensembl_ids

# Write a list of all genes where no Ensembl ID could be assigned
na_features <- my_features[is.na(my_features[,1]),2, drop = FALSE]
colnames(na_features) <- "Genes_without_Ensembl_ID"
no_match_file <- gsub("matrix\\.original\\.mtx", "unmatched_genes\\.txt", opt$matrix)
write.table(na_features, no_match_file, row.names = FALSE, quote = FALSE)

# Filter out rows in the feature barcode matrix where the feature is NA
non_na_features <- my_features[!is.na(my_features[, 1]),1]
my_non_na_features <- my_features[!is.na(my_features[, 1]),]
print(head(non_na_features))
print("length non_na_features")
print(length(non_na_features))

print("length my_non_na_features")
print(length(my_non_na_features[,1]))
print("dimensions before removing nas")
print(dim(my_matrix))
print(str(my_matrix))

my_matrix <- my_matrix[!is.na(my_features[, 1]), , drop = FALSE]
print("dimensions after removing nas")
print(dim(my_matrix))
print(str(my_matrix))

print("Duplicated Entries:")
is_duplicated_features <- duplicated(my_non_na_features[,1]) | duplicated(my_non_na_features[,1], fromLast = TRUE)
#print(duplicated_features )
print(dim(my_matrix[!is_duplicated_features, , drop = FALSE]))


# Step 3: Sum up the rows in mtx file corresponding to duplicated features
duplicated_features <- my_non_na_features[duplicated(my_non_na_features[,1]) | duplicated(my_non_na_features[,1], fromLast = TRUE),]
print(duplicated_features)

for (feature_name in unique(duplicated_features[,1])) {
#  print(feature_name)
  duplicated_rows <- my_non_na_features[my_non_na_features[,1] == feature_name, ]
  print(duplicated_rows)
  print(as.numeric(rownames(duplicated_rows))[-1])
  sum_row <- colSums(my_matrix[as.numeric(rownames(duplicated_rows)), , drop = FALSE])
  my_matrix[as.numeric(rownames(duplicated_rows))[1], ] <- sum_row
  print(dim(my_matrix))
  my_matrix <- my_matrix[-as.numeric(rownames(duplicated_rows))[-1], ]
 

}
print(dim(my_matrix))
writeMM(my_matrix, file = gsub("\\.original", "", opt$matrix))


# Write the modified feature table back to a file
modified_file <- gsub("\\.original", "", opt$features)
my_non_na_features <- my_non_na_features[!duplicated(my_non_na_features[,1], fromLast = TRUE),]
print(length(my_non_na_features[,1]))
print(head(my_non_na_features))
write.table(my_non_na_features, file = modified_file, sep = "\t", row.names = FALSE, quote = FALSE, col.names=FALSE)


