########################################

####    File name: celltype_summary.R
####    Author: Irene Keller
####    April 2024
####    R Version: 4.3

########################################

# This script takes as input the final cell type assignments for each cell
# and outputs a summary table with number of cells per cell type and sample

# R libraries
library("optparse")
library("dplyr")

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")


# parse command line arguments
option_list <- list(
  make_option("--input", type = "character", help = "List of all input files", multiple = TRUE),
  make_option("--output", type = "character", help = "Output filename")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


infiles <- opt$input
print("these are my input files")
print(infiles)

cells_per_sample <- lapply(infiles, 
                            read.table,
                            header = TRUE,
                            sep = "\t")

names(cells_per_sample) <- lapply(cells_per_sample, function(x) {
  filename <- basename(x)
  filename <- sub("\\.cts_final.txt$", "", filename)
  return(filename)
})

print(names(cells_per_sample))

# Summarize the counts
count_cell_types <- function(df) {
  count(df, celltype_final)
}
counts_per_sample <- lapply(cells_per_sample, count_cell_types)
names(counts_per_sample) <- names(cells_per_sample)

# Add the sample name to each table
for(sample in names(counts_per_sample)){
  counts_per_sample[[sample]] <- counts_per_sample[[sample]] %>%
    mutate(sample = sample)
}

# merge all summary tables and convert to wide format
summary_table_long <- counts_per_sample %>%
  do.call("rbind", .)
summary_table_wide <- reshape(summary_table_long, idvar = "celltype_final", timevar = "sample", direction = "wide")
write.table(summary_table, opt$output, row.names = FALSE, quote = FALSE, sep = "\t")


