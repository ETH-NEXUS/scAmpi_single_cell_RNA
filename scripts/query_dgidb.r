#######################################
### Prepare session, load package
#######################################
library(rDGIdb)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

### Check input arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript Query_DGIDB.R path/to/input path/to/output #DBs")
} else {
    srcFile <- args[1]
    outFile <- args[2]
    MIN_DB_SUPPORT <- as.integer(args[3])
    colName_genes <- args[4]
}

### Read input file
if (!file.exists(srcFile)) {
    stop("Input file does not exist!")
} else {
    if (file.info(srcFile)$size == 0) {
		# write empty files, necessary for snakemake pipeline
		interactionData_file=paste(outFile, 'CompleteTable', 'txt', sep='.')
		categories_file=paste(outFile, 'GeneCategories', 'txt', sep='.')
		interactionData_cols = c('Gene','Drug','PMID','Score','Type')
		write(interactionData_cols, file=interactionData_file, ncolumns=length(interactionData_cols), sep='\t')

		categories_cols = c("Gene","GeneName","Category")
		write(categories_cols, file=categories_file, ncolumns=length(categories_cols), sep='\t')

		output_cols = c('Gene', 'Drug')
		write(output_cols, file=outFile, ncolumns=length(output_cols), sep='\t')
        quit(save = "default", status = 0)
    } else {
        input <- read.table(srcFile, sep='\t', header = TRUE, stringsAsFactors = FALSE)
    }
}

### Map to unique official gene symbols (limma package)
#if (ncol(input) != 1) stop("Wrong input format, single column gene list expected!")
genes <- unique(input[[colName_genes]])
geneCategories2 <- c("CLINICALLY ACTIONABLE", "DRUGGABLE GENOME", "DRUG RESISTANCE")

result = NULL  # initialize
emptyQuery = FALSE
if (length(genes) == 0) {
	emptyQuery = TRUE
} else {
	result <- queryDGIdb(genes, geneCategories = geneCategories2)
	if(length(result@matchedTerms) == 0) emptyQuery = TRUE
	if (length(resultSummary(result)) == 0) emptyQuery = TRUE
}

if(emptyQuery == TRUE) {	
	# write empty files, necessary for snakemake pipeline
	interactionData_file=paste(outFile, 'CompleteTable', 'txt', sep='.')
	categories_file=paste(outFile, 'GeneCategories', 'txt', sep='.')
	interactionData_cols = c('Gene','Drug','PMID','Score','Type')
    write(interactionData_cols, file=interactionData_file, ncolumns=length(interactionData_cols), sep='\t')

    categories_cols = c("Gene","GeneName","Category")
    write(categories_cols, file=categories_file, ncolumns=length(categories_cols), sep='\t')

    output_cols = c('Gene', 'Drug')
    write(output_cols, file=outFile, ncolumns=length(output_cols), sep='\t')
	
	quit(save = "default", status = 0)
}

#######################################
### Prepare and write results
######################################

### Complete interaction table
interactionData <- resultSummary(result)
interactionData$Type <- apply(interactionData, 1, function(x, details) {
    idx <- which(x['Gene'] == details$Gene & x['Drug'] == details$Drug)
    if (nchar(details$InteractionType[idx]) == 0){'.'}
    else{
    paste(sort(unique(details$InteractionType[idx])), collapse = ",")
    }
}, result@detailedResults)
interactionData$PMID <- apply(interactionData, 1, function(x, details) {
    idx <- which(x['Gene'] == details$Gene & x['Drug'] == details$Drug)
    if (nchar(details$PMIDs[idx]) == 0){'.'}
    else{
        paste(sort(unique(details$PMIDs[idx])), collapse = ",")
    }
}, result@detailedResults)
#interactionData <- interactionData[interactionData$Score >= MIN_DB_SUPPORT,]
interactionData <- interactionData[as.integer(interactionData$Score) >= MIN_DB_SUPPORT,]

# if results are found but never with enough support, interactionData will be empty after MIN_DB_SUPPORT filtering and script crashes. Avoid by writing dummy output files:
emptyQuery = FALSE
if(length(interactionData$Drug) == 0) emptyQuery = TRUE

if(emptyQuery == TRUE) {	
	# write empty files, necessary for snakemake pipelinie
    print("Create empty output file (search results are empty)")
	interactionData_file=paste(outFile, 'CompleteTable', 'txt', sep='.')
	categories_file=paste(outFile, 'GeneCategories', 'txt', sep='.')
	interactionData_cols = c('Gene','Drug','PMID','Score','Type')
    write(interactionData_cols, file=interactionData_file, 
          ncolumns=length(interactionData_cols), sep='\t')

    categories_cols = c("Gene","GeneName","Category")
    write(categories_cols, file=categories_file, 
          ncolumns=length(categories_cols), sep='\t')

    output_cols = c('Gene', 'Drug')
    write(output_cols, file=outFile, 
          ncolumns=length(output_cols), sep='\t')
	
	quit(save = "default", status = 0)
}
write.table(x = interactionData,
            file = paste(outFile, 'CompleteTable', 'txt', sep='.'),
            sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

### Gene categories file
categories <- byGene(result)
categories <- categories[,c('Gene','GeneName','DruggableGeneCategories')]
# Remove genes that dropped out from MIN_DB_SUPPORT test
categories <- categories[categories$Gene %in% interactionData$Gene,]
write.table(x = categories,
            file = paste(outFile, 'GeneCategories', 'txt', sep='.'),
            sep="\t", row.names = FALSE, 
            col.names = c("Gene","GeneName","Category"), quote = FALSE)


### Write summary interaction table
collapseInteractionTable <- function(gene, interactionData) {
    geneDrug <- c(gene, paste(interactionData$Drug[interactionData$Gene == gene], collapse = ","))
}
interactionData$Drug <- paste(interactionData$Drug, " ", "(", interactionData$Score, ")", sep = "")
geneDrugTable <- t(sapply(unique(interactionData$Gene), 
                          collapseInteractionTable, 
                          interactionData))
colnames(geneDrugTable) <- c("Gene", "Drug")

write.table(x = geneDrugTable, file = outFile, sep="\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)
