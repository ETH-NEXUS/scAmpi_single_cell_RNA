#######################################
### Prepare session, load package
#######################################
#library(rDGIdb)
library(httr)
library(jsonlite)
library(R6)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

#######################################
# Functions queryDGIdb.R
#######################################
queryDGIdb <- function(genes,
                       sourceDatabases = NULL,
                       geneCategories = NULL,
                       interactionTypes = NULL) {

    if (missing(genes)) stop("Need to specify a vector of genes to query.")

    if (is.null(genes) || length(genes) == 0 || !is.character(genes)) {
        stop("Need to specify a non-empty vector of genes names.")
    }

    if (missing(sourceDatabases) || is.null(sourceDatabases)) {
        databases <- NULL
    } else {
        databases <- match.arg(arg = sourceDatabases,
                               choices = sourceDatabases(),
                               several.ok = TRUE)
        databases <- paste(databases, collapse = ",")
    }
    if (missing(geneCategories) || is.null(geneCategories)) {
        categories <- NULL
    } else {
        categories <- match.arg(arg = geneCategories,
                                choices = geneCategories(),
                                several.ok = TRUE)
        categories <- paste(categories, collapse=",")
    }
    if (missing(interactionTypes) || is.null(interactionTypes)) {
        interactions <- NULL
    } else {
        interactions <- match.arg(arg = interactionTypes,
                                  choices = interactionTypes(),
                                  several.ok = TRUE)
        interactions <- paste(interactions, collapse = ",")
    }

    # if (missing(curatedOnly)) {
    #     trustLevel <- NULL
    # } else if (!is.logical(curatedOnly)) {
    #     stop("Argument curatedOnly has to be logical (TRUE or FALSE)")
    # } else if (curatedOnly) {
    #     trustLevel <- "Expert%20cureated"
    # } else {
    #     trustLevel <- NULL
    # }

    # Check internet connection
    tryCatch({
        msg <- ""
        r <- GET("https://dgidb.org/api/v2/interaction_types.json")
        if (status_code(r) != 200) {
            msg <- "DGIdb service not available."
        }
    }, error = function(err) {
        msg <- "Check internet connection"
    })
    if (msg != "") stop(msg)

    # Query DGIdb
    cat("Querying DGIDB...")
    queryResult <- queryDgidbPost(genes,
                                  interactionSources = databases,
                                  geneCategories = categories,
                                  interactionTypes = interactions)
    #,trustLevel = trustLevel)
    cat("done!\n")

    # Init result class: rDGIdbResult
    result <- new(Class = "rDGIdbResult")

    # Set unmatched terms if any
    if (length(queryResult$unmatchedTerms) != 0) {
      result <- setUnmatchedTerms(result, queryResult$unmatchedTerms)
    }

    # Set matched terms and populate different formats of result tables
    if (!is.null(queryResult$matchedTerms$searchTerm)) {

        # Set result data
        result <- setData(result, queryResult$matchedTerms)

        # Populate result summary
        result <- setResultSummary(result)

        # Populate by gene table
        result <- setByGene(result)

        # Populate search term summary
        result <- setSearchTermSummary(result)

        #Populate detailed results
        if (nrow(result@resultSummary) > 0) result <- setDetailedResults(result)
    }

    return(result)
    # End of function queryDGIdb()
}

# Uses httr POST function to query DGIdb. Post instead of get allows
# long list of genes to be queried.
queryDgidbPost <- function(genes, interactionSources, geneCategories,
                           interactionTypes) {
    url <- "https://dgidb.org/api/v2/interactions.json"
    body <- list(genes = paste(genes, collapse = ","),
                 interaction_sources = paste(interactionSources, collapse = ","),
                 gene_categories = paste(geneCategories, collapse = ","),
                 interaction_types = paste(interactionTypes, collapse = ","))
    #source_trust_levels = trustLevel)
    body <- body[!sapply(body, is.null)]
    postRequest <- POST(url = url, body = body, encode = 'multipart')
    text <- content(postRequest, as = "text", encoding = "ISO-8859-1")
    if (grepl('error|DOCTYPE', text)) stop("Oops, badly formatted query.")
    if (identical(text, "")) stop("Query response was emtpy.")
    result <- fromJSON(text, simplifyVector = TRUE)
    return(result)
}

sourceDatabases <- function() {
    url <- "https://dgidb.org/api/v2/interaction_sources.json"
    result <- queryAPIget(url)
    return(result)
}

geneCategories <- function() {
    url <- "https://dgidb.org/api/v2/gene_categories.json"
    result <- queryAPIget(url)
    return(result)
}

interactionTypes <- function() {
    url <- "https://dgidb.org/api/v2/interaction_types.json"
    result <- queryAPIget(url)
    return(result)
}

queryAPIget <- function(url) {
    getRequest <- GET(url = url)
    text <- content(getRequest, as = "text", encoding = "ISO-8859-1")
    if (grepl('error|DOCTYPE', text)) stop("Oops, badly formatted query.")
    if (identical(text, "")) stop("Query response was emtpy.")
    result <- fromJSON(text, simplifyVector = TRUE)

    return(result)
}

###########################################
# Functions rDGIdbResult.R
###########################################
rDGIdbResult <- setClass(

    "rDGIdbResult",

    slots = c(
        matchedTerms = "character",
        unmatchedTerms = "character",
        resultSummary = "data.frame",
        detailedResults = "data.frame",
        byGene = "data.frame",
        searchTermSummary = "data.frame",
        data = "data.frame" # Matched terms raw data returned from DGIdb
    )
)

setGeneric(name = "setData",
           def = function(theObject, data)
               standardGeneric("setData")
)
setMethod(f = "setData",
          definition = function(theObject, data) {
              theObject@data <- data
              theObject@matchedTerms <- data$searchTerm
              return(theObject)
          })

setGeneric(name = "getData",
           def = function(theObject)
               standardGeneric("getData")
)
setMethod(f = "getData",
          definition = function(theObject) {
              return(theObject@data)
          })

setGeneric(name = "setResultSummary",
           def = function(theObject)
               standardGeneric("setResultSummary")
)
setMethod(f = "setResultSummary", definition = function(theObject) {
    # List of sources for which interactions were found with the genes given
    uniqueSources <- unique(unlist(sapply(
        theObject@data$interactions,
        function(x) {unique(x$source)}, simplify = FALSE)))

    # Summarize interactions similar to web interface
    interactionList <- lapply(theObject@data$geneName[
        sapply(theObject@data$interactions, length) != 0],
        getResultSummary, theObject@data, uniqueSources)
    interactionData <- data.frame(do.call(rbind, interactionList),
                                  stringsAsFactors = FALSE)
    if(nrow(interactionData) > 0) {
        # Convert DB and other columns to numeric
        interactionData[,3:(ncol(interactionData) - 2)] <-
            sapply(3:(ncol(interactionData) - 2),
                   function(x) as.numeric(interactionData[,x]))
        interactionData <- interactionData[
            order(as.numeric(interactionData$Score),
                  decreasing = TRUE),]
        rownames(interactionData) <- seq_len(nrow(interactionData))
    }
    theObject@resultSummary <- interactionData
    return(theObject)
})

setGeneric(name = "setByGene",
           def = function(theObject)
               standardGeneric("setByGene")
)
setMethod(f = "setByGene",
          signature = "rDGIdbResult",
          definition = function(theObject) {
              nrow <- length(theObject@matchedTerms)
              tmp <- data.frame(matrix(nrow = nrow, ncol = 5,
                                       dimnames = list(NULL, c('SearchTerm','Gene','GeneName',
                                                               'DistinctDrugCount','DruggableGeneCategories'))),
                                stringsAsFactors = FALSE)
              tmp[,c('SearchTerm','Gene','GeneName')] <-
                  theObject@data[,
                                 c('searchTerm','geneName','geneLongName')]
              tmp$DistinctDrugCount <-
                  sapply(theObject@data$interactions, function(x) {
                      length(unique(x$drugName)) } )
              tmp$DruggableGeneCategories <-
                  sapply(theObject@data$geneCategories, function(x) {
                      paste(x$name, collapse=",")
                  })
              tmp <- tmp[order(tmp$SearchTerm),]
              theObject@byGene <- tmp
              return(theObject)
          })

setGeneric(name = "setSearchTermSummary",
           def = function(theObject)
               standardGeneric("setSearchTermSummary")
)
setMethod(f = "setSearchTermSummary",
          signature = "rDGIdbResult",
          definition = function(theObject) {

              names <- c('SearchTerm','MatchType','Matches')
              tmp <- data.frame(
                  matrix(nrow = 0, ncol = 3, dimnames = list(NULL, names)),
                  stringsAsFactors = FALSE)
              if (length(theObject@matchedTerms) > 0) {
                  matchedTermSummary <- data.frame(cbind(theObject@matchedTerms,
                                                         rep('Definite', length(theObject@matchedTerms)),
                                                         theObject@data$geneName), stringsAsFactors = FALSE)
                  colnames(matchedTermSummary) <- names
                  tmp <- rbind(tmp, matchedTermSummary)
              }
              # if (length(theObject@unmatchedTerms) > 0) {
              #     noneId <- is.null(theObject@unmatchedTerms$suggestions) |
              #         sapply(theObject@unmatchedTerms$suggestions, length) == 0
              #     unmatchedTermsVector <- unlist(strsplit(
              #       theObject@unmatchedTerms$searchTerm[noneId], split = ', '))
              #     # Unmatched terms without suggestions
              #     unmatchedTermsNone <- cbind(unmatchedTermsVector,
              #         rep('None', length(unmatchedTermsVector)),
              #             rep('None', length(unmatchedTermsVector)))
              #     colnames(unmatchedTermsNone) <- names
              #     tmp <- rbind(tmp, unmatchedTermsNone)
              #     # Unmatched terms with suggestions
              #     unmatchedTermsSuggestions <-
              #         cbind(theObject@unmatchedTerms$searchTerm[!noneId],
              #             rep('Ambiguous', sum(!noneId)),
              #             sapply(theObject@unmatchedTerms$suggestions[!noneId],
              #                 function(x) { paste(x, collapse=',') } ))
              #     colnames(unmatchedTermsSuggestions) <- names
              #     tmp <- rbind(tmp, unmatchedTermsSuggestions)
              # }
              tmp <- tmp[order(tmp$SearchTerm),]
              theObject@searchTermSummary <- tmp
              return(theObject)
          })

setGeneric(name = "setDetailedResults",
           def = function(theObject)
               standardGeneric("setDetailedResults")
)
setMethod(f = "setDetailedResults",
          signature = "rDGIdbResult",
          definition = function(theObject) {
              tmp <- do.call(rbind, apply(theObject@data, 1, function(x) {
                  nrow <- nrow(x$interactions)
                  if (nrow > 0) {
                      data.frame(cbind(rep(x$searchTerm, times = nrow),
                                       rep(x$geneName, times = nrow),
                                       x$interactions[,
                                                      c('drugName', 'interactionTypes',
                                                        'sources','pmids')]),
                                 stringsAsFactors = FALSE)
                  }
              }))
              tmp$interactionTypes <-
                  sapply(tmp$interactionTypes, paste, collapse = ",")
              tmp$sources <- sapply(tmp$sources, paste, collapse = ",")
              tmp$pmids <- sapply(tmp$pmids, paste, collapse = ",")
              colnames(tmp) <-
                  c('SearchTerm', 'Gene', 'Drug',
                    'InteractionType', 'Source','PMIDs')
              tmp <- tmp[order(tmp$SearchTerm),]
              rownames(tmp) <- seq_len(nrow(tmp))
              theObject@detailedResults <- tmp
              return(theObject)
          })

setGeneric(name = "setUnmatchedTerms",
           def = function(theObject, unmatchedTerms)
               standardGeneric("setUnmatchedTerms")
)
setMethod(f = "setUnmatchedTerms",
          signature = "rDGIdbResult",
          definition = function(theObject, unmatchedTerms) {
              theObject@unmatchedTerms <- unmatchedTerms
              # data.frame(cbind(unmatchedTerms,
              #                  rep(NULL, length(unmatchedTerms))),
              #            stringsAsFactors = FALSE)
              return(theObject)
          })
setGeneric(name = "resultSummary",
           def = function(object)
               standardGeneric("resultSummary")
)
setMethod(f = "resultSummary",
          signature = "rDGIdbResult",
          definition = function(object) {
              return(object@resultSummary)
          })
setGeneric(name = "detailedResults",
           def = function(object)
               standardGeneric("detailedResults")
)
setMethod(f = "detailedResults",
          signature = "rDGIdbResult",
          definition = function(object) {
              return(object@detailedResults)
          })
setGeneric(name = "byGene",
           def = function(object)
               standardGeneric("byGene")
)
setMethod(f = "byGene",
          signature = "rDGIdbResult",
          definition = function(object) {
              return(object@byGene)
          })
setGeneric(name = "searchTermSummary",
           def = function(object)
               standardGeneric("searchTermSummary")
)
setMethod(f = "searchTermSummary",
          signature = "rDGIdbResult",
          definition = function(object) {
              return(object@searchTermSummary)
          })


getResultSummary <- function(gene, output, sources) {
    # Row index with gene interaction information
    idx <- which(output$geneName == gene)
    # Prepare table of interactions (drug vs. interaction DB)
    foundSources <- unique(unlist(output[idx,]$interactions[[1]]$sources))
    result <- #data.frame(cbind(output[idx,]$interactions[[1]]$drugName,
        as.data.frame(t(sapply(output[idx,]$interactions[[1]]$sources,
                               function(x, y) { return(y %in% x) }, foundSources)),
                      stringsAsFactors = FALSE)
    if (length(foundSources) == 1) {result <- t(result)} # result is vector
    # dimnames(result) <-
    #     list(output[idx,]$interactions[[1]]$drugName, foundSources)
    # Expand matrix to all possible DBs, set multiple occurances of
    # interactions to one, and add gene and drug names
    tmp <- data.frame(matrix(0, nrow = nrow(result), ncol = 4 + length(sources),
                             dimnames = list(NULL,c('Gene', 'Drug', sources, 'PMID','Score'))),
                      stringsAsFactors = FALSE, check.names = F)
    #tmp[, foundSources] <- 1*result
    if(length(foundSources) > 0) { # no sources or source="NULL"
      tmp[, foundSources] <- 1*result
    }
    tmp[tmp > 1] <- 1 # Remove double counts
    tmp$Score <- rowSums(tmp) +
        sapply(output[idx,]$interactions[[1]]$pmids, length)
    tmp$Gene <- rep(gene, nrow(tmp))
    tmp$Drug <- output[idx,]$interactions[[1]]$drugName
    tmp$PMID <- sapply(output[idx,]$interactions[[1]]$pmids, paste, collapse=",")
    # Determine type of interaction
    #resultType <- table(output[idx,]$interactions[[1]]$drugName,
    #                    output[idx,]$interactions[[1]]$interactionType)
    #if (nrow(tmp) == 1) {
    #    tmp$Type <- paste(colnames(resultType), collapse = ",")
    #} else {
    #    listResult <- lapply(split(resultType, seq(nrow(resultType))),
    #                         function(x, names) { names[x>0] },
    #                         colnames(resultType))
    #    tmp$Type <- sapply(listResult, paste, collapse=',')
    #}
    return(as.matrix(tmp))
}



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


