# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

args <- commandArgs(trailingOnly = TRUE)

inputFile <- args[1]
mapFile <- args[2]
outFile <- args[3]

cat("\n")
print(paste0("Processing file ", inputFile, "..."))
input <- read.table(inputFile, sep="\t", header=F, stringsAsFactors=F)
print(head(input,1))
print(paste0("Dimension: ",nrow(input)))

cat("\n")
print(paste0("Using as mapper file ", mapFile, "..."))
map <- read.table(mapFile, sep="\t", header=F, stringsAsFactors=F)
print(head(map,1))
print(paste0("Dimension: ",nrow(map)))

cat("\n")
print(paste0("Identical dim?: ",identical(nrow(map),nrow(input))))

## Tables have been loaded as single-column dataframes
input_ids <- input[,1]
map_ids <- map[,1]
interGenes <- intersect(input_ids,map_ids)
print(paste0("No. intersecting gene ids: ", length(interGenes)))
if (length(interGenes) < nrow(input)){
    cat("\n")
    stop("Error! Not all EnsemblIDs are available in provided mapper file.")
}

cat("\n")
print("Mapping gene symbols...")
mapper <- merge(data.frame(Ids=input_ids),data.frame(Ids=map_ids, Symbols=map[,2], Type=map[,3]), by="Ids", sort=F)
print(head(mapper,1))
mapGene <- intersect(mapper$Ids,map_ids)
print(paste0("Mapped gene ids: ",length(mapGene)))
if (length(mapGene) < nrow(input)){
    cat("\n")
    stop("Error! Something went wrong with the mapping.")
}
# splitFile <- strsplit(inputFile, "\\/")[[1]]
# splitFile <- paste(splitFile[1:(length(splitFile)-1)], collapse="/")
# outFile <- paste0(splitFile,"/features.tsv")
cat("\n")
print(paste0("Saving to new file ", outFile, "..."))
write.table(mapper, file=outFile, sep="\t", col.names=F, row.names=F, quote=F)
cat("\n")
