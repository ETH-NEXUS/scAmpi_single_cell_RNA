#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")

# assumes object “res” is a dataframe containing the DEA results
# res must contain the following columns
#   “test_statistic” – statistics for ranking DEGs.
#   “gene_names” – gene symbols for DE genes


library(limma)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript GeneSetEnrichmentAnalysis.R path/to/inputTable path/to/outputFile path/to/geneSetDB.gmt columnName_testStatistic")
} else {
  inputTab <- args[1]
  outFile <- args[2]
  geneSetDB <- args[3]
  nameStat <- args[4]
}

dfTab = read.table(inputTab, header = TRUE, sep = '\t')
# dfTab = data.frame(myTab)
# dfTab = na.omit(dfTab) # MAST compatibility

if(endsWith(geneSetDB,'.gmt')){
  tmp = readLines(geneSetDB)
  tmp = lapply(tmp, function(x) strsplit(x, "\\\t")[[1]])
  names(tmp) = sapply(tmp, function(x) x[1])
  gset_temp = sapply(tmp, function(x) x[-1])
  gset = sapply(gset_temp, function(x) x[-1])
}else{
  print('Error. File type of cell type list must be .gmt.')
  quit(status = 1)
}

idx = ids2indices(gene.sets = gset, identifiers = dfTab$gene_names)
nr_genes_per_set = sapply(idx, length)
idx = idx[which(nr_genes_per_set > 4)]
if(length(idx) > 0){
  dat = cameraPR(dfTab[, nameStat], idx, sort=F)
  # non-directional
  dat$PValue.Mixed = cameraPR(abs(dfTab[[nameStat]]), idx, sort = F)$PValue
  dat$FDR.Mixed = p.adjust(dat$PValue.Mixed, method="BH")
  dat$name = rownames(dat)
  
  # define the direction of the gene set effect ("Up", "Down", "Mixed")
  dat$Direction = as.character(dat$Direction)
  dat$Direction[dat$FDR > 0.05] = "Mixed"
  dat$Direction[dat$FDR > 0.05 & dat$FDR.Mixed > 0.05] = "NOT"
  
  idx = which(dat$Direction=="Mixed")
  if(length(idx)>0) {
    dat$FDR[idx] = dat$FDR.Mixed[idx]
    dat$PValue[idx] = dat$PValue.Mixed[idx]
  }
  #if(length(idx)>0) dat$FDR[idx] = dat$FDR.Mixed[idx]
  dat = dat[,-grep("\\.Mixed", names(dat))]
  dat = dat[dat$Direction != "NOT", ]
  dat$Direction = factor(dat$Direction, levels=c("Up", "Down", "Mixed"))
  
  write.table(dat,file = outFile, sep="\t", row.names = FALSE, quote = FALSE)
}
