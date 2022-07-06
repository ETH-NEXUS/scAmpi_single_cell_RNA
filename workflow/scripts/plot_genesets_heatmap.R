library(reshape2)
library(ggplot2)
library(gridExtra)
library(gtable)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

###################
### FUNCTIONS
###################

## Extract and return the legend from a given ggplot
##  'https://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram'
g_legend<-function(a.gplot){
    ## Work around bug that opens a new plot device
    pdf(file=NULL)
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    dev.off()
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}


###################
### SCRIPT
###################

### Check input arguments

args <- commandArgs(trailingOnly = TRUE)

#if (length(args) < 2) {
#  # stop("Usage: Rscript plot_genesets_heatmap.R path/to/directory inFileTag path/to/outFile")
#  stop("Usage: Rscript plot_genesets_heatmap.R path/to/directory path/to/outFile")
#}

# assumes first parameter = output file
# second to nth parameter = input files

outfile = args[1]

## List all input files

allFiles = NULL
for (file_temp in 2:length(args)){
    allFiles = c(allFiles,args[file_temp])
}

#allFiles <- list.files(path=args[1], pattern="*.\\d+.enrichedGeneSets.txt",full.names=T)
# allFiles <- list.files(path=args[1], pattern=paste("*.\\d+.",args[2],sep=""),full.names=T)
cat("\nInput files:\n")
cat(paste(" ",allFiles,collapse="\n"), "\n")

## Load all input files into a list
data.files <- sapply(allFiles, function(f) read.table(f, sep="\t", header=T, quote="", check.names=F, fill=T, stringsAsFactors=F), simplify=F, USE.NAMES=T)

## Generate appropriate names that allow keeping track of each cluster ID
clusterIDs <- names(data.files)
clusterIDs <- lapply(clusterIDs, function(x){ 
                nameSplit <- strsplit(strsplit(x,split=".enrichedGeneSets.txt")[[1]],split="\\.")[[1]]
                # nameSplit <- strsplit(strsplit(x,split=paste(".",args[2],sep=""))[[1]],split="\\.")[[1]]
                nameSplit <- nameSplit[length(nameSplit)]
              })
clusterIDs <- paste("Cluster",clusterIDs,sep="_")
names(data.files) <- clusterIDs


## Retrieve all unique geneset names in the overall dataset
## This information is contained in the column "name"
genesetNames <- as.character(unique(unlist(sapply(data.files, function(x) subset(x,select=name)[,1], USE.NAMES=F))))

## Sanity check for no cluster having anything significant
## Save empty plot
if (length(genesetNames)<1){
  p_empty <- ggplot() + theme_void()
  pdf(file=NULL)
  ## Arrange main plot and 3 legends
  ggsave(outfile, plot=p_empty, width=20, height=20)
  garbage <- dev.off()
  cat("Saving empty plot...")

## When at least 1 cluster and 1 pathway were significant, continue
} else{
  ## Create 2 matrices gathering all data contained in the input files (FDR and Direction separately)
  gs_FDR <- matrix(data=NA, ncol=length(clusterIDs), nrow=length(genesetNames))
  colnames(gs_FDR) <- clusterIDs
  rownames(gs_FDR) <- genesetNames
  gs_Direction <- gs_FDR
  
  ## For each cluster, make sure that ordering of genesets is correct. Introduce NA when a geneset is not present
  for (i in seq_along(colnames(gs_FDR))){
    x <- colnames(gs_FDR)[i]
    gsCurrent <- subset(data.files[[x]], select=c(name,Direction,FDR))
    dummyTarget <- data.frame(name=rownames(gs_FDR), Indx=seq_along(rownames(gs_FDR)), stringsAsFactors=F, row.names=NULL)
    gsTarget <- merge(dummyTarget, gsCurrent, by="name", all=T)
    gsTarget <- gsTarget[order(gsTarget$Indx),,drop=F]
    gs_FDR[,i] <- as.numeric(gsTarget$FDR)
    gs_Direction[,i] <- gsTarget$Direction
  }
  
  ## Melt both dataframes to comply with ggplot input format
  m.gs_FDR <- melt(gs_FDR, na.rm=FALSE)
  m.gs_Direction <- melt(gs_Direction, na.rm=FALSE)
  colnames(m.gs_FDR) <- c("Genesets","Clusters","FDR")
  colnames(m.gs_Direction) <- c("Genesets","Clusters","Direction")
  
  ## Merge both molten dataframes to join columns FDR and Direction
  m.gs <- merge(m.gs_FDR, m.gs_Direction, by=c("Genesets","Clusters"), all=T)
  
  ## Remove prefix tag 'HALLMARK' from geneset names
  m.gs$Genesets <- gsub("HALLMARK_", "", m.gs$Genesets)
  
  ## Bin the FDR values to generate discrete categories (ie. FDR intervals)
  m.gs$FDR_bin <- cut(m.gs$FDR, breaks = c(0.05,0.035,0.02,0.005,0), include.lowest = TRUE)
  
#   ## Remove all FDRs outside the range of interest (ie. [0,0.05]), if any
#   ## These will have been assigned <NA> for variable FDR_bin after using cut() to bin the data
#   m.gs <- m.gs[!is.na(m.gs$FDR_bin),,drop=F]
  
  
  ### Cluster FDRs within heatmap
  
  ## Perform this step only when nclusters>1 and npaths>1 (ie. when there is an actual matrix available)
  nClusters <- unique(m.gs$Clusters)
  nPaths <- unique(m.gs$Genesets)
  if (length(nClusters) > 1 & length(nPaths) > 1){
    
    ## Retrieve subset of the data to run clustering algorithm on
    ## Non-existing combinations of cluster+pathway (due to lack of significance) will be assigned NA
    dend.df <- dcast(subset(m.gs,select=c(Genesets,Clusters,FDR)), Genesets ~ Clusters)
    rownames(dend.df) <- dend.df$Genesets
    dend.df <- subset(dend.df, select=-Genesets)
    ## Assign FDR=-1 to all NA (ie. non-significant) combinations of cluster+pathway
    dend.df[is.na(dend.df)] <- -1
    dend.df[] <- sapply(dend.df, as.numeric)
    dend.matrix <- as.matrix(dend.df)
    
    ## Run clustering algorithm over columns and rows, separately
    distance.row = dist(dend.matrix, method = "euclidean")
    cluster.row = hclust(distance.row, method = "average")
    distance.col = dist(t(dend.matrix), method = "euclidean")
    cluster.col = hclust(distance.col, method = "average")
    
    ## Extract the order of the tips in the corresponding dendrograms
    ord.row <- order.dendrogram(as.dendrogram(cluster.row))
    ord.col <- order.dendrogram(as.dendrogram(cluster.col))
    ## Order the levels of Genesets (rows) and Clusters (cols) according to their corresponding clustering
    m.gs$Genesets <- factor(x = m.gs$Genesets, levels = unique(m.gs$Genesets)[ord.row], ordered = TRUE)
    m.gs$Clusters <- factor(x = m.gs$Clusters, levels = unique(m.gs$Clusters)[ord.col], ordered = TRUE)
  }
  
  
  ## Assign color scales to each pathway category
  cols_down <- colorspace::diverge_hcl(12, h=c(246,40), c=96)[c(1:3,5)]
  cols_up <- c("red3","red","lightcoral","#FFDEDE")
  cols_mix <- c("#E89331","#EAB12A","#FEE17F","#F9E9BD")
  
  ## If all 3 scales are discrete, then we can generate 4x3 different labels (4 intensity values x 3 color scales) to make coloring easier
  m.gs$colorTag <- paste(m.gs$FDR_bin, m.gs$Direction, sep="_")
  interNames <- c("[0,0.005]","(0.005,0.02]", "(0.02,0.035]", "(0.035,0.05]")
  all_colors <- setNames(c(cols_up,cols_down,cols_mix,"white"),c(apply(expand.grid(interNames, c("Up","Down","Mixed")), 1, paste, collapse="_"),"NA_NA"))
  
  
  ### Plot the heatmap
  
  ## Plot main heatmap (no legends)
  
  p <- ggplot(data=m.gs, aes(Clusters, Genesets)) +
    geom_tile(aes(fill=colorTag)) +
    scale_fill_manual(values=all_colors, drop=F) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, vjust=1, size=20, hjust=1, color="black", face="bold"),
          axis.text.y = element_text(vjust=0.6, size=18, hjust=1, color="black", face="bold")) +
    coord_fixed(ratio=0.5) +
    theme(legend.title = element_text(colour="black", size=16, face="bold")) +
    theme(legend.text = element_text(colour="black", size = 16, face = "bold")) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.background = element_blank(), axis.ticks = element_blank()) +
    guides(fill=FALSE,alpha=FALSE)
  
  
  ## Generate dummy plots to retrieve separate legends for each color scale
  
  ## Dummy dataframe containing all the directions and possible FDR categories
  ## Avoids errors due to lack of all the categories and directions
  dummy_df <- data.frame(Genesets=NA,Clusters=NA,Direction=NA,FDR_bin=NA)
  for (dir in c("Down","Mixed","Up")){
    for (fdr_categ in interNames){
      dummy_df <- rbind(dummy_df, c("DummyGeneset","DummyCluster",dir,fdr_categ))
    }
  }
  ## Order the FDR categories appropriately
  dummy_df$FDR_bin <- factor(dummy_df$FDR_bin, levels=interNames)
  
  g_dummy1 <- ggplot(subset(dummy_df, Direction=="Down"), aes(Clusters, Genesets,fill=FDR_bin)) + geom_tile() +
    scale_fill_manual(name="FDR Down", values=setNames(cols_down, interNames), drop=F) +
    theme(legend.text=element_text(size=25), legend.title=element_text(size=25)) +
    guides(fill = guide_legend(override.aes = list(size = 10)))
  g_dummy2 <- ggplot(subset(dummy_df, Direction=="Mixed"), aes(Clusters, Genesets,fill=FDR_bin)) + geom_tile() +
    scale_fill_manual(name = "FDR Mixed", values=setNames(cols_mix, interNames), drop=F) +
    theme(legend.text=element_text(size=25), legend.title=element_text(size=25)) +
    guides(fill = guide_legend(override.aes = list(size = 10)))
  g_dummy3 <- ggplot(subset(dummy_df, Direction=="Up"), aes(Clusters, Genesets,fill=FDR_bin)) + geom_tile() +
    scale_fill_manual(name = "FDR Up", values=setNames(cols_up, interNames), drop=F) +
    theme(legend.text=element_text(size=25), legend.title=element_text(size=25)) +
    guides(fill = guide_legend(override.aes = list(size = 10)))
  
  ## Extract legends
  legend1 <- g_legend(g_dummy1)
  legend2 <- g_legend(g_dummy2)
  legend3 <- g_legend(g_dummy3)
  
  ## Arrange legends properly for plotting
  ## From: 'https://stackoverflow.com/questions/50191091/trouble-with-ggplot2-multiple-legend-alignment-with-grid-arrange'
  leg1 <- legend1$grobs[[1]]
  leg2 <- gtable_add_rows(legend2, pos = nrow(legend2) - 1, heights = sum(leg1$heights))
  leg2 <- gtable_add_grob(leg2, leg1, t = nrow(leg2) - 1, l = 3)
  leg3 <- gtable_add_rows(legend3, pos = nrow(legend3) - 1, heights = sum(leg2$heights))
  leg3 <- gtable_add_grob(leg3, leg2, t = nrow(leg3) - 1, l = 3)
  
  ## Save plot to output file
  ## Work around bug that opens a new plot device
  pdf(file=NULL)
  ## Arrange main plot and 3 legends
  ggsave(outfile, plot=arrangeGrob(p, leg3, widths=c(16,4)), width=20, height=20)
  garbage <- dev.off()  
}
