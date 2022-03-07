# This script is given the drug-sets targeting DE genes of each cluster in specific format preprocessed for Upset
# (Input file of type drug,cluster1,cluster2..clustern; tab separated)
# It plots the overlaps of those drug-sets in an Upset plot and returns a png of the plot

# Anne Richter, April 2018


library(UpSetR)
library(optparse)
library(ggplot2)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# parse command line arguments
option_list = list(
  make_option("--inFile", type = "character", help = "/path/to/inputTable of type drug,cluster1..clustern tab separated"),
  make_option("--outFile", type = "character", help = "/path/to/outputPlot enter full absolute path, file must have ending '.png'"),
  make_option("--nsets", type = "numeric", help = "number of sets (e.g. clusters) displayed in plot, default is 100", default = "100"),
  make_option("--nintersects", type = "numeric", help = "number of intersects displayed in plot, default is 20", default = "20")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# exit if outFile has no .png ending
if (!(grepl(".png$", opt$outFile, perl=TRUE))){
  print_help(opt_parser)
  stop("Output file must have '.png' file ending")
}

# exit if header seems to be corrupted
druglist <- read.csv(opt$inFile, header=T, sep="\t" )
header <- names(druglist)
if (!(grepl("^drug$", header[1], perl=TRUE))){
  print_help(opt_parser)
  stop("Header of input file seems to be missing")
}

names(druglist) <- gsub("(\\d+)", "cluster_\\1", names(druglist))
names(druglist) <- gsub("X", "", names(druglist))
print("head(druglist):")
print(head(druglist))


# read command line options into variables
opt_nintersects <- opt$nintersects
opt_nsets <- opt$nsets

# if less than two malignant clusters are predicted, generate substitute bar plot instead of upsetR plot 
if (ncol(druglist)==2){
  names <- colnames(druglist)
  name_cluster <- names[2]
  print(name_cluster)
  # open PNG file for writing, generate bar plot, close and export PNG
  png(opt$outFile, width = 2200, height = 1800, res = 300)
  plot1 <- ggplot(druglist, aes(x=druglist[,2]))+
    geom_bar()+
    ylab("Number of drugs targeting malignant cell cluster")+
    xlab(name_cluster)+
    scale_x_continuous(breaks = c(1, 0))+
#    scale_x_discrete(breaks = c(1, 0))+
    ggtitle("In this sample only one malignant cell cluster was found.
No Upset plot could be generated.")+
    geom_text(stat='count', aes(label=..count..), vjust=-1)
  print(plot1)
  dev.off()
  #stop("Only one malignant cluster was found!")
}


# if several malignant clusters are found, open PNG plot for writing, generate upset plot, close and export PNG
if (ncol(druglist)>2){
  # open PNG file for writing, generate upset plot, close and export PNG
  png(opt$outFile, width = 2200, height = 1800, res = 300)
  upset(data=druglist, nsets=opt_nsets, nintersects=opt_nintersects, order.by=c("freq"), 
                   #intersections = list("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7"),
                   matrix.color = 3,        # color of dots and lines
                   shade.color = 3,        # shading in table
                   shade.alpha = 0.25,        # Transparency of shading in matrix, 0 is no shading
                   keep.order = FALSE,
                   main.bar.color = "grey20",
                   line.size = 0.5,        # lines connecting dots
                   point.size = 3.5,
                   mainbar.y.label = "Number of drugs",
                   #mainbar.y.max = "100" #?
                   sets.x.label = "Number of drugs targeting cluster",        # label of sets barplots
                   #att.color = 5,    #?
                   #group.by = "degree", #?
                   matrix.dot.alpha = 1,        # transparency of dots not in set, 1 is strongest
                   sets.bar.color = "grey38",
                   mb.ratio = c(0.7, 0.3),        # ratio between main bar plot and matrix plot
                   text.scale = c(1.5,1.5,1,1,1.5,1.5)
                   # Can be a universal scale, or a vector containing individual scales in the following format: 
                   # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
                   )
  dev.off()
}


# if the input file has less than two columns something is wrong. A warning is written into the output .png file.
if (ncol(druglist)<2){
  png(opt$outFile, width = 2200, height = 1800, res = 300)
  plot1 <- ggplot()+
    ggtitle("Warning! The input file has less than two columns and looks therefore different than expected.
Something went wrong during the analysis or no cluster was classified as 'malignant'.")+
    theme(plot.title = element_text(size=12))
  print(plot1)
  dev.off()
}
