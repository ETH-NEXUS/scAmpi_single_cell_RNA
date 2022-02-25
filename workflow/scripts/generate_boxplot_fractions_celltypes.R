###################################################
#### File name: generate_boxplot_fractions_celltypes.R
#### Author: Anne Bertolini
#### Date created: August 2019
#### R Version: 3.5.1
###################################################

# make boxplot with cell type fractions of all samples with the current sample highlighted
library(optparse)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# parse command line arguments
option_list = list(
make_option("--current_sample", type = "character", help = "Path to the *phenograph_celltype_association.txt table of the current sample."),
make_option("--previous_samples", type = "character", help = "Path to text file that contains the paths to the *phenograph_celltype_association.txt tables of previous samples. One path per line."),
make_option("--sampleName", type = "character", help = "Sample name that will be added to the names of all output files."),
make_option("--sampleName_short", type = "character", help = "Short version of the sample name that is given to the cellranger run."),
make_option("--outDir", type = "character", help = "Full path to output directory.")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


print(opt$previous_samples)
print("Current sample:")
print(opt$sampleName_short)
ct_tables <- read.csv(opt$previous_samples, header = FALSE, stringsAsFactors = FALSE)
print("ct_tables$V1, all cohort samples:")
print(ct_tables$V1)

# check if current sample is already in list and remove it
mask_current_sample <- grepl(opt$sampleName_short, ct_tables$V1)
print("Check if current sample is in cohort:")
print(table(mask_current_sample))
ct_tables <- as.list(ct_tables$V1[!(mask_current_sample)])
print("str(ct_tables), current sample is removed if it was already in the cohort list:")
print(str(ct_tables))

# add current sample to list
ct_tables <- c(ct_tables, opt$current_sample)
print("str(ct_tables) after adding current sample separately:")
print(str(ct_tables))

list_tables <- lapply(seq_along(ct_tables), function(x){
  fraction_table <- read.table(ct_tables[[x]],stringsAsFactors = FALSE)
  #print(head(fraction_table))
  trans <- as.data.frame(t(fraction_table[1,]))
  names(trans) <- "celltypes"
  last_col <- length(fraction_table[,1])
  percent <- t(fraction_table[last_col,])
  trans$percentage <- percent
  trans <- trans[-c(1),]
  trans <- head(trans, -5)
  trans$percentage_cleaned <- gsub(pattern = "(\\S+)\\s+.+", "\\1", trans$percentage)
  trans$percentage_numeric <- as.numeric(trans$percentage_cleaned)
  trans$ID <- unlist(rep(basename(ct_tables[[x]]), length(trans$percentage_numeric)))
  trans
})

all_samples <- do.call(rbind, list_tables)
#all_samples$location <- as.factor(gsub(pattern = "MAHACEB.*|MATIWAQ.*|MIGOFIW.*|MIJUCYK.*", "subcutaneous", all_samples$ID))
#all_samples$location <- as.factor(gsub(pattern = "MAPOXUB.*", "brain", all_samples$location))
#all_samples$location <- as.factor(gsub(pattern = "MENYTEK.*|MEXUXEH.*|MIDOBOL.*", "lymph_node", all_samples$location))
#all_samples$location <- as.factor(gsub(pattern = "MIKYCYN.*", "lung", all_samples$location))
#all_samples$location <- as.factor(gsub(pattern = "MOFYCAG.*", "paracardial", all_samples$location))
print("str(all_samples):")
print(str(all_samples))

# mask for sample of interest
mask_current_sample_only <- grepl(opt$sampleName, all_samples$ID)
current_sample_only <- all_samples[mask_current_sample_only,]
# For showing the outliers in the boxplot:
#maxVal = 90
# all other outliers
#labelDD = subset(all_samples,percentage_numeric>maxVal & !(ID == current_sample_only$ID))
#outlier_txt=paste(format(round(labelDD$percentage_numeric, 0), nsmall = 1),collapse=", ")
# outliers of current sample
#label_current_outlier <- subset(all_samples, percentage_numeric>maxVal & ID == current_sample_only$ID)
#current_outlier_txt=paste(format(round(label_current_outlier$percentage_numeric, 0), nsmall = 1),collapse=", ")
# sort cell types alphabetically for the plot
all_samples <- all_samples[order(all_samples$celltypes),]

boxplot_melanoma = ggplot(all_samples, aes(x = celltypes, y = percentage_numeric)) +
  geom_boxplot(notch=FALSE, fill = "grey80") +
#  geom_boxplot(outlier.colour="red", outlier.shape=5, outlier.size=1, notch=FALSE, fill = "grey80") +
  geom_boxplot(notch=FALSE, fill = "grey80") +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=7, binwidth = 0.1, fill = "black") +
  geom_dotplot(data = current_sample_only, binaxis='y', stackdir='center', dotsize=18, binwidth = 0.1, fill = "red", colour = "red") +
  labs(title="Fractions of cell types from previous samples, current sample highlighted",x="cell types", y = "percentage") +
  theme(legend.position="none")+
  theme(axis.title.x = element_text(face="bold", size=12),
        axis.text.x  = element_text(angle=50, size=12,  hjust = 1, face = "bold"),
        axis.title.y = element_text(face="bold", size=12))
#  coord_cartesian(ylim = c(0, maxVal+3), expand = TRUE, default = FALSE, clip = "on") +
#  scale_y_continuous(limits=c(0, maxVal+3)) +
#  geom_text(data=labelDD,aes(y=maxVal,label=outlier_txt),size=3.5,vjust=-0.5,hjust=0.5) +
#  geom_text(data=labelDD,aes(y=maxVal,label=current_outlier_txt),size=4.3,vjust=-1.8,hjust=0.5, colour = "red", fontface = "bold") +
#  geom_segment(data=labelDD,aes(y=maxVal*0.97,yend=maxVal, xend=factor(celltypes)), arrow = arrow(length = unit(0.1,"cm")))

boxplot_melanoma

filename = paste0(opt$outDir, opt$sampleName, ".boxplot_cell_types_cohort.png")
print("Output file:")
print(filename)
ggsave(filename = filename, width = 30, height = 20,  dpi = 600, units = "cm")
