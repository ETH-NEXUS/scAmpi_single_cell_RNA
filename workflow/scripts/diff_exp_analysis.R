################################################################################
#### File name: diff_exp_analysis.R
#### Author: Anne Bertolini
#### Co-author: Michael Prummer
#### Date created: February 2020
#### R Version: 4.0
################################################################################

# Run linear model  on single TP sample,
# Compare each malignant cluster against the average of all non-malignant clusters of the same sample.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(ggplot2)
  library(optparse)
  library(multcomp)
  library(broom)
  library(reshape2)
})

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

option_list = list(
  make_option("--sample_data", type = "character", help = "Path to RDS input file that contains SCE object."),
  make_option("--sampleName", type = "character", help = "Character string that will be prefix to all output files"),
  make_option("--cluster_table", type = "character", help = "Table with overview over cluster and the cell types of the cells in each cluster. '*.phenograph_celltype_association.txt'"),
  make_option("--malignant_tag", type = "character", help = "Name or substring of malignant cell type"),
  make_option("--threshold_comparison", type = "character", help = "Minumum number of cells with a cell type to be tested against the malignant clusters for DE genes."),
  make_option("--fdr_cut", type = "character", help = "Cut-off for p-value."),
  make_option("--fc_cut", type = "character", help = "Cut-off for fold change."),
  make_option("--mindiff2second", type = "character", help = "Minimum difference in mean expression value between malignant cluster and highest/lowest non-malignant cluster."),
  make_option("--minNumberNonMalignant", type = "character", help = "Minimum number of non-malignant clusters in the sample for DE malignant vs. 'rest' to be performed."),
  make_option("--outdir", type = "character", help = "Full path to output directory")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# convenience function for string concatenation
'%&%' = function(a,b) paste(a,b,sep="")


# Input parameters
fdr_cut <- as.numeric(opt$fdr_cut)
fc_cut <- as.numeric(opt$fc_cut)
mindiff2second <- as.numeric(opt$mindiff2second)
threshold_comparison <- as.numeric(opt$threshold_comparison)
minNumberNonMalignant <- as.numeric(opt$minNumberNonMalignant)

######################   Load SCE object  #####################

# Read in SCE object from sample
my_sce <- readRDS(opt$sample_data)
print(my_sce)
writeLines("\n")

###############   Read in cluster table with dominant cell type   ##############

# Read in cell type composition of phenograph clusters.
ct_composition <- read.csv(opt$cluster_table, header=TRUE,
                           stringsAsFactors=FALSE, sep = "\t")

# Get the dominant cell types of the clusters into a list.
# Here, it is assumed that there is no cluster 0 present in the data!
dominant_types <- as.list(ct_composition$Dominant.celltype)
number_clusters <- length(dominant_types) - 2
length(dominant_types) <- number_clusters
dominant_types <- unlist(dominant_types)
print("dominant_types:")
print(dominant_types)

# Sort the cluster IDs to the respective dominant cell type.
cluster_ids <- as.list(ct_composition$Cluster)
length(cluster_ids) <- number_clusters
cluster_ids <- unlist(cluster_ids)
names(dominant_types) <- cluster_ids
print("str(dominant_types):")
print(str(dominant_types))

# remove all clusters that are too small
cluster_sizes <- table(as.numeric(as.character(my_sce$phenograph_clusters)))
too_small_mask <- cluster_sizes < threshold_comparison
clusters_too_small <- names(too_small_mask)[too_small_mask, drop = F]
cat("\n\n\nWarning : cluster(s)", clusters_too_small, "is/are too small and will not be included in the DE analysis\n")
cat("dominant celltype(s) of small cluster(s):", dominant_types[too_small_mask], "\n\n\n")
dominant_types <- dominant_types[!names(dominant_types) %in% clusters_too_small, drop = FALSE]

# Get all malignant clusters of the current sample.
malignant_mask <- grepl(opt$malignant_tag, dominant_types,
                        ignore.case = TRUE, perl=TRUE)
print("Print malignant mask table: FALSE = non-malignant, TRUE = malignant")
print(table(malignant_mask))
malignant_clusters <- dominant_types[malignant_mask]
mal_clusterIDs <- names(malignant_clusters)
print("mal_clusterIDs:")
print(mal_clusterIDs)
non_mal_clusters <- dominant_types[!malignant_mask]
non_mal_clusterIDs <- names(non_mal_clusters)
print("non_mal_clusterIDs:")
print(non_mal_clusterIDs)

# column number of malignant cells in ColData
malig_colid = which(malignant_mask)
nonmal_colid = which(!malignant_mask)

# Get all non-malignant cell types of the cohort
# nonmal_cts <- sort(as.character(unique(my_sce$celltype_final)))
# # Remove cell types from the list with the malignant tag.
# mask_nonmal <- grepl(opt$malignant_tag, nonmal_cts, ignore.case = TRUE, perl=TRUE)
# nonmal_cts <- nonmal_cts[!mask_nonmal]
# mask_unc <- grepl("uncertain|unknown", nonmal_cts, ignore.case = TRUE, perl=TRUE)
# nonmal_cts <- nonmal_cts[!mask_unc]

################################   DE analysis   ###############################

######################################################################################################
## MIPR
###############   DEA of Pearson Residuals one cluster against average of all others   ###############
## Initiate list with DE results.
list_res = list()
if (length(malignant_clusters) > 0) {
  # use SCE object
  mresid = assay(my_sce, "pearson_resid", min.cells = 0, min.features = 0)
  mnorm = assay(my_sce, "normcounts", min.cells = 0, min.features = 0)
  cat("\n\n\n Before excluding small clusters:\n")
  print("dim(mnorm):")
  print(dim(mnorm))
  print("dim(mresid):")
  print(dim(mresid))
  # get cluster ids
  clu = sprintf("c%02.0f", as.numeric(as.character(my_sce$phenograph_clusters)))
  names(clu) = my_sce$barcodes
  print("str(clu):")
  print(str(clu))
  # exclude small clusters
  clu.size = table(clu)
  print("clu.size:")
  print(clu.size)
  print("str(threshold_comparison):")
  print(str(threshold_comparison))
  clu.exclude = which(clu.size < threshold_comparison)
  print("str(clu.exclude):")
  print(str(clu.exclude))
  if(length(clu.exclude)>0) {
    idx = which(clu %in% names(clu.exclude))
    print("str(idx):")
    print(str(idx))
    clu = clu[-idx]
    mresid = mresid[, -idx]
    mnorm = mnorm[, -idx]
  }
  #print("str(clu):")
  cat("\n\n\n After excluding small clusters:\n")
  print("dim(mresid):")
  print(dim(mresid))
  print("dim(mnorm):")
  print(dim(mnorm))
  fit.model = formula("y ~ clu")
  # Find DE genes with regular linear model.
  dea_fit = function(y, clu, fit.model){
    dlm = try(lm(fit.model, data=data.frame(y=y, clu=clu[names(y)])), silent=T)
  }
  # Model of expression vs cluster_id for all genes and all clusters
  fitobject = apply(mresid, 1, dea_fit, clu, fit.model)
  # cl = makeCluster(getOption("cl.cores", 7))
  # fitobject = parApply(cl=cl, mresid, 1, dea_fit, clu, fit.model)
  # stopCluster(cl)
  print("head(sapply(fitobject, class))")
  print(head(sapply(fitobject, class)))
  # y = mresid[1, 1:100]
  # idx = match(names(y), colnames(mresid))
  # fitobject = lm(fit.model, data=data.frame(y=y, clu=clu[idx]))
  
  dea_results  = function(fitobject, contrastmatrix){
    if(class(fitobject)[1] != "try-error" & all(!is.na(coef(fitobject)))){
      t.cm = contrastmatrix[, names(coef(fitobject)), drop=F]
      idy = apply(contrastmatrix, 1, function(x) names(which(x!=0)))
      #idy = sapply(idy, function(x) all(x %in% names(coef(fitobject))))
      tout = try(glht(fitobject, t.cm), silent=T)
      tout = try(tidy(summary(tout, type="none")), silent=T)
      if(class(tout)[1] != "try-error") {
        names(tout)[6] = "p.value"
        tout[, -2]
      } else {
        return("-1")
      }
    } else {
      return("-1")
    }
  }
  
  malig_clusterid = as.numeric(mal_clusterIDs)
  nonmal_clusterid = as.numeric(non_mal_clusterIDs)

  # only perform DE malignant vs. "rest" if a minimum number of non-malignant clusters was found
  if(length(non_mal_clusters) >= minNumberNonMalignant){
    writeLines("\n")
    print("###### DEA between all malignant cluster and the average over all non-malignant clusters ######")
    # build contrast matrix for malig vs rest
    contrastmatrix = matrix(0, nrow=length(malig_clusterid), ncol=length(clu.size))
    colnames(contrastmatrix) = "clu" %&% names(clu.size)
    rownames(contrastmatrix) = sprintf("c%02.0f.vs.rest", as.numeric(names(malignant_clusters)))
    id = malig_colid + seq(nrow(contrastmatrix)) - 1 + (malig_colid-1)*(nrow(contrastmatrix)-1) 
    contrastmatrix[id] = 1
    contrastmatrix[, nonmal_colid] = -1/length(nonmal_colid)
    colnames(contrastmatrix)[1] = "(Intercept)"
    print("contrastmatrix:")
    print(contrastmatrix)
    # DEA contrasts between each malignant cluster and the average over all nonmalignant clusters.
    list_res = c(list_res, lapply(fitobject, dea_results, contrastmatrix))
  }

  # if present, compare current malignant cluster with the average overr all other malignant clusters
  if(length(malignant_clusters) > 1){
    writeLines("\n")
    print("###### DEA between each malignant cluster and the average over all other malignant clusters ######")
    # build contrast matrix for malig vs other malig clusters
    malig_cm = matrix(0, nrow=length(malig_clusterid), ncol=length(clu.size))
    colnames(malig_cm) = "clu" %&% names(clu.size)
    rownames(malig_cm) = sprintf("c%02.0f.vs.malig", as.numeric(names(malignant_clusters)))
    malig_cm[, malig_colid] = -1/(length(malig_colid)-1)
    id = malig_colid + seq(nrow(malig_cm)) - 1 + (malig_colid-1)*(nrow(malig_cm)-1) 
    malig_cm[id] = 1
    colnames(malig_cm)[1] = "(Intercept)"
    print("malig_cm:")
    print(malig_cm)
    # DEA contrasts between each malignant cluster and the average over all other malignant clusters.
    list_res = c(list_res, lapply(fitobject, dea_results, malig_cm))
  } 


  df_res = melt(list_res, id.var=1:5)
  attr(df_res$p.value, "error") = NULL 
  names(df_res)[c(1,6)] = c("comparison", "gene_names")
  print("str(df_res)")
  print(str(df_res))
  
  df_res$padj = p.adjust(df_res$p.value, method="BH")
  df_res$padj[df_res$padj==0] = min(df_res$padj[df_res$padj > 0])
  df_res$malig = gsub("(^.*)\\.vs.*", "\\1", df_res$comparison)
  
  # rename headers for backwards consistency
  print("str(df_res) before renaming the column headers:")
  print(str(df_res))
  colnames(df_res) <- gsub("estimate", "diff", colnames(df_res))
  colnames(df_res) <- gsub("statistic", "test_statistic", colnames(df_res))
  print("str(df_res) after renaming the column headers:")
  print(str(df_res))

  ggplot(df_res, aes(p.value, ..density..)) + geom_histogram(bins=100, fill="seagreen") + 
    geom_hline(yintercept = 1, col="black") + xlab("P-value") + ylab("Density") + 
    scale_fill_brewer(palette="Dark2") + theme(legend.position="none") + 
    coord_cartesian(ylim=c(0, 5)) + facet_wrap(~comparison, ncol=2) + theme_bw()
  ggsave(opt$outdir %&% opt$sampleName %&% ".pvd_pearson_residuals.png",
         width = 30, height = 20, dpi = 600, units = "cm")
  

  # Volcano plot
  t.plot = df_res
  t.plot$hit = c("not","signif")[as.numeric(t.plot$padj < fdr_cut)+1]
  t.plot$hit[abs(t.plot$diff) > fc_cut & t.plot$hit== "signif"] = "signif & strong"
  t.plot$hit[abs(t.plot$diff) > fc_cut & t.plot$hit== "not"] = "strong"
  t.plot$hit = factor(t.plot$hit, levels=c("not", "strong", "signif", "signif & strong"))
  cbPalette = c("#000000", "#009E73", "#56B4E9", "#D55E00", 
                "#F0E442", "#0072B2", "#E69F00", "#CC79A7")
  ggplot() +
    geom_point(data = t.plot, aes(x=diff, y=-log10(padj), color=hit), shape=16, size=1) + 
    geom_hline(yintercept = -log10(fdr_cut)) + 
    geom_vline(xintercept = c(-1,1)*fc_cut) +
    xlab(expression(paste("Effect size (model coefficient)"))) + 
    ylab(expression(paste("Significance [", -log[10](padj), "]"))) +
    scale_color_manual(values = cbPalette, name="Hit class", drop=F) +
    coord_cartesian(xlim=c(-5,5)) + theme_bw() +
    guides(colour = guide_legend(override.aes = list(size=2, shape=15))) +
    theme(legend.position="top") + facet_wrap(~comparison, ncol=2)  
  
  # one volcano plot per comparison "malig.vs.rest"
  tt.plot = t.plot[grep("rest", t.plot$comparison), ]
  by(tt.plot, tt.plot$comparison, function(x){
    ggplot(x, aes(x=diff, y=-log10(padj), color=hit)) +
      geom_point(shape=16, size=1) + 
      geom_hline(yintercept = -log10(fdr_cut)) + 
      geom_vline(xintercept = c(-1,1)*fc_cut) +
      xlab(expression(paste("Effect size (model coefficient)"))) + 
      ylab(expression(paste("Significance [", -log[10](padj), "]"))) +
      scale_color_manual(values = cbPalette, name="Hit class", drop=F) +
      coord_cartesian(xlim=c(-5,5)) + theme_bw() +
      guides(colour = guide_legend(override.aes = list(size=2, shape=15))) +
      theme(legend.position="top")
    malig = as.numeric(substr(x$malig[1], 2, 3))
    fname = sprintf("%s%s.%i.volcanoPlot.png", opt$outdir, opt$sampleName, malig)
    ggsave(fname, width = 30, height = 20, dpi = 600, units = "cm")
    }
  )
  
  # per cluster average for each gene:
  cluster_mean = t(apply(mnorm, 1, function(x) tapply(x, clu, mean))) # gene x cluster matrix
  # What is the cluster mean for each gene in the malignant cluster relevant for each comparison?
  malig_cluster_mean = cluster_mean[, malig_colid, drop=F]
  malig_cluster_mean = melt(malig_cluster_mean)
  names(malig_cluster_mean) = c("gene_names", "malig", "malig_mean")
  print("str(malig_cluster_mean)")
  print(str(malig_cluster_mean))
  # merge with df_res
  df_res = merge(df_res, malig_cluster_mean, sort=F)

  # What is the maximal (minimal) cluster mean for each gene across all nonmalignant clusters?
  nonmal_mean_minmax = data.frame(gene_names=rownames(mnorm), nonmal_min = NA, nonmal_max = NA)
  if(length(nonmal_colid) >= minNumberNonMalignant){
    nonmal_cluster_mean = cluster_mean[, nonmal_colid, drop=F]
    nonmal_mean_minmax$nonmal_min = rowMin(nonmal_cluster_mean)
    nonmal_mean_minmax$nonmal_max = rowMax(nonmal_cluster_mean)
  }
  # merge with df_res
  df_res = merge(df_res, nonmal_mean_minmax, sort=F)
  # for each gene, the percent non-zero cells for the target malignant cluster
  pct_nonzero = t(apply(mnorm, 1, function(x) tapply(x, clu, function(x) 100*length(which(x > 0))/length(x))))
  pct_nonzero = pct_nonzero[, malig_colid, drop=F]
  pct_nonzero = melt(pct_nonzero)
  names(pct_nonzero) = c("gene_names", "malig", "pct_nonzero")
  print("str(pct_nonzero)")
  print(str(pct_nonzero))
  df_res = merge(df_res, pct_nonzero, sort=F)

  print("str(df_res)")
  print(str(df_res))
  
  
  # Write out table with results.
  # output_filename = opt$outdir %&% opt$sampleName %&% ".DEgenes.tsv"
  # write.table(df_res, output_filename, quote = FALSE, sep = "\t", row.names = FALSE)
  
  # Write to disk malig.vs.rest
  df_out = split(df_res, grepl("rest", df_res$comparison))
  print("str(df_out):")
  print(str(df_out))
  if(length(non_mal_clusters) >= minNumberNonMalignant){
    by(df_out[["TRUE"]], df_out[["TRUE"]]$malig, function(x){
      malig = as.numeric(substr(x$malig[1], 2, 3))
      output_filename = sprintf("%s%s.%i.DEgenes.tsv", opt$outdir, opt$sampleName, malig)
      write.table(x, output_filename, quote = FALSE, sep = "\t", row.names = FALSE)
      print("Result written to file: <" %&% output_filename %&% ">.")
      writeLines("\n")
    })
  }
  if(length(malignant_clusters) > 1){
    outputPath_vsMalignant = paste(opt$outdir,"vs_other_malignant/",sep='')
    dir.create(outputPath_vsMalignant, showWarnings = FALSE)
    by(df_out[["FALSE"]], df_out[["FALSE"]]$malig, function(x){
      # x=df_out[["FALSE"]][df_out[["FALSE"]]$malig==df_out[["FALSE"]]$malig[1], ]
      malig = as.numeric(substr(x$malig[1], 2, 3))
      output_filename = sprintf("%svs_other_malignant/%s.DEmalignant.%i.DEgenes.tsv", opt$outdir, opt$sampleName, malig)
      write.table(x, output_filename, quote = FALSE, sep = "\t", row.names = FALSE)
      print("Result written to file: <" %&% output_filename %&% ">.")
      writeLines("\n")
    })
  }
  
  ################################################################################
  ##### Plot the Pearson residuals for the top 8 upregulated and 
  ##### the top 8 downregulated genes (with respect to any comparison).
  ################################################################################
  # filter for the desired target genes:
  idx = which(df_res$malig_mean < df_res$nonmal_min - mindiff2second & df_res$padj < fdr_cut  |
                df_res$malig_mean > df_res$nonmal_max + mindiff2second & df_res$padj < fdr_cut)
  if(length(idx)>0){
    # expression plot of top 16 genes of interest.
    goi = df_res[idx, c("gene_names", "diff"), drop=F]
    goi = unique(goi$gene_names[order(goi$diff, decreasing = T)])
    # for the plot, select half up-regulated, half down-regulated.
    goi = c(head(goi, 8), tail(goi, 8))
    goi <- unique(goi)
    t.plot = mresid[rownames(mresid) %in% goi, , drop=F]
    t.plot = melt(t.plot)
    names(t.plot) = c("gene_names", "barcodes", "Pearson_residual")
    t.plot$gene_names = factor(as.character(t.plot$gene_names), levels=goi)
    t.plot = merge(t.plot, as.data.frame(colData(my_sce)), sort=F)
    
    ggplot(t.plot, aes(as.factor(phenograph_clusters), Pearson_residual)) + 
      geom_boxplot() + facet_wrap(~gene_names, nrow=4) +
      coord_cartesian(ylim=c(-1,15)) + 
      xlab("Cluster ID") + theme_bw()
    ggsave(opt$outdir %&% opt$sampleName %&% ".candidate_genes_pearson_residuals.png",
           width = 30, height = 20, dpi = 600, units = "cm")
  }
}

