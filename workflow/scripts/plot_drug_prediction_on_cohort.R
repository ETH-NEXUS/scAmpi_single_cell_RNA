# #################################################

#### File name: show_drugPrediction_on_cohort.R
#### Author: Lourdes Rosano
#### Date created: June 2020
#### R Version: 3.5.1

# #################################################

################################################################################
## functions
################################################################################

retrieve_malig_ids = function(files, drugPredEnd){
  
  malig_ids = sapply(files, function(x) {
    firstSplit = strsplit(x, drugPredEnd)[[1]]
    secondSplit = strsplit(firstSplit, "\\.")[[1]]
    secondSplit[length(secondSplit)]
  })
  malig_ids = suppressWarnings(as.integer(malig_ids))
  
  ## Sanity check for finding integer malignant ids
  stopifnot(sum(is.na(malig_ids))==0)
  # NOTE: it is critical that cluster ids are of type="character" for selecting objects based on cluster id and not on position
  malig_ids = as.character(malig_ids) # sanity check
  
  return(malig_ids)  
}


process_dgidb_single_sample = function(all_clin, name_DGIDB){
  
  ## Sanity check for finding DGIDB column name in all clinical files for a given sample
  notFound_dgidb = FALSE
  for (i in seq_along(all_clin)){
    file = all_clin[[i]]
    # Retrieve DGIDB drug prediction column within each clinical annotations file
    tmp_dgidb = which(names(file) %in% name_DGIDB)
    if (length(tmp_dgidb)==0){ notFound_dgidb = TRUE }
  }
  ## Check if column was not found for some file
  stopifnot(identical(notFound_dgidb,FALSE))
  
  
  #### Process DGIDB drug prediction data
  
  ## Extract DGIDB drug predictions from clinical annotation files (tumor clusters only)
  # The same drug can occur several times for the same cluster if the associated DGIDB score was different depending on the gene involved
  # DGIDB does not make predictions for drug combinations
  all_dgidb_data = sapply(all_clin, function(x){
    # Retrieve drug prediction column within each clinical annotations file
    colIndx = which(names(x) %in% name_DGIDB)
    cluster_drugs = x[,colIndx,drop=F]
    drugMat_dgidb = data.frame(Drug=character(), Support=character())
    # Retrieve available drug predictions (if any)
    # When at least 1 available, process annotations further
    if (nrow(cluster_drugs)>0){
      # Split multiple annotations by ";" (if any)
      # Eg: HEXADECANESULFONYL FLUORIDE(4,.);PALMITIC ACID(4,.);
      drugArr = as.character(unlist(apply(cluster_drugs, 1, function(y) strsplit(y, ";")[[1]])))
      # Retrieve relevant annotations (drug name + score, ignore type info) for each malignant cluster
      # NOTE: Careful with drug names containing "(" in their names!
      drugMat_dgidb = data.frame(matrix(
        sapply(drugArr, function(z){
          tmp_row = strsplit(z, "\\(")[[1]]
          drug = tmp_row[1]
          tmp_score = tmp_row[length(tmp_row)]
          # Make sure to remove trailing and leading spaces in the score
          score = trimws(strsplit(tmp_score, "\\,")[[1]][1])
          # Rejoin drug names containing "(" (if any)
          if (length(tmp_row)>2){
            drug = paste(tmp_row[1:(length(tmp_row)-1)], collapse="(")
          }
          # Make sure to remove trailing and leading spaces in the drug name as well
          # Also, turn to uppercase to avoid mismatching drugs due to case differences
          drug = toupper(trimws(drug))
          c(drug,score)
        }), 
        ncol=2, byrow=T))
      names(drugMat_dgidb) = c("Drug","Score")
    }
    # Remove duplicated pairs of drug + score (if any)
    drugMat_dgidb = unique(drugMat_dgidb)
  }, simplify=FALSE)
  
  # Extract DGIDB predicted drugs for each detected tumor cluster
  all_tumor_drugs_dgidb = sapply(names(all_dgidb_data), function(id) {
    drugMat_dgidb = data.frame(Drug=character(), phenograph_clusters=character()) 
    tumor_drugs = unique(all_dgidb_data[[id]]$Drug)
    if (length(tumor_drugs)>0){
      drugMat_dgidb = data.frame(Drug=tumor_drugs, phenograph_clusters=id)  
    }
    drugMat_dgidb
  }, simplify=FALSE)
  all_tumor_data_dgidb = Reduce(rbind, all_tumor_drugs_dgidb)
  
  return(all_tumor_data_dgidb)
}


translate_civic_drugs = function(civic_mat, civic_dict){
  
  all_civic_drugs = civic_mat$Drug
  # CIVIC drugs should only contain "+" when they correspond to a combination treatment
  isCombi = grepl("\\+",all_civic_drugs)
  
  civic_mat_single = civic_mat
  
  ## TODO: also translate to synonyms within CIVIC combination drugs
  ## TODO: test all single and all combi
  ## TODO: make sure that this still works for no civic drug prediction (nrow=0)
  
  ## Split civic drug matrix into 2 parts: single drugs and combi drugs
  # civic_mat_single = civic_mat[!isCombi,,drop=F]
  # civic_mat_combi = civic_mat[isCombi,,drop=F]
  # stopifnot(nrow(civic_mat_single)+nrow(civic_mat_combi)==nrow(civic_mat))
  
  if (nrow(civic_mat_single)>0){
    all_single_drugs = civic_mat_single$Drug
    to_translate = sapply(all_single_drugs, function (x) x %in% civic_dict[,1])
    to_translate = unique(names(to_translate[to_translate]))
    # Will only enter the loop if there is at least 1 drug to translate in civic_mat_single
    # Otherwise, civic_mat_single remains unchanged (nothing is translated)
    for (z in seq_along(to_translate)){
      this_drug = to_translate[z]
      this_syn = subset(civic_dict, subset = CIVIC == this_drug, select=Synonym)[,1]
      civic_mat_single$Drug[civic_mat_single$Drug %in% this_drug] = this_syn
    }
  }
  
  ## TODO:
  # if (nrow(civic_mat_combi)>0){
  # }
  
  new_civic_mat = civic_mat_single
  ## TODO:
  # new_civic_mat = rbind(civic_mat_single, civic_mat_combi)
  stopifnot(nrow(new_civic_mat)==nrow(civic_mat))
  
  return(new_civic_mat)
}


process_civic_single_sample = function(all_clin, name_CIVIC, civic_dict){
  
  ## Sanity check for finding CIVIC column name in all clinical files for a given sample
  notFound_civic = FALSE
  for (i in seq_along(all_clin)){
    file = all_clin[[i]]
    # Retrieve CIVIC drug prediction column within each clinical annotations file
    tmp_civic = which(names(file) %in% name_CIVIC)
    if (length(tmp_civic)==0){ notFound_civic = TRUE }
  }
  ## Check if column was not found for some file
  stopifnot(identical(notFound_civic,FALSE))
  
  #### Process CIVIC drug prediction data
  
  ## Extract CIVIC drug predictions from clinical annotation files (tumor clusters only)
  # Opposite to DGIDB, CIVIC makes predictions for drug combinations (separated with "+")
  all_civic_data = sapply(all_clin, function(x){
    # Retrieve drug prediction column within each clinical annotations file
    colIndx = which(names(x) %in% name_CIVIC)
    cluster_col = x[,colIndx,drop=F]
    # Remove rows where there was no CIVIC prediction (annotated as '.')
    cluster_drugs = cluster_col[cluster_col[,1]!=".",,drop=F]
    drugMat_civic = data.frame(Drug=character(), Support=character())
    # Retrieve available drug predictions (if any)
    # When at least 1 available, process annotations further
    if (nrow(cluster_drugs)>0){
      # Split multiple annotations by ";" (if any)
      # Eg: DOCETAXEL:NCT:CIVIC_RESISTANCE;TOPOTECAN+CARBOPLATIN+CYCLOPHOSPHAMIDE:CT:CIVIC_SUPPORT
      drugArr = as.character(unlist(apply(cluster_drugs, 1, function(y) strsplit(y, ";")[[1]])))
      # Retrieve relevant annotations (drug name + support, ignore ct info) for each malignant cluster
      # Also, alphabetically order individual drugs of combination treatments to avoid mismatches
      drugMat_civic = data.frame(matrix(
        sapply(drugArr, function(z){
          tmp_row = strsplit(z, ":")[[1]][c(1,3)]
          # Reorder original drugs for combination treatments
          if (grepl("\\+",tmp_row[1])){
            comb = tmp_row[1]
            combArr = strsplit(comb,"\\+")[[1]]
            s.comb = paste(sort(combArr), collapse="+")
            tmp_row[1] = s.comb
          }
          # Turn drug and civic support to uppercase to avoid mismatchings due to case differences
          tmp_row[1] = toupper(tmp_row[1])
          tmp_row[2] = toupper(tmp_row[2])
          tmp_row
        }), 
        ncol=2, byrow=T))
      names(drugMat_civic) = c("Drug","Support")
    }
    # Remove duplicated pairs of drug + civic support (if any), eg. when several genes have the same drug+support associated
    drugMat_civic = unique(drugMat_civic)
  }, simplify=FALSE)
  

  ### Translate CIVIC drug names that have synonyms (if any)
  
  # NOTE: it is critical to do this step before cleaning the civic data, as we could miss drugs with gene-dependent support due to name mismatches
  # This step works even if no CIVIC drug predictions are available (nrow=0)
  all_civic_data_translated = sapply(names(all_civic_data), function(id) {
    mat_civic = all_civic_data[[id]]
    mat_civic_translated = translate_civic_drugs(mat_civic, civic_dict)
    unique(mat_civic_translated) # in case translation causes drug name duplicates
  }, simplify=FALSE)
  
  # Only if no drug was translated will both dataframes be identical
  if (identical(identical(all_civic_data_translated, all_civic_data),FALSE)){
    print("At least one CIVIC drug name was translated to an available synonym.")
    print("Original civic drug predictions:")
    print(all_civic_data)
    print("Translated civic drug predictions:")
    print(all_civic_data_translated)
  }
  
  # Remove drugs associated to non-confident CIVIC evidence: CIVIC_UNKNOWN, CIVIC_CONFLICT
  # Also, mark drugs associated to conflicting CIVIC predictions due to several DE genes: CIVIC_GENE_DEP
  # Again, here it is critical that drug synonyms are already translated to be able to match drug names
  all_civic_data_clean = sapply(names(all_civic_data_translated), function(id) {
    mat_civic = all_civic_data_translated[[id]]
    
    ## 1) Remove rows associated to non-confident evidence (unknown or conflicting)
    keep_rows = !mat_civic$Support %in% c("CIVIC_UNKNOWN","CIVIC_CONFLICT")
    new_mat_civic = mat_civic[keep_rows,,drop=F]
    new_mat_civic = unique(new_mat_civic) # sanity check
    
    ## 2) Mark drugs associated to conflicting predictions as so
    # All remaining CIVIC support should correspond to either CIVIC_RESISTANCE or CIVIC_SUPPORT,
    # so duplicated drug names mean that they are associated to conflicting predictions
    drug_freqs = table(new_mat_civic$Drug)
    to_mark = drug_freqs > 1
    conflict_drugs = names(drug_freqs[to_mark])
    # Extra sanity check that duplicated drug is associated to conflicting predictions
    if (length(conflict_drugs)>0){
      checked_conflict_drugs = sapply(conflict_drugs, function(x){
        submat = subset(new_mat_civic, subset = Drug == x)
        logi_expr = (("CIVIC_RESISTANCE" %in% submat$Support) & ("CIVIC_SUPPORT" %in% submat$Support))
      })
      checked_conflict_drugs = names(checked_conflict_drugs[checked_conflict_drugs])
      # Sanity check
      stopifnot(conflict_drugs==checked_conflict_drugs)
      # Mark these drugs as associated to conflicting support (CIVIC_GENE_DEP)
      # final_mat_civic = subset(new_mat_civic, subset = Drug != checked_conflict_drugs)
      final_mat_civic = new_mat_civic
      final_mat_civic$Support[final_mat_civic$Drug %in% checked_conflict_drugs] = "CIVIC_GENE_DEP"
      final_mat_civic = unique(final_mat_civic)
    } else{
      final_mat_civic = new_mat_civic
    }
    final_mat_civic
  }, simplify=FALSE)
  
  # Extract CIVIC predicted drugs for each detected tumor cluster
  all_tumor_drugs_civic = sapply(names(all_civic_data_clean), function(id) {
    mat_civic = all_civic_data_clean[[id]]
    if (nrow(mat_civic)>0){
      mat_civic$phenograph_clusters = id  
    } else{
      mat_civic$phenograph_clusters = character()
    }
    mat_civic
  }, simplify=FALSE)
  all_tumor_data_civic = Reduce(rbind, all_tumor_drugs_civic)
 
  return(all_tumor_data_civic) 
}


reduce_allSamples_dataset = function(dataList){
  annotated = sapply(names(dataList), function(sampleName){
    df = dataList[[sampleName]]
    if (nrow(df)>0){
      df$Sample_name = sampleName
      stopifnot(unique(df$Sample_name)==sampleName)  
    } else{
      df$Sample_name = character()  
    }
    df
  }, USE.NAMES=T, simplify=F)
  
  # Sanity check that all individual dataframes have the same number of columns
  stopifnot(length(unique(sapply(annotated, ncol)))==1)
  # Join all dataframes into a single table
  reduced = Reduce(rbind,annotated)
  # Sanity check that dimensions match
  stopifnot(sum(sapply(dataList,nrow))==nrow(reduced))
  
  return(reduced)
}

  
process_pred_data = function(tumor_pred_data, malignant_ids, normal_ids){
  these_drugs = unique(tumor_pred_data$Drug)
  
  # Switch from long to wide format
  drug_mat = dcast(tumor_pred_data, unique_cluster_id ~ Drug, value.var="unique_cluster_id", length)
  
  # Add tumor clusters that do not have a drug prediction (if any) for correctly merging later on
  # Due to either: empty clinical file (only header), no drug predictions at all, or drugs not contained in input list
  missing_ids = malignant_ids[!malignant_ids %in% unique(tumor_pred_data$unique_cluster_id)]
  if (length(missing_ids)>0){
    missing_data = data.frame(matrix(c(missing_ids,rep(0,length(missing_ids)*length(these_drugs))), nrow=length(missing_ids)))
    colnames(missing_data) = c("unique_cluster_id",these_drugs)
    drug_mat = rbind(drug_mat,missing_data)
  }
  # Add normal clusters to the drug prediction dataframe for correctly merging later on
  norm_data = data.frame(matrix(c(normal_ids,rep(0,length(normal_ids)*length(these_drugs))), nrow=length(normal_ids)))
  # Sanity check for having no normal clusters (ie all malignant)
  if (nrow(norm_data)>0){
    colnames(norm_data) = c("unique_cluster_id",these_drugs)
    drug_mat = rbind(drug_mat,norm_data)
  }
  
  # Sanity check that all clusters (tumor and normal) are present in drug_mat
  stopifnot(all.equal(sort(drug_mat$unique_cluster_id),sort(cluster_ids)))
  # Sanity check that only 1/0 values occur in drug_mat
  stopifnot(sum(!names(table(as.matrix(subset(drug_mat, select = -unique_cluster_id)))) %in% c(0,1)) == 0)
  
  # TODO is all zeros possible?
  # TODO is all ones possible?
  
  return(drug_mat)
}


################################################################################
## initialize script
################################################################################

lby = c("optparse", "Seurat", "SingleCellExperiment", "reshape2", "ggplot2", "RColorBrewer")
resp = lapply(lby, require, character.only=T, warn.conflicts=F, quietly=T)
if(!all(unlist(resp))) stop("Could not load one or more packages")
rm(resp, lby)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

#options(error=function()traceback(2))
options(stringsAsFactors = FALSE)
fontsize = theme(axis.text=element_text(size=9), axis.title=element_text(size=11))
theme_set(theme_bw(12) + fontsize)

## Define prediction color scheme
# nnn,nnr*/nns*/nnd*,nyn*,ynn,nyr*/nys*/nyd*,yyn,ynr/yns/ynd,yyr/yys/yyd
pred_labels = c("nnn" = "Normal", "ynn" = "Tumor", "yyn" = "DGIDB prediction", 
                "yys" = "CIVIC support", "yyr" = "CIVIC resistance", "yyd" = "CIVIC gene-dependent",
                "yns" = "CIVIC support", "ynr" = "CIVIC resistance", "ynd" = "CIVIC gene-dependent")
pred_col = c("nnn" = "gray80", "ynn" = "gray55", "yyn" = "deepskyblue3", "yys" = "limegreen", "yyr" = "red", "yyd" = "orange", "yns" = "limegreen", "ynr" = "red", "ynd" = "orange")

## Define dummy labels and breaks for the legend
dummy_labels = c("Normal", "Tumor", "DGIDB prediction", "CIVIC support", "CIVIC resistance", "CIVIC gene-dependent")
dummy_breaks = c("nnn", "ynn", "yyn", "yys", "yyr", "yyd")

# convenience function for string concatenation
'%&%' = function(a,b) paste(a,b,sep="")

## FIXME: write updated documentation for parameters

# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to integrated sce object file containig cohort data (sce_celltypes_noatypical.RDS)."),
  make_option("--drugPredMap", type = "character", help = "File listing the path to each sample's CIVIC directory, which should contain a tab-separated file per detected *malignant* cell cluster in the given sample, each with clinical annotations for the corresponding cluster (with header). Assumes maps has no header"),
  make_option("--drugPredEnd", type = "character", help = "Complete file ending of desired clinical annotation files (eg. 'clinicalAnnotation.civic.txt'). The specified file ending must allow retrieval of malignant cluster IDs directly from the file name (eg. '[..].8.clinicalAnnotation.civic.txt')."),
  make_option("--name_DGIDB", type = "character", help = "Name of column containing DGIDB drug predictions in the clinical annotation files (eg. 'DGIDB-drugs(Score,Type)')."),
  make_option("--name_CIVIC", type = "character", help = "Name of column containing CIVIC drug predictions in the clinical annotation files (eg. 'CIViC_Drug_Support')."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--drugList", type = "character", help = "Path to tab separated file with clinically relevant priority drugs in the first column (with header)."),
  make_option("--combiList", type = "character", help = "Path to tab separated file with clinically relevant drug combinations in the first column (no header). Individual drugs are separated with '+', and they can correspond to combinations themselves (indicated with drugs listed between parenthesis and separated with '+')."),
  make_option("--civicDict", type = "character", help = "Path to tab separated file listing drug synonym pairs in CIVIC (with header): CIVIC name in first column and desired synonym in second column. The same drug can have more than 1 row corresponding to different synonyms."),
  make_option("--outTag", type = "character", help = "Cohort identifier. Attached to each output file name.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

path = opt$outputDirec %&% opt$outTag

print(opt$SCE)


################################################################################
## Main code starts here
################################################################################

## Load input data
sce_data = readRDS(opt$SCE)

# TODO: confirm with Anne this corresponds to single-sample phenograph cluster and not integrated ids (are they somewhere, anyway?)
sample_phenograph = cbind(sce_data[["sample_name"]],sce_data[["phenograph_clusters"]])
table(sample_phenograph)

sample_names = unique(sample_phenograph$sample_name)

## Retrieve input clinical annotation files corresponding to each integrated sample:
# Assumes no header and 2 columns: sample name (matching with integrated sce object) and path to CIVIC directory
drugPredMap = read.table(opt$drugPredMap, sep="\t", header=FALSE, stringsAsFactors=FALSE, quote="", check.names=FALSE)

# TODO include more extensive sanity check and error messages to user
stopifnot(identical(sort(drugPredMap[,1]),sort(sample_names)))

civicDirs = drugPredMap[,2]
names(civicDirs) = drugPredMap[,1]

# Sanity check that all provided CIVIC dirs exist
stopifnot(sum(sapply(civicDirs, dir.exists))==length(civicDirs))

# Generate right pattern for retrieving necessary input files
endPattern = gsub("\\.", "\\\\.", opt$drugPredEnd) %&% "$"

# FIXME
# All preprocessing should be done in a function, so that many clinical files can be integrated into 1 table per sample
# This way, it can be applied to many samples and return a list of nSamples integrated dataframes

# Retrieve drug prediction files from each directory
# For each sample, files correspond to all detected *malignant* cell clusters (1 file per cluster)
allFiles = sapply(civicDirs, function(x) list.files(path=x, pattern=endPattern, full.names=T) )
print(sapply(allFiles, length))

## Sanity check for finding some files in each directory
## NOTE: this should never happen as it is directly handled via snakemake
stopifnot(sum(sapply(allFiles, length)==0)==0)

## For each sample, retrieve their malignant cluster IDs detected as malignant (directly from file names)
allMaligIds = sapply(allFiles, function(x) retrieve_malig_ids(x,opt$drugPredEnd))

allFilesAnn = sapply(1:length(allFiles), function(i){
  files = allFiles[[i]]
  ids = allMaligIds[[i]]
  names(files) = ids
  files
})
names(allFilesAnn) = names(allFiles)

## Load clinical info for each sample (can be many files per sample)
allData = sapply(allFilesAnn, function(files) sapply(files, function(x) read.table(x, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", check.names=FALSE), simplify=FALSE))
print(str(allData, max.level=2))  

## Process DGIDB data separately for each sample
DGIDB = sapply(allData, function(all_clin) process_dgidb_single_sample(all_clin, opt$name_DGIDB), simplify=F, USE.NAMES=T)
## Join all single-sample dataframes into a single table
allDGIDB = reduce_allSamples_dataset(DGIDB)
print(str(allDGIDB))

## Substitute CIVIC drug names for their accepted synonyms when necessary
# Load dictionary table
input_civic_dict = read.table(opt$civicDict, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", check.names=FALSE)
# Sanity check that input dictionary has at least 2 columns, assumes: first = CIVIC drug names, second = drug synonym
stopifnot(ncol(input_civic_dict)>=2)
# Avoid missing a drug name because of case differences
civic_dict_tmp = data.frame(CIVIC = toupper(input_civic_dict[,1]), Synonym = toupper(input_civic_dict[,2]))
civic_dict = unique(civic_dict_tmp)
# Inform about duplicated entries in the dictionary
if (nrow(civic_dict_tmp) != nrow(civic_dict)){
  print("Detected duplicates in input CIVIC synonym dictionary.")
}
# Sanity check that CIVIC drug names are listed only once (ie. unique synonym per CIVIC drug)
n_civic_names = length(unique(civic_dict$CIVIC))
if (n_civic_names!=nrow(civic_dict)){
  print("Detected non-unique CIVIC names in first column of input synonyms dictionary. Need 1:1 correspondence to synonyms.")  
}
stopifnot(n_civic_names==nrow(civic_dict))


## Process CIVIC data separately for each sample
CIVIC = sapply(allData, function(all_clin) process_civic_single_sample(all_clin, opt$name_CIVIC, civic_dict), simplify=F, USE.NAMES=T)
## Join all single-sample dataframes into a single table
allCIVIC = reduce_allSamples_dataset(CIVIC)
print(str(allCIVIC))


## Load table of clinically-relevant drugs (individual list)
drugList = read.table(opt$drugList, sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="", check.names=FALSE)
# Avoid missing a drug because of case differences
tmp_drugs = data.frame(Drug = toupper(drugList[,1]))
# Remove duplicate drugs (if any)
inputDrugs = unique(tmp_drugs)

if (nrow(tmp_drugs) != nrow(inputDrugs)){
  print("Detected duplicates in input drug list.")
}

## For now, ignore drug combinations
########################################################################################################################
# ## Load table of drug combinations
# combiList = read.table(opt$combiList, sep="\t", header=FALSE, stringsAsFactors=FALSE, quote="", check.names=FALSE)
# # Avoid missing a drug because of case differences
# tmp_combs = data.frame(Drug = toupper(combiList[,1]))

# # Sort names of individual drugs in combination treatments alphabetically
# # Avoid mismatching treatments due to different order
# s.tmp_combs = data.frame(Drug = apply(tmp_combs, 1, function(x){
#   drugArr = strsplit(x,"\\+")[[1]]
#   paste(sort(drugArr), collapse="+")
# }))
# # Remove duplicate treatments after sorting (if any)
# inputCombs = unique(s.tmp_combs)

# if (nrow(tmp_combs) != nrow(inputCombs)){
#   print("Detected duplicates in input combiation treatments upon sorting individual drugs alphabetically.")
# }

# ## Distinguish between individual drugs and combination treatments when plotting
# inputDrugs$isCombi = F
# inputCombs$isCombi = T

allDrugs = inputDrugs   # for now, ignore combis
# allDrugs = rbind(inputDrugs,inputCombs)


## Drugs to be displayed in the final visualization are mostly determined by DGIDB prediction results
input_and_dgidb = Reduce(intersect, list(allDGIDB$Drug,allDrugs$Drug))
## But those with that have CIVIC info but no DGIDB info should also be included
## By intersecting to the input priority drugs, we ensure that CIVIC drugs overlapping with DGIDB priority drugs are also there
input_and_civic = Reduce(intersect, list(allCIVIC$Drug,allDrugs$Drug))
all_drugs = unique(c(input_and_dgidb,input_and_civic))

expected_overlap = intersect(input_and_dgidb,input_and_civic)       # both having dgidb and civic
expected_dgidb_only = setdiff(input_and_dgidb,input_and_civic)      # in dgidb and not in civic
expected_civic_only = setdiff(input_and_civic,input_and_dgidb)      # in civic and not in dgidb

## The final set of selected drugs should correspond to all DGIDB drugs + extra drugs from CIVIC
stopifnot(length(input_and_dgidb)+length(expected_civic_only)==length(all_drugs))
stopifnot(length(expected_overlap)+length(expected_dgidb_only)+length(expected_civic_only)==length(all_drugs))

ndrugs = length(all_drugs)

## Generate unique malignant cluster IDs by creating a new ID containing the originating sample name
## TODO: how to make sure that sample names in integrated sce and input map match? The key is the sample map!
malig_ids = as.character(unlist(sapply(names(allMaligIds), function(x) paste(x, allMaligIds[[x]], sep="_"))))

########################################################################################################################

## Sanity check for no drug predictions at all (or nothing being matched to our priority drug list)
## Save empty plot
if (ndrugs == 0){
  # if (ndrugs == 0 & ncombs == 0){
  p_empty <- ggplot() + theme_void()
  pdf(file=NULL)
  ggsave(path %&% ".drug_prediction_umap.png", plot=p_empty, dpi = 600, width = 50, height = 50, units = "cm")
  garbage <- dev.off()
  cat("No predictions found for clinically relevant treatments. Saving empty plot...")
  
  
  ## When at least 1 drug prediction was found
} else{
  
  ### Drug prediction painting on tSNE plots (phenograph)

  ## Retrieve data relevant for plotting (cell coordinates, cluster id)
  tomerge = as.data.frame(Embeddings(object = sce_data[["umap"]]))
  tomerge$barcodes = rownames(tomerge)
  nCells = nrow(tomerge)
  tomerge = merge(tomerge, as.data.frame(sce_data[[]]))
  # Sanity check that merging of both tables respects the number of cells
  stopifnot(nrow(tomerge)==nCells)
  
  ## Generate unique cluster IDs across samples by creating a new ID containing the source sample name
  tomerge$unique_cluster_id = paste(tomerge$sample_name, tomerge$phenograph_clusters, sep="_")
  print(table(tomerge$unique_cluster_id))
  
  ## FIXME: introduce permanent sanity check for tumor-only integrated cohorts?
  
  ## Retrieve tumor clusters (assign 'y'/'n' in variable 'tumor')
  cluster_ids = unique(tomerge$unique_cluster_id)
  normal_clusters = cluster_ids[!cluster_ids %in% (malig_ids)]
  cluster_dat = data.frame(unique_cluster_id = cluster_ids, tumor = ifelse(cluster_ids %in%  malig_ids, "y", "n"))
  # Sanity check that the same tumor clusters are found
  stopifnot(identical(sort(subset(cluster_dat, subset = tumor == "y")[,1]),sort(malig_ids)))
  
  tomerge = merge(tomerge, cluster_dat)
  # unique(subset(tomerge, select=c(unique_cluster_id,tumor)))
  stopifnot(length(setdiff(normal_clusters, cluster_dat$unique_cluster_id[cluster_dat$tumor=="n"]))==0)
  
  ## 1) Process DGIDB drug prediction data
  
  # Keep only predictions for drugs included in clinically-relevant input lists
  tumor_data_dgidb = merge(allDGIDB,allDrugs)
  # Retrieve all available DGIDB drugs
  dgidb_drugs = unique(tumor_data_dgidb$Drug)
  stopifnot(identical(sort(input_and_dgidb),sort(dgidb_drugs)))
  
  # Make sure that at least 1 priority drug from the input list had a DGIDB prediction
  if (nrow(tumor_data_dgidb)>0){
    # Generate unique cluster IDs across samples by creating a new ID containing the source sample name
    tumor_data_dgidb$unique_cluster_id = paste(tumor_data_dgidb$Sample_name, tumor_data_dgidb$phenograph_clusters, sep="_") 
    
    # Process drug prediction data across all clusters to a binary matrix of the form clusters (rows) x drugs (cols)
    drug_mat_dgidb = process_pred_data(subset(tumor_data_dgidb, select=c(Drug,unique_cluster_id)), malig_ids, normal_clusters)
    
    # Get data in required format for plotting (for each drug+tumor cluster, use "y"/"n" to indicate DGIDB predictions)
    td.plot = melt(drug_mat_dgidb, id.vars="unique_cluster_id")
    colnames(td.plot) = c("unique_cluster_id","Drug","DGIDB")
    td.plot$DGIDB[td.plot$DGIDB == 0] = "n"
    td.plot$DGIDB[td.plot$DGIDB == 1] = "y"
    stopifnot(sum(as.numeric(as.matrix(subset(drug_mat_dgidb, select=-unique_cluster_id))))==sum(td.plot$DGIDB == "y"))
    
  } else{
    ## Create empty dataframe so it can still be merged with the CIVIC one
    td.plot = data.frame(unique_cluster_id=character(), Drug=character(), DGIDB=character())
  }
  
  
  ## 2) Process CIVIC drug prediction data
  
  # Keep only predictions for drugs included in clinically-relevant input lists
  # We automatically ensure that CIVIC drugs overlapping with predicted DGIDB drugs are also contained (they are all in the input list)
  tumor_data_civic = merge(allCIVIC,allDrugs)
  # Retrieve all available CIVIC drugs
  civic_drugs = unique(tumor_data_civic$Drug)
  stopifnot(identical(sort(input_and_civic),sort(civic_drugs)))
  # Retrieve only-CIVIC drugs
  civic_only_drugs = setdiff(civic_drugs, dgidb_drugs)
  
  ## Sanity check of identical previous expected predictions and currently existing predictions
  dgidb_only_drugs = setdiff(dgidb_drugs, civic_drugs)
  both_drugs = intersect(dgidb_drugs, civic_drugs)
  stopifnot(identical(sort(both_drugs), sort(expected_overlap)))
  stopifnot(identical(sort(dgidb_only_drugs), sort(expected_dgidb_only)))
  stopifnot(identical(sort(civic_only_drugs), sort(expected_civic_only)))
  
  ## Report relevant drug predictions found across DGIDB and CIVIC
  print(paste0("Available DGIDB-only drugs: ", paste(dgidb_only_drugs, collapse=", ")))
  print(paste0("Available CIVIC-only drugs: ", paste(civic_only_drugs, collapse=", ")))
  print(paste0("Available DGIDB+CIVIC overlapping drugs: ", paste(both_drugs, collapse=", ")))
  
  # Make sure that all existing CIVIC support at this point correspond to either CIVIC_SUPPORT, CIVIC_RESISTANCE or CIVIC_GENE_DEP
  # This sanity check works even if there were no available or overlapping CIVIC drug predictions(ie 0 == 0)
  stopifnot(sum(tumor_data_civic$Support %in% c("CIVIC_SUPPORT", "CIVIC_RESISTANCE", "CIVIC_GENE_DEP")) == length(tumor_data_civic$Support))
  
  ## Further process CIVIC prediction data to include support information of listed DGIs (ie "sensitivity","resistance","gene-dependent")
  # Make sure that at least 1 priority drug from the input list had a CIVIC prediction
  if (nrow(tumor_data_civic)>0){
    # Generate unique cluster IDs across samples by creating a new ID containing the source sample name
    tumor_data_civic$unique_cluster_id = paste(tumor_data_civic$Sample_name, tumor_data_civic$phenograph_clusters, sep="_") 
    
    # Process drug prediction data across all clusters to a binary matrix of the form clusters (rows) x drugs (cols)
    drug_mat_civic = process_pred_data(subset(tumor_data_civic, select=c(Drug,unique_cluster_id)), malig_ids, normal_clusters)
    
    # Get data in required format for plotting (for each drug+tumor cluster, use "y"/"n" to indicate CIVIC predictions)
    tc.plot = melt(drug_mat_civic, id.vars="unique_cluster_id")
    colnames(tc.plot) = c("unique_cluster_id","Drug","CIVIC")
    tc.plot$CIVIC[tc.plot$CIVIC == 0] = "n"
    tc.plot$CIVIC[tc.plot$CIVIC == 1] = "y"
    stopifnot(sum(as.numeric(as.matrix(subset(drug_mat_civic, select=-unique_cluster_id))))==sum(tc.plot$CIVIC == "y"))
    
    # Get copy of CIVIC tumor prediction info for merging with CIVIC support info
    dummy_mat = tc.plot
    # Important to keep track of row indexes in original CIVIC dataframe
    dummy_mat$RowIndx = 1:nrow(dummy_mat)
    
    # Merge to find corresponding CIVIC rows in the original CIVIC support data (ie correct combinations of drug + tumor cluster)
    # Only look for actual CIVIC predictions ("y") as these are the only ones that will have support information
    stopifnot(nrow(dummy_mat[dummy_mat$CIVIC %in% "y",,drop=F])==nrow(tumor_data_civic))
    dummy_merge = merge(dummy_mat[dummy_mat$CIVIC %in% "y",,drop=F], tumor_data_civic, by=c("Drug","unique_cluster_id"))
    stopifnot(nrow(dummy_merge)==nrow(tumor_data_civic))
    
    # Retrieve row indexes in the original CIVIC dataframe
    # Distinguish between CIVIC support of "sensitivity", "resistance" or "gene-dependent"
    dummy_support = subset(dummy_merge, subset = Support == "CIVIC_SUPPORT")
    dummy_resist = subset(dummy_merge, subset = Support == "CIVIC_RESISTANCE")
    dummy_dep = subset(dummy_merge, subset = Support == "CIVIC_GENE_DEP")
    civic_rows_sens = unique(dummy_support$RowIndx)
    civic_rows_res = unique(dummy_resist$RowIndx)
    civic_rows_dep = unique(dummy_dep$RowIndx)
    
    ## Annotate CIVIC support for drugs + tumor clusters found in CIVIC (if any)
    ## Use 'r', 's' or 'd' depending on whether a given drug+cluster had CIVIC_RESISTANCE, CIVIC_SUPPORT or CIVIC_GENE_DEP, respectively
    tc.plot$CIVIC[civic_rows_sens] = "s"
    tc.plot$CIVIC[civic_rows_res] = "r"
    tc.plot$CIVIC[civic_rows_dep] = "d"
    
  } else{
    ## Create empty dataframe so it can still be merged with the DGIDB one
    tc.plot = data.frame(unique_cluster_id=character(), Drug=character(), CIVIC=character())
  }
  
  ## 3) Merge all drug prediction data (DGIDB and CIVIC)
  ## NOTE: here it is critical to use "all=T" to ensure that CIVIC-only and DGIDB-only predictions are also kept
  t.plot = merge(td.plot, tc.plot, by=c("Drug","unique_cluster_id"), all=T)
  
  # Sanity check that at least something could be merged
  stopifnot(nrow(t.plot)>0)
  
  ## Make use of NAs to infer what predictions come from each database (DGIDB, CIVIC or both)
  # Sanity check that there are no drug predictions having NA for both DGIDB and CIVIC (not possible)
  stopifnot(sum(is.na(t.plot$CIVIC) & is.na(t.plot$DGIDB))==0)
  table(t.plot$DGIDB, t.plot$CIVIC, useNA="always")
  
  ## Predictions that have both DGIDB and CIVIC
  both_preds = nrow(t.plot[!is.na(t.plot$CIVIC) & !is.na(t.plot$DGIDB),,drop=F])
  
  ## Predictions that only have DGIDB info
  dgidb_preds = nrow(t.plot[is.na(t.plot$CIVIC) & !is.na(t.plot$DGIDB),,drop=F])
  # Substitute corresponding NAs
  t.plot$CIVIC[is.na(t.plot$CIVIC) & !is.na(t.plot$DGIDB)] = "n"
  
  ## Predictions that only have CIVIC info
  civic_preds = nrow(t.plot[is.na(t.plot$DGIDB) & !is.na(t.plot$CIVIC),,drop=F])
  # Substitute corresponding NAs
  t.plot$DGIDB[is.na(t.plot$DGIDB) & !is.na(t.plot$CIVIC)] = "n"
  
  # Sanity check that no predictions were lost along the processing steps
  # Number of current predictions should be equal to number of (overlapping + dgidb-only + civic-only) predictions
  stopifnot(both_preds+dgidb_preds+civic_preds==nrow(t.plot))
  table(t.plot$DGIDB, t.plot$CIVIC, useNA="always")
  # table(t.plot$DGIDB %&% t.plot$CIVIC, useNA="always")
  
  print("Drugs + clusters with DGIDB drug predictions (DGIDB=y):")
  print(subset(t.plot, subset = DGIDB != "n"))
  print("Drugs + clusters with CIVIC drug predictions (CIVIC=s/r/d):")
  print(subset(t.plot, subset = CIVIC != "n"))
  
  ## Merge all data necessary for plotting
  dfPlot = merge(t.plot, tomerge)
  
  # Make sure that all clusters (tumor and normal) are present for plotting
  stopifnot(length(unique(dfPlot$unique_cluster_id)) == length(unique(paste(sce_data[["sample_name"]][,1], sce_data[["phenograph_clusters"]][,1], sep="_"))))
  # Make sure that all drugs are present for plotting
  stopifnot(length(unique(dfPlot$Drug))==ndrugs && length(intersect(unique(dfPlot$Drug), all_drugs))==ndrugs)
  
  # Generate drug prediction values for plotting
  # nnn,nnr*/nns*/nnd*,nyn*,ynn,nyr*/nys*/nyd*,yyn,ynr/yns/ynd,yyr/yys/yyd
  dfPlot$Prediction = dfPlot$tumor %&% dfPlot$DGIDB %&% dfPlot$CIVIC
  
  # FIXME
  # subPlot = subset(dfPlot, select=c(unique_cluster_id,Drug,DGIDB,CIVIC,barcodes,UMAP_1,UMAP_2,phenograph_clusters,sample_name,tumor,Prediction))
  # identical(subPlot$unique_cluster_id, paste(subPlot$sample_name,subPlot$phenograph_clusters,sep="_"))
  # table(unique(subset(subPlot, select=c(unique_cluster_id,tumor)))[,2])
  # print(table(dfPlot$Prediction, useNA="always"))
  # subDrug = subset(dfPlot, select=c(Drug,Prediction))
  # table(unique(subDrug))
  
  # Make sure that only tumor clusters have a drug prediction associated (no "nyr"/"nys"/"nyd", "nyn" or "nnr"/"nns"/"nnd")
  # NOTE: CIVIC predictions can occur even when there is no DGIDB prediction available ("ynr"/"yns"/"ynd")
  stopifnot(sum(unique(dfPlot$Prediction) %in% c("nyr","nys","nyd","nyn","nnr","nns","nnd"))==0)
  dfPlot$Prediction = factor(dfPlot$Prediction, levels=names(pred_labels))
  
  # Sanity check that no NAs (would occur after factoring when non-allowed categories were present)
  # This should not happen anymore because of above stopifnot() clause
  stopifnot(sum(is.na(dfPlot$Prediction))==0)
  print("Distribution of drug prediction categories (tumor?/dgidb?/civic?):")
  print(table(dfPlot$Prediction, useNA="always"))
  
  # NOTE: for cohort plots, DGIDB-only drug predictions are not shown, only CIVIC predictions
  # Remove DGIDB-only predictions (ie. "yyn") and assign them as Tumor (ie. "ynn")
  oldPrediction = dfPlot$Prediction
  dfPlot$Prediction[dfPlot$Prediction == "yyn"] = "ynn"
  print("Distribution after removing DGIDB-only predictions (tumor?/dgidb?/civic?):")
  print(table(dfPlot$Prediction, useNA="always"))
  # Sanity check that the numbers match
  stopifnot(sum(table(dfPlot$Prediction, useNA="always"))==sum(table(oldPrediction, useNA="always")))
  
  ## Remove drugs that only had DIGDB-only drug predictions
  newDf = subset(dfPlot, select=c(Drug, Prediction), subset = !Prediction %in% c("nnn","ynn","yyn"))
  # table(newDf[,2])
  # table(newDf)
  display_drugs = as.character(unique(newDf$Drug))
  
  dfPlot = subset(dfPlot, subset = Drug %in% display_drugs)
  
  # Create new drug labels for facet titles in case combination treatments exist (split at "+")
  dfPlot$Drug_wrap = gsub("\\+", "+\n", dfPlot$Drug)
  # Force alphabetical sorting of drug names
  dfPlot$Drug_wrap = factor(dfPlot$Drug_wrap, levels=sort(unique(dfPlot$Drug_wrap)))
  
  ## For having more or less squared facets, use this rule of thumb
  nPanels = sqrt(length(display_drugs))
  str(dfPlot)

  # What samples had CIVIC infor for drug SUNITINIB?
  # print(unique(subset(dfPlot, select=sample_name, subset = Drug == "SUNITINIB" & !Prediction %in% c("nnn","ynn","yyn"))))
  
  ## Remove DGIDB-only category from legend
  ## TODO
  
  ## Add transparency only to those points that have no prediction
  dfPlot$alpha = 1
  dfPlot$alpha[dfPlot$Prediction %in% c("nnn","ynn","yyn")] = 0.5 # "yyn" added for consistency; should not occur after filtering
  
  pp = ggplot(dfPlot, aes(x = UMAP_1, y = UMAP_2) ) +
    geom_point(aes(color=Prediction, alpha=alpha), size=1) + xlab("UMAP_1") + ylab("UMAP_2") + 
    scale_colour_manual(values=pred_col, labels=dummy_labels, breaks=dummy_breaks, drop=FALSE) +
    facet_wrap(~Drug_wrap, ncol=nPanels) +
    theme(legend.text=element_text(size=17),
          legend.title=element_text(size=17)) +
    guides(colour = guide_legend(override.aes = list(size=5)), alpha=FALSE)
  ggsave(path %&% ".drug_prediction_umap.png", pp,
         dpi = 600, width = 50, height = 50, units = "cm")
}
