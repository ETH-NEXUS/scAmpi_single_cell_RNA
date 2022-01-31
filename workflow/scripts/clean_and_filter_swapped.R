suppressPackageStartupMessages(library(DropletUtils))

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

args <- commandArgs(trailingOnly = TRUE)

## assumes first parameter = output RData object for cleaned (and unfiltered) data
novaFile <- args[1]
## assumes second parameter = root output folder, also for saving statistic files (cell calling comparison, summary tables, etc.)
outTag_nova <- args[2]
## assumes third parameter = FDR threshold for cell calling
fdr_thres <- args[3]
## assumes fourth parameter = library size threshold for cell calling
libSize_thres <- args[4]
## assumes fifth to nth parameter = input files for clean and filtered data (will derive output filenames from there)
allFiles_nova <- args[5:length(args)]


### Check input threshold parameters
fdr_thres <- as.double(fdr_thres)
libSize_thres <- as.double(libSize_thres)
if (is.na(fdr_thres) | is.na(libSize_thres)){
    stop("\nInvalid input parameters introduced. Please use valid threshold values!\n")
}

cat("\n")
cat(paste0("Input FDR threshold: ", fdr_thres, "\n"))
cat(paste0("Input library size threshold: ", libSize_thres, "\n"))


### 1) Perform swapping removal on NovaSeq data
cat("\nInput files:\n")
cat(paste(" ",allFiles_nova,collapse="\n"), "\n\n")

print("Performing swapping removal on input data...")
time_nova <- system.time(cleaned_nova <- swappedDrops(allFiles_nova, get.diagnostics = TRUE, get.swapped = TRUE))
print("Run time:")
print(time_nova)

sampleNames_nova <- sapply(allFiles_nova, function(x) {
  firstSplit <- strsplit(x, ".molecule_info.h5")[[1]][1]
  secondSplit <- strsplit(firstSplit, "\\/")[[1]]
  secondSplit[length(secondSplit)]
})
outFiles_nova <- sapply(sampleNames_nova, function(x) paste0(outTag_nova, "/", x, "/"))

## Important to assign sample names to their corresponding matrix before saving
names(cleaned_nova$cleaned) <- names(cleaned_nova$swapped) <- colnames(cleaned_nova$diagnostics) <- sampleNames_nova
## Save .RData object containing the results from applying the swapping removal
print("Saving clean unfiltered data as .RData object...")
time_save_nova <- system.time(save(cleaned_nova, file=novaFile))
print("Save time:")
print(time_save_nova)
cat("\n\n")


### 2) Call cells on clean (ie. unswapped) data

print("Cell calling without swapping removal...")
cat("\n")
tab_cells_nova <- data.frame(Samples = sampleNames_nova,
                             Cells_called_unswapped = rep(NA, length(sampleNames_nova)))

for (i in seq_along(cleaned_nova$cleaned)){
    print(paste0("Processing sample ", sampleNames_nova[i], "..."))
    outFile_nova <- outFiles_nova[i]
    cleanMat_nova <- cleaned_nova$cleaned[[i]]
    print("Dataset dimensions:")
    print(dim(cleanMat_nova))

    ## Use custom function for calling cells (as CellRanger's version will not work over filtered data)
    print("Applying cell filtering with EmptyDrops()...")
    ## The p-values are calculated by permutation testing, hence the need to set a seed.
    set.seed(100)
    ## Perform test for empty cells
    time_cell_nova <- system.time(test_nova <- emptyDrops(cleanMat_nova, niters = 30000, ignore = NULL))
    print("Run time:")
    print(time_cell_nova)
    ## And apply filtering to resulting FDR and library size
#     isCell_nova <- test_nova$FDR <= 0.01 & test_nova$Total >= 1000
    isCell_nova <- test_nova$FDR <= fdr_thres & test_nova$Total >= libSize_thres
    print(table(isCell_nova, useNA="always"))
    isCell_nova[is.na(isCell_nova)] <- FALSE
    print(table(isCell_nova, useNA="always"))
    print("Diagnostics table:")
    print(table(Limited=test_nova$Limited, Significant=isCell_nova))
    ## Calculate final number of called cells
    ncells_nova <- sum(isCell_nova, na.rm = TRUE)
    print(paste0("Total no. cells: ", ncells_nova))

    ## Call cells and filter UMI count matrix accordingly
    print("Filtering UMI count matrix...")
    filtered_mat_nova <- cleanMat_nova[,isCell_nova]

    ## Save cleaned and filtered UMI count matrices
    print(paste0("Saving processed files to ", outFile_nova, "..."))
    write10xCounts(filtered_mat_nova, path=outFile_nova, version="3", type="sparse", overwrite=TRUE)

    ## Do diagnostics plot
    temp_cellFile_nova <- paste0(outFile_nova, "cell_call_diagnostics.pdf")
    print(paste0("Saving corresponding diagnostics plot in file ", temp_cellFile_nova,"..."))
    pdf(file=temp_cellFile_nova)
    plot(test_nova$Total, -test_nova$LogProb, col=ifelse(isCell_nova, "red", "black"), xlab="Total UMI count", ylab="-Log Probability")
    dev.off()
    cat("\n")

    ## Store results for each sample
    tab_cells_nova$Cells_called_unswapped[i] <- ncells_nova
}
write.table(tab_cells_nova, file=paste0(outTag_nova, "/summary_cells.txt"), sep="\t", quote=F, row.names=F, col.names=T)


### 3) Calculate statistics

cat("\n\n")
print("Calculating summary statistics...")
cat("\n")

## 3a) Swapped fraction and excluded percent over whole dataset

## Use diagnostics matrix: total reads for molecules (rows) x samples (cols), where molecules correspond to unique tuples UMI+barcode+gene.

## Calculate the overall percent of removed molecules
## 1 - (total no. molecules (after) / total no. molecules before)
excluded_nova <- format( ( 1 - sum(sapply(cleaned_nova$cleaned, sum))/nrow(cleaned_nova$diagnostics) ) * 100, digits = 3)

## Get total number of reads across all molecules and samples
totalReads <- sum(cleaned_nova$diagnostics)
## Get total number of molecules
totalMols <- nrow(cleaned_nova$diagnostics)

## Get swapped fraction, ie. fraction of molecules observed in more than 1 sample
## For each molecule, get no. samples with reads > 0
sampleSum <- rowSums(cleaned_nova$diagnostics>0)

## Which molecules have no.samples>1? (ie. swapping) As each distinct tuple should only correspond to one sample (hence, no.reads>0 should only occur for one sample)
swapped_tuples <- sampleSum>1
## Get fraction of molecules (identical UMI+cell barcode+aligned gene) observed in more than 1 sample
swap_fraction <- (sum(swapped_tuples)/totalMols) * 100

tab_nova <- data.frame(Reads = format(totalReads, big.mark = ","),
                      Molecules = format(totalMols, big.mark = ","),
                      Swapped_molecules = format(sum(swapped_tuples), big.mark = ","),
                      Swapped_percent = paste0(format(swap_fraction, digits = 3), " %"),
                      Removed_percent = paste0(excluded_nova, " %"))

tab_nova_num <- data.frame(Reads = totalReads,
                          Molecules = totalMols,
                          Swapped_molecules = sum(swapped_tuples),
                          Swapped_percent = swap_fraction,
                          Removed_percent = excluded_nova)

print("Summary statistics over whole dataset:")
print(tab_nova)
cat("\n\n")
summary_all_nova <- paste0(outTag_nova, "/summaryStats.txt")
write.table(tab_nova_num, file=summary_all_nova, sep="\t", quote=F, row.names=F, col.names=T)


## 3b) Swapped fraction across individual samples

tab_swap_samples <- data.frame(Samples = sampleNames_nova,
                               Reads = rep(NA, length(sampleNames_nova)),
                               Molecules = rep(NA, length(sampleNames_nova)),
                               Swapped_molecules = rep(NA, length(sampleNames_nova)),
                               Swapped_percent = rep(NA, length(sampleNames_nova)),
                               Removed_percent = rep(NA, length(sampleNames_nova)))

tab_swap_samples_num <- data.frame(Samples = sampleNames_nova,
                               Reads = rep(NA, length(sampleNames_nova)),
                               Molecules = rep(NA, length(sampleNames_nova)),
                               Swapped_molecules = rep(NA, length(sampleNames_nova)),
                               Swapped_percent = rep(NA, length(sampleNames_nova)),
                               Removed_percent = rep(NA, length(sampleNames_nova)))

for (i in seq_along(sampleNames_nova)){
  sample <- sampleNames_nova[i]
  cat("\n")
  print(paste0("Calculating swapped fraction for NovaSeq sample ", sample, "..."))
  ## Subset the column in the diagnostics matrix corresponding to the corresponding sample
  k <- which(sampleNames_nova %in% sample)
  subSample <- cleaned_nova$diagnostics[,k]
  
  ## Get UMI count matrices before and after swapping removal for current sample
  postMat <- cleaned_nova$cleaned[[i]]
  priorMat <- cleaned_nova$cleaned[[i]] + cleaned_nova$swapped[[i]]
  ## Divide total no. of UMIs (ie. distinct molecules) after removal by total no. of UMIs before removal
  ## This gives the fraction of remaining molecules, then get percent of excluded molecules for the current sample
  excluded_sample <- (1 - (sum(postMat)/sum(priorMat))) * 100
  tab_swap_samples$Removed_percent[i] <- paste0(format(excluded_sample, digits=3), " %")
  tab_swap_samples_num$Removed_percent[i] <- excluded_sample
  
  ## Get total number of reads for the current sample
  totalReads_sample <- sum(subSample)
  tab_swap_samples$Reads[i] <- format(totalReads_sample, big.mark = ",")
  tab_swap_samples_num$Reads[i] <- totalReads_sample
  
  ## Get swapped fraction, ie. fraction of molecules observed in more than 1 sample
  ## Get molecules present in current sample (ie. with reads > 0)
  molsInSample <- subSample>0
  ## Intersect molecules present in current sample with swapped molecules (ie. no. samples with reads>0 is > 1)
  swappedInSample <- sum(molsInSample & swapped_tuples)
  molsInSample <- sum(molsInSample)
  tab_swap_samples$Molecules[i] <- format(molsInSample, big.mark = ",")
  tab_swap_samples$Swapped_molecules[i] <- format(swappedInSample, big.mark = ",")
  tab_swap_samples_num$Molecules[i] <- molsInSample
  tab_swap_samples_num$Swapped_molecules[i] <- swappedInSample
  
  ## Get fraction of swapped molecules (identical UMI+cell barcode+aligned gene) observed in the current sample
  swap_fraction_sample <- (swappedInSample/molsInSample) * 100
  tab_swap_samples$Swapped_percent[i] <- paste0(format(swap_fraction_sample, digits=3), " %")
  tab_swap_samples_num$Swapped_percent[i] <- swap_fraction_sample
  cat("\n")
}

print("Summary statistics across individual samples:")
print(tab_swap_samples)
cat("\n\n")
summary_samples_nova <- paste0(outTag_nova, "/summaryStats_samples.txt")
write.table(tab_swap_samples_num, file=summary_samples_nova, sep="\t", quote=F, row.names=F, col.names=T)


### Check whether both 3pr and 5pr libraries are present

samples_5pr <- grepl("_5pr_", sampleNames_nova)
samples_3pr <- grepl("_3pr_", sampleNames_nova)

if (sum(samples_5pr) < 1 | sum(samples_3pr) < 1){
    cat("\n\nOutput files and directories:\n")
    cat(paste(" ",c(novaFile, summary_all_nova, summary_samples_nova, outFiles_nova),collapse="\n"), "\n\n")

} else{

    ## 3c) Swapped fraction across 3pr and 5pr libraries (if any)

    library_types <- c("3pr", "5pr")
    tab_swap_libs <- data.frame(Libraries = library_types,
                                Reads = rep(NA, length(library_types)),
                                Molecules = rep(NA, length(library_types)),
                                Swapped_molecules_overall = rep(NA, length(library_types)),
                                Swapped_percent_overall = rep(NA, length(library_types)),
                                Removed_percent_overall = rep(NA, length(library_types)))

    tab_swap_libs_num <- data.frame(Libraries = library_types,
                                Reads = rep(NA, length(library_types)),
                                Molecules = rep(NA, length(library_types)),
                                Swapped_molecules_overall = rep(NA, length(library_types)),
                                Swapped_percent_overall = rep(NA, length(library_types)),
                                Removed_percent_overall = rep(NA, length(library_types)))


    for (i in seq_along(library_types)){
      lib <- library_types[i]
      cat("\n")
      print(paste0("Calculating overall swapped fraction for ", lib, " samples..."))
      ## Subset the column in the diagnostics matrix corresponding to the corresponding sample
      lib_samples <- grepl(lib, sampleNames_nova)
      subLib <- cleaned_nova$diagnostics[,lib_samples,drop=F]

      ## Get UMI count matrices before and after swapping removal for current sample
      postLib <- sum(sapply(cleaned_nova$cleaned[lib_samples],sum))
      priorLib <- sum(sapply(which(lib_samples), function(i) sum(cleaned_nova$cleaned[[i]] + cleaned_nova$swapped[[i]])))
      ## Divide total no. of UMIs (ie. distinct molecules) after removal by total no. of UMIs before removal
      ## This gives the fraction of remaining molecules, then get percent of excluded molecules for the current sample
      excluded_lib <- (1 - (postLib/priorLib)) * 100
      tab_swap_libs$Removed_percent_overall[i] <- paste0(format(excluded_lib, digits=3), " %")
      tab_swap_libs_num$Removed_percent_overall[i] <- excluded_lib

      ## Get total number of reads for the current sample
      totalReads_lib <- sum(subLib)
      tab_swap_libs$Reads[i] <- format(totalReads_lib, big.mark = ",")
      tab_swap_libs_num$Reads[i] <- totalReads_lib

      ## Get swapped fraction, ie. fraction of molecules observed in more than 1 sample
      ## Get molecules present in current library (ie. with reads > 0)
      molsInLib <- rowSums(subLib>0)>0
      ## Intersect molecules present in current sample with swapped molecules (ie. no. samples with reads>0 is > 1)
      swappedInLib <- sum(molsInLib & swapped_tuples)
      molsInLib <- sum(molsInLib)
      tab_swap_libs$Molecules[i] <- format(c(molsInLib), big.mark = ",")
      tab_swap_libs$Swapped_molecules_overall[i] <- format(c(swappedInLib), big.mark = ",")
      tab_swap_libs_num$Molecules[i] <- molsInLib
      tab_swap_libs_num$Swapped_molecules_overall[i] <- swappedInLib

      ## Get fraction of swapped molecules (identical UMI+cell barcode+aligned gene) observed in the current library
      swap_fraction_lib <- (swappedInLib/molsInLib) * 100
      tab_swap_libs$Swapped_percent_overall[i] <- paste0(format(swap_fraction_lib,digits=3), " %")
      tab_swap_libs_num$Swapped_percent_overall[i] <- swap_fraction_lib
      cat("\n")
    }

    print("Library summary statistics:")
    print(tab_swap_libs)
    cat("\n\n")
    summary_libs_nova <- paste0(outTag_nova, "/summaryStats_libraries.txt")
    write.table(tab_swap_libs_num, file=summary_libs_nova, sep="\t", quote=F, row.names=F, col.names=T)


    ## 3d) Get statistics for particular swapping events

    swap_categs <- c("3pr_to_3pr", "5pr_to_5pr", "mixed")
    tab_swap_split <- data.frame(Libraries = swap_categs,
                                 Total_swapped_molecules = rep(NA, length(swap_categs)),
                                 Swapped_molecules = rep(NA, length(swap_categs)),
                                 Percent_over_all_swapped = rep(NA, length(swap_categs)))

    tab_swap_split_num <- data.frame(Libraries = swap_categs,
                                     Total_swapped_molecules = rep(NA, length(swap_categs)),
                                     Swapped_molecules = rep(NA, length(swap_categs)),
                                     Percent_over_all_swapped = rep(NA, length(swap_categs)))


    ### Which swapped molecules correspond to 3pr-only swapping (3->3, internal), 5pr-only swapping (5->5, internal) or mixed swapping (3->5 or 5->3 or mixed, external)?
    ## Subset diagnostics matrix to only contain swapped molecules
    subDiagn <- cleaned_nova$diagnostics[swapped_tuples,]

    ## Which swapped molecules are present in which samples?
    which_samples <- subDiagn>0

    ## Sanity check
    # table(rowSums(which_samples)>1) # all should be present in >1 samples
    # table(rowSums(cleaned_nova$diagnostics[!swapped_tuples,]>0)>1) # none should be present in >1 samples

    ## For each swapped molecule, sum up the number of 3pr and 5pr samples that contain it
    ## Subset 3pr samples and sum
    n_3pr <- rowSums(which_samples[,samples_3pr,drop=F])
    ## Subset 5pr samples and sum
    n_5pr <- rowSums(which_samples[,samples_5pr,drop=F])
    ## Get swapped molecules for each category
    molecules_3pr <- sum(n_3pr>0 & n_5pr<1)   # present in 3pr samples and no 5pr samples
    molecules_5pr <- sum(n_3pr<1 & n_5pr>0)   # present in 5pr samples and no 3pr samples
    mixed_molecules <- sum(n_3pr>0 & n_5pr>0) # present in 3pr samples and 5pr samples
    total_swapped <- sum(swapped_tuples)

    ## Swapped fraction for each category
    swapped_fraction_3pr <- (molecules_3pr/total_swapped) * 100
    swapped_fraction_5pr <- (molecules_5pr/total_swapped) * 100
    swapped_fraction_mixed <- (mixed_molecules/total_swapped) * 100

    ## Report in table
    tab_swap_split[1,2:4] <- c(format(total_swapped, big.mark=","), format(molecules_3pr, big.mark=","), paste0(format(swapped_fraction_3pr, digits=3), " %"))
    tab_swap_split[2,2:4] <- c(format(total_swapped, big.mark=","), format(molecules_5pr, big.mark=","), paste0(format(swapped_fraction_5pr, digits=3), " %"))
    tab_swap_split[3,2:4] <- c(format(total_swapped, big.mark=","), format(mixed_molecules, big.mark=","), paste0(format(swapped_fraction_mixed, digits=3), " %"))

    tab_swap_split_num[1,2:4] <- c(total_swapped, molecules_3pr, swapped_fraction_3pr)
    tab_swap_split_num[2,2:4] <- c(total_swapped, molecules_5pr, swapped_fraction_5pr)
    tab_swap_split_num[3,2:4] <- c(total_swapped, mixed_molecules, swapped_fraction_mixed)

    print("Library summary statistics for individual swapping events:")
    print(tab_swap_split)
    cat("\n\n")
    summary_events_nova <- paste0(outTag_nova, "/summaryStats_events.txt")
    write.table(tab_swap_split_num, file=summary_events_nova, sep="\t", quote=F, row.names=F, col.names=T)

    cat("\n\nOutput files and directories:\n")
    cat(paste(" ",c(novaFile, summary_all_nova, summary_samples_nova, summary_libs_nova, summary_events_nova, outFiles_nova),collapse="\n"), "\n\n")
}

