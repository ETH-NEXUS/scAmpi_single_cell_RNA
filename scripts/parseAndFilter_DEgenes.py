#!/usr/bin/env python

'''
Parse a table with differentially expressed genes.
Franziska Singer, March 2018
'''

import sys
import os
import re


if (len(sys.argv) <= 1) or ("-h" in sys.argv[1]):
    print("Parse and filter differentially expressed genes.")
    print("Usage: python parseAndFilter_DEgenes.py [inputTable] [outfile] [pvalue-threshold] [columnName_geneNames] [columnName_pvalues] [columnName_diff] [columnName_testStatistic] [colName_nonmalMax] [colName_nonmalMin] [colName_maligMean] [diffThreshold] [diff2secondThreshold] [optionalColumn1,optionalColumn2...]")
    sys.exit(1)


inputTable = sys.argv[1]
outName = sys.argv[2]
pvalThreshold = float(sys.argv[3])
colName_gene = sys.argv[4]
colName_pval = sys.argv[5]
colName_diff = sys.argv[6]
colName_testStat = sys.argv[7]
colName_nonmalMax = sys.argv[8]
colName_nonmalMin = sys.argv[9]
colName_maligMean= sys.argv[10]
diffThreshold = float(sys.argv[11])
diff2secondThreshold = float(sys.argv[12])
optColumns = sys.argv[13]

print("Input:\nPvalue threshold: %s\nColumn name genes: %s\nColumn name pvalue: %s\nColumn name estimate of difference: %s\nColumn name test statistic: %s\nColumn name maximum mean of non-malignant cell type: %s\nColumn name minimum mean of non-malignant cell type: %s\nColumn name mean malignant cluster: %s\nDifference between means threshold: %s\nThreshold for difference between means: %s\nOptional columns: %s\n" %(str(pvalThreshold),colName_gene,colName_pval,colName_diff,colName_testStat,colName_nonmalMax,colName_nonmalMin,colName_maligMean,diffThreshold,diff2secondThreshold,optColumns))

infile = open(inputTable,'r')
outfile = open(outName,'w')

geneCount = 0
geneFilterIn = 0

# check if any optional columns have been provided
optSplit = optColumns.strip().split(",")
# opt will be empty when no optional columns have been provided
opt = []    # list of names of the provided optional columns
if optSplit == ['']:
    opt = []
else:
    for tmpOpt in optSplit:
        # Skip potential empty strings; eg. if "pct.1," was provided
        if not tmpOpt.strip():
            continue
        if tmpOpt.strip() not in opt:
            opt.append(tmpOpt.strip())
        else:
            print("Warning! Optional column '%s' was provided twice as an input." %(tmpOpt.strip()))
            continue


# first line is used to get column indices
index_gene = -1
index_pval = -1
index_diff = -1
index_testStat = -1
index_nonmalMax = -1
index_nonmalMin = -1
index_maligMean = -1
# if any optional columns were provided, search for their indices as well
# order of indices in index_opt follows a 1-1 correspondance with order of optional columns as provided in input
if opt:
    index_opt = [-1] * len(opt)

firstLineSplit = infile.readline().strip().split("\t")
for pos in range(0,len(firstLineSplit)):
    if colName_gene == firstLineSplit[pos]:
        index_gene = pos
        print("gene found")
    if colName_pval == firstLineSplit[pos]:
        index_pval = pos
        print("pval found")
    if colName_diff == firstLineSplit[pos]:
        index_diff = pos
        print("foldChange found")
    if colName_testStat == firstLineSplit[pos]:
        index_testStat = pos
        print("testStat found")
    if colName_nonmalMax == firstLineSplit[pos]:
        index_nonmalMax = pos
        print("nonmalMax found")
    if colName_nonmalMin == firstLineSplit[pos]:
        index_nonmalMin = pos
        print("nonmalMin found")
    if colName_maligMean == firstLineSplit[pos]:
        index_maligMean = pos
        print("maligMean found")
        # when optional columns were provided, check and retrieve their indices as well
    if opt:
        # need to iterate through all provided columns
        for i,colName_opt in enumerate(opt):
            if colName_opt == firstLineSplit[pos]:
                print("i:")
                print(i)
                index_opt[i] = pos
                print("%s found" %(colName_opt))

if index_gene == -1 or index_pval == -1 or index_diff == -1 or index_testStat == 1 or index_nonmalMax == -1 or index_nonmalMin == -1 or index_maligMean == -1:
    print("Error! Could not find all necessary columnNames in header:\n%s" %(firstLineSplit))
    sys.exit(1)
# when optional columns were provided, check that all their corresponding indices were found as well
if opt:
    if -1 in index_opt:
        print("Error! Could not find all necessary optional columnNames in header:\n%s" %(firstLineSplit))
        sys.exit(1)

# Now that we have the indices, write header to outfile and then parse rest of the table
# distinguish between cases where optional columns were provided or not
if opt:
    outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(colName_gene,colName_diff,colName_pval,colName_testStat,colName_maligMean,colName_nonmalMax,colName_nonmalMin,"\t".join(opt)))
else:
    outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(colName_gene,colName_diff,colName_pval,colName_testStat,colName_maligMean,colName_nonmalMax,colName_nonmalMin))

filtered_NA = 0
for line in infile:
    lineSplit = line.strip().split('\t')
    geneCount += 1
    gene = lineSplit[index_gene]
    if "NA" in lineSplit[index_pval]:
        filtered_NA += 1
        continue
    pval = float(lineSplit[index_pval])
    if "NA" in lineSplit[index_diff]:
        filtered_NA += 1
        continue
    diffEst = float(lineSplit[index_diff])
    if "NA" in lineSplit[index_testStat]:
        filtered_NA += 1
        continue
    testStat = float(lineSplit[index_testStat])
    malig_mean = float(lineSplit[index_maligMean])
    nonmal_max = float(lineSplit[index_nonmalMax])
    nonmal_min = float(lineSplit[index_nonmalMin])
    if opt:
        # optVals will contain all the values from the current line for the provided optional columns
        optVals = []
        # Skip current line if at least 1 of the columns has value 'NA'
        hasNA = False
        # order of indices in index_opt follows a 1-1 correspondance with order of optional columns in output
        for optIndx in index_opt:
            if "NA" in lineSplit[optIndx]:
                hasNA = True
                continue
            optVal = float(lineSplit[optIndx])
            optVals.append(str(optVal))
        if hasNA:
            filtered_NA += 1
            continue
    # perform the actual filtering
    # check if p-value is below threshold
    if pval <= pvalThreshold:
        # check if estimated difference (logFC) is above threshold
        if abs(diffEst) >= diffThreshold:
            # check if either mean of malignant cluster is above nonmal_max
            if malig_mean > nonmal_max:
                diff2second = malig_mean - nonmal_max
                # and if the difference between the two is above the threshold
                if diff2second >= diff2secondThreshold:
                    geneFilterIn += 1
                    # then write to file
                    # distinguish between cases where optional columns were provided or not
                    if opt:
                        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,diffEst,pval,testStat,malig_mean,nonmal_max,nonmal_min,"\t".join(optVals)))
                    else:
                        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,diffEst,pval,testStat,malig_mean,nonmal_max,nonmal_min))
            # or if mean of malignant cluster is below nonmal_min
            elif malig_mean < nonmal_min:
                # and if the difference between the two is above the threshold
                diff2second = nonmal_min - malig_mean
                if diff2second >= diff2secondThreshold:
                    geneFilterIn += 1
                    # then write to file
                    # distinguish between cases where optional columns were provided or not
                    if opt:
                        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,diffEst,pval,testStat,malig_mean,nonmal_max,nonmal_min,"\t".join(optVals)))
                    else:
                        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,diffEst,pval,testStat,malig_mean,nonmal_max,nonmal_min))

infile.close()
outfile.close()

print("Parsed %s genes. %s genes passed the filter thresholds. Filtered %s NA genes." %(geneCount,geneFilterIn,filtered_NA))
