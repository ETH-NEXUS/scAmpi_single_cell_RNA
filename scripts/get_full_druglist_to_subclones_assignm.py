#!/usr/bin/env python

'''

This gets two inputs:
    - {sample}.coding_region_only.cell_cycle_removed.drugToCluster.allDrugs.txt, which has three tab separated columns: drug,clusters,weight and one line for each drug from the list where a drug-gene interaction was found to a differentially expressed gene in one of the malignant cell clusters
    - the list of drugs that are especially relevant for the given cancer indication, given by clinicians (e.g. /cluster/work/bewi/ngs/projects/tumorProfiler/code/all/scTumorProfiler/required_files/melanoma/melanoma_drug_list_converted_corrected.txt)

The output is a full list of the clinically relevant drugs (rows) with one column for each malignant cluster. With "1" or "0" it is indicated if a drug-gene interaction was found between the drug and any differentially expressed gene of the respective malignant cluster.

File name: get_full_druglist_to_cluster_assignm.py
Author: Anne Richter
November 2018

'''

import argparse

parser = argparse.ArgumentParser(description='Script that generates full list of clinically relevant genes and assignment to malignant clusters, "1" or "0" depending on if the drug interacts with any differentially expressed gene of this cluster.')
parser.add_argument('--in_drugToCluster', dest='in_drugToCluster', required=True, help='Input file of type {}.genes_cells_filtered.cell_cycle_removed.drugToCluster.filteredDrugs.txt')
parser.add_argument('--in_drugList', dest='in_drugList', required=True, help='Input list of relevant drugs')
parser.add_argument('--outFile', dest='outFile', required=True, help='Output file.')
args = parser.parse_args()


# open and read in drug list
drugListFile = open(args.in_drugList,'r')
header_drugListFile = drugListFile.readline()

drugList = []

for line in drugListFile:
    lineSplit = line.strip().split("\t")
    drugFromList = lineSplit[0].upper()
    drugList.append(drugFromList)
drugListFile.close()

# open and read in drugToCluster assignment
drugToClusterFile = open(args.in_drugToCluster,'r')
header_drugToClusterFile = drugToClusterFile.readline()

drugAssignmentCluster = {}
allMalignantClusterIDs = []

for line in drugToClusterFile:
    lineSplit = line.strip().split("\t")
    drugWithAssignm = lineSplit[0].upper()
    clusters = lineSplit[1]
    individualClusters = clusters.split(",")
    drugAssignmentCluster[drugWithAssignm] = individualClusters
    allMalignantClusterIDs = allMalignantClusterIDs + individualClusters
drugToClusterFile.close()

#for key,val in drugAssignmentCluster.items():
#    print(key, ":", val)

# get header names for output file
allMalignantClusterIDs_unique = list(set(allMalignantClusterIDs))
allMalignantClusterIDs_unique.sort()
outFile_headerNames = ['cluster_' + malignantClusterID for malignantClusterID in allMalignantClusterIDs_unique]
outFile_headerNames.insert(0, "drug")
outputHeader = "\t".join(outFile_headerNames)

# write header into outputFile
outFile = open(args.outFile,'w')
outFile.write(outputHeader + "\n")

# write lines into outputFile
for drug in drugList:
    outputLineItems = []
    outputLineItems.append(drug)
    # one column for each cluster
    for clusterID in allMalignantClusterIDs_unique:
        if drug in drugAssignmentCluster and clusterID in drugAssignmentCluster[drug]:
            # if drug-gene interaction between drug and diff. expr. gene in cluster: "1", else "0"
            outputLineItems.append("1")
        else:
            outputLineItems.append("0")
    outputLine = "\t".join(outputLineItems)
    outFile.write(outputLine + "\n")

outFile.close()
