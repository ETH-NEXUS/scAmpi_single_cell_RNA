#!/usr/bin/env python

'''
Given a table (of type {sample}.coding_region_only.cell_cycle_removed.phenograph_celltype_association.txt) that sums up information about celltype classification and phenograph clustering this script per cluster determines for each malignant cluster: the number of cells in respective cluster, the total number of malignant cells and the percentage of malignant cells in respective cluster
Anne Richter, March 2018
'''

import argparse

parser = argparse.ArgumentParser(description='Gives out number of cells per malignant cluster and percentage of total cells per cluster.')
parser.add_argument('--inputTable', dest='inputFile', required=True, help='Input table {sample}.coding_region_only.cell_cycle_removed.phenograph_celltype_association.txt that sums up information of celltype classification and phenograph clustering.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name of output file.')
parser.add_argument('--malignant', dest='malignant', required=True, help='Term that is given as label to clusters with predominantly malignant celltype (e.g. melanoma)')

args = parser.parse_args()
infile = open(args.inputFile,'r')
outfile = open(args.outFile,'w')

# getting term for malignant from command line argument
mal = args.malignant

# extracting header line of input table
headerLine = infile.readline()

# find column indices that specify dominant celltype of cluster and sum of cells in it
sumCol = 0
celltypeCol = 0
headerSplit = headerLine.strip().split("\t")
for headPos in range(0, len(headerSplit)):
    if headerSplit[headPos] == '"Sum"':
        sumCol = headPos
    elif headerSplit[headPos] == '"Dominant.celltype"':
        celltypeCol = headPos

# introducing variables
sumCells = 0
dictClusters = {}

# parsing all other lines of input file
for line in infile:
    lineSplit = line.strip().split("\t")
    clusterID = lineSplit[0].strip('\"')
    # only consider rows with info about one cluster, not "Sum" or "Percent" row
    if clusterID.isdigit():
        celltype = lineSplit[celltypeCol].strip('\"')
        numberCells = int(lineSplit[sumCol].strip('\"'))
        # only consider rows of malignant clusters
        #if celltype == mal:
        if mal in celltype:  # Less strict matching to ensure that different tumor subtypes are recognized (e.g. melanoma_mesenchymal)
            dictClusters[clusterID] = numberCells
            sumCells += numberCells

outHeader = "cluster_ID\tnumber_cells_in_cluster\ttotal_number_malignant_cells\tpercentage_of_cells_in_cluster"
outfile.write(outHeader + "\n")

# write output line for each malignant cluster
for key in sorted(dictClusters):
    perc = ((float(dictClusters[key])) / sumCells)*100.0
    percForm = '{:.2f}'.format(perc)
    print("Cluster " +  key + " contains " + percForm + " percent of the cells in clusters labelled as malignant.")
    outfileLine = '{0}\t{1}\t{2}\t{3}'.format(key,dictClusters[key],sumCells,percForm)
    outfile.write(outfileLine + "\n")

infile.close()
outfile.close()

# print warning if none of the clusters in the input table has the tag given at --malignant
if sumCells == 0:
    print("\n========================================================\n" +
            "WARNING: the output file is empty except for the Header\n" +
            "This is probably because none of the clusters in the input table is declared to be of the type that was specified with '--malignant'\n" +
            "========================================================\n"
            )
