#!/usr/bin/env python

'''
Given a table of type drug,clusters,weight (with columns being tab separated and clusters being comma separated)
this script reformats the data to be of type drug,cluster1,cluster2...clustern with the entries being 1 or 0
depending on if the respective drug targets any gene in cluster1..n.
This preprocessing is needed for a visualisation of the data with upsetR plot

Anne Richter, April 2018
'''

import argparse
import re

parser = argparse.ArgumentParser(description='Script reformats table from type drug,clusters,weight to type drug,cluster1,cluster2...clustern so that it is fit for plotting in R')
parser.add_argument('--inFile', dest='inFile', required=True, help='Path to input table of type drug,clusters,weight')
parser.add_argument('--outFile', dest='outFile', required=True, help='Path to output table')

args = parser.parse_args()

# dict is filled with all lines of the input file with the drugName as key and list of targeted clusters as value
# this is done to generate a full list of occurring clusters in allclusters
dictDrugs = {}
allclusters = []

infile = open(args.inFile, 'r')
headerline = infile.readline()

for line in infile:
    lineSplit = line.strip().split('\t')
    # special case if drug column is empty (or weight column, one of the two outer lines)
    if len(lineSplit) < 3:
        continue
    drugName = lineSplit[0].upper()
    clusters = lineSplit[1]
    # special case if clusters column in the middle is empty
    if clusters:
        clustersSplit = clusters.split(',')
        # dictionary is filled with each line
        dictDrugs[drugName] = clustersSplit
        # list clustersSplit is integrated fully to allclusters list
        allclusters.extend(clustersSplit)
    else:
        dictDrugs[drugName] = []

infile.close()

# allclusters list is converted to set to remove redundancy
allclusters = sorted(set(allclusters))

# outfile opened for writing
outfile=open(args.outFile, 'w')
# header is written into outfile, one column for the drugName and one for each cluster
outfile.write('drug\t' + '\t'.join(allclusters) + '\n')

# for each drug a '1' or '0' entry is created corresponding to each cluster
# depending on if exists any DE gene in cluster that is targetable by drug
for key in dictDrugs:
    outline = []
    for element in allclusters:
        if element in dictDrugs[key]:
            outline.append(str(1))
        elif element not in dictDrugs[key]:
            outline.append(str(0))
    outfile.write(key + '\t' + '\t'.join(outline) + '\n')

outfile.close()
