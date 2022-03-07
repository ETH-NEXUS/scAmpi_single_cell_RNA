#!/usr/bin/env python

'''
Given a list of drug-gene interactions and a list of preferred Melanoma drugs this script
filters the first list and gives out only those drug-gene interactions of drugs 
that are also found in the second list

Anne Richter, April 2018
'''

import argparse
import re

parser = argparse.ArgumentParser(description='Filters drug-gene interaction entries after occurrence in the melanoma drug list.')
parser.add_argument('--inFile', dest='inFile', required=True, help='Input file of type {sample}.dgidb.txt.CompleteTable.ClinicalTrials.txt containing drug-gene interactions of DE genes')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name of output file.')
parser.add_argument('--drugList', dest='drugList', required=True, help='Drug list for cancer drugs that are clinically preferred, tab delimited file of type "orig-name,rank (comment),approved for melanoma (CH) ?"')

args = parser.parse_args()

# will be filled with drugs of clinical list as keys, other columns of clinical list as values
dictDrugs = {}
# clinical list of drugs
druglist = open(args.drugList, 'r')
headerLineDrug = druglist.readline()

for line in druglist:
    lineSplit = line.strip().split('\t')
    drugName = lineSplit[0]
    drugName = drugName.upper()
    dictDrugs.setdefault(drugName,[]).append(lineSplit[1])
    #dictDrugs[drugName].append(lineSplit[2])
#for key in dictDrugs:
#    print(key)
druglist.close()

# infile is list of drug-gene interactions of DE genes
infile = open(args.inFile, 'r')
outfile = open(args.outFile, 'w')
headerLine = infile.readline()
outfile.write(headerLine)

for line in infile:
    lineSplit = line.strip().split('\t')
    drugID = lineSplit[1]
    drugID = drugID.upper()
    # filtered list of drug-gene interactions including only drugs that are included in the clinical list as well is written to outfile
    if drugID in dictDrugs:
        outfile.write(line)
        continue
    else:
        # special case for Cisplatin that is as 'Cisplatin (Chembl2068237) in dgidb
        for key in dictDrugs:
            pattern = str(key) + '\s'
            pattern2 = '\s' + str(key)
            if re.search(pattern, drugID):
                outfile.write(line)
            elif re.search(pattern2, drugID):
                outfile.write(line)


infile.close()
outfile.close()
