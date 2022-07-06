#!/usr/bin/env python

'''
Given a gmt file from MSigDB, convert format to a maaping of each gene to its associated pathways
Franziska Singer, July 2018
'''

import sys
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description='Format gmt gene sets file.')
parser.add_argument('--inputGMT', dest='inputFile', required=True, help='Input gmt file from MSigDB.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name of the output file.')

args = parser.parse_args()

infile = open(args.inputFile,'r')
# lines are of the form: HALLMARK_ADIPOGENESIS	> Genes up-regulated during adipocyte differentiation (adipogenesis).	ABCA1	ABCB8	ACAA2	ACADL	ACADM	ACADS	ACLY	ACO2	ACOX1	ADCY6	ADIG	ADIPOQ	ADIPOR2	AGPAT3	AGPAT6	AIFM1	AK2	ALDH2
dictGenes = {}
allSets = []
allGenes = []

for line in infile:
	lineSplit = line.strip().split("\t")
	setName = lineSplit[0].strip()
	if setName in allSets:
		print("Warning! Set %s already seen before!" %(setName))
	allSets.append(setName)

	myGenes = lineSplit[2:]

	for gene in myGenes:
		if len(gene) == 0:
			continue
		if gene not in allGenes:
			allGenes.append(gene)
		if gene not in dictGenes.keys():
			dictGenes[gene] = []
		dictGenes[gene].append(setName)
	
infile.close()

allGenes = np.sort(allGenes)
outfile = open(args.outFile,'w')

for myGene in allGenes:
	mySets = ", ".join(dictGenes[myGene])
	outfile.write(myGene + "\t" + mySets + "\n")

outfile.close()

print("Read in %s gene sets and %s genes." %(len(allSets),len(allGenes)))
