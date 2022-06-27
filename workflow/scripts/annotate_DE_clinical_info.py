#!/usr/bin/env python

'''
Given a table with genes and variants (or expression levels etc), extend with columns regarding clinical information, e.g. dgidb drugs, clinical trials, pathways etc.
Franziska Singer, March 2018
'''

import sys
import numpy as np
import os
import re
import argparse

'''
function definitions
'''

def getPathwayInfo(pathwayDB):
    infilePathways = open(pathwayDB,'r')
    dictPathways = {}

    for linePW in infilePathways:
        lineSplit = linePW.strip().split()
        geneName = lineSplit[0].strip("\"")  # strip to remove the " symbols
        allPWs = "".join(lineSplit[1:])
        dictPathways[geneName] = allPWs.strip("\"")

    infilePathways.close()
    return dictPathways

def getDGIDB_categoryInfo(dgidbFile_categ):
    infileDGIDB_categories = open(dgidbFile_categ,'r')
    infileDGIDB_categories.readline()  # skip first line

    dictDGIDB_categ = {}
    dgidbGenes = 0
    for lineDGIDB in infileDGIDB_categories: # format: Gene    GeneName    Category
        lineSplit = lineDGIDB.strip().split("\t")
        gene = lineSplit[0]

        if gene in dictDGIDB_categ.keys():
            print("Error! %s already contained!" %(gene))
            continue

        dictDGIDB_categ[gene] = [lineSplit[1],lineSplit[2]]
        dgidbGenes += 1

    print("Read in %s genes from DGIDB query result." %(dgidbGenes))
    infileDGIDB_categories.close()
    return dictDGIDB_categ

def getClinicalTrialInfo(clinicalTrialsFile):
    infileTrials = open(clinicalTrialsFile,'r')
    dictClinTrials = {}

    drugNum = 0
    headerLine = infileTrials.readline()
    scoreColumn = 0
    clinTrialColumn = 0
    drugTypeColumn = 0

    headerSplit = headerLine.strip().split("\t")
    for headPos in range(0,len(headerSplit)):
        if headerSplit[headPos] == "Score":
            scoreColumn = headPos
            continue
        if headerSplit[headPos] == "Type":
            drugTypeColumn = headPos
            continue
        if headerSplit[headPos] == "ClinicalTrials (ID,isRecruiting y/n,Phase,isCancerTypeSpecific yes/no)":
            clinTrialColumn = headPos
            continue
    for lineTrials in infileTrials:
        lineSplit = lineTrials.strip().split("\t")
        gene = lineSplit[0]
        drug = lineSplit[1]
        score = lineSplit[scoreColumn]
        drugType = lineSplit[drugTypeColumn]
        clinicalTrialsInfo = lineSplit[clinTrialColumn]

        if gene not in dictClinTrials.keys():
            dictClinTrials[gene] = {}
        if drug not in dictClinTrials[gene].keys():
            dictClinTrials[gene][drug] = []
        else:
            print("Warning! Frug %s already contained for gene %s!" %(drug,gene))
            continue

        dictClinTrials[gene][drug].append(score)
        dictClinTrials[gene][drug].append(drugType)
        dictClinTrials[gene][drug].append(clinicalTrialsInfo)
        drugNum += 1

    infileTrials.close()
    print("Found %s drugs for dgidb genes in clinical trials file.\n" %(drugNum))
    return dictClinTrials

def getGeneColumnIndex(firstInputLine):
    firstLineSplit = firstInputLine.split('\t')
    index_geneCol = -1
    for pos in range(0,len(firstLineSplit)):
        if args.colName_gene in firstLineSplit[pos]:
            index_geneCol = pos
            break
    if index_geneCol == -1:
        print("Error! Did not find a column matching %s in input header %s." %(args.colName_gene,firstInputLine))
        sys.exit(1)

    return index_geneCol


parser = argparse.ArgumentParser(description='Annotate variants with clincial information.')
parser.add_argument('--inputTable', dest='inputFile', required=True, help='Input table with variants, expression levels, copy numbers.... Columns need to be tab separated.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name of the output file.')
parser.add_argument('--colName_gene', dest='colName_gene', required=True, help='Column name of column containing gene names')
parser.add_argument('--pathwayDB', dest='pathwayDB', required=False, help='Pathway database to annotae pathways to gene names. Optional.')
parser.add_argument('--clinTrials', dest='clinTrials', required=False, help='Clinical trials information annotated based on dgidb query. Optional.')
parser.add_argument('--dgidb_categ', dest='dgidb_categ', required=False, help='DGIdb gene categories, result of DGIdb query. Optional.')

args = parser.parse_args()

# pathway information
dictPathways = {}
if args.pathwayDB is not None:
    dictPathways = getPathwayInfo(args.pathwayDB)

# gene category information
dictDGIDB_categ = {}  # gene names to information
if args.dgidb_categ is not None:
    dictDGIDB_categ = getDGIDB_categoryInfo(args.dgidb_categ)

# gene drug and trials information
dictTrialInfo = {}  # gene names to information on drugs and clinical trials, including support of drug gene interaction
                # [score,type,clinicalTrials]
if args.clinTrials is not None:
    dictTrialInfo = getClinicalTrialInfo(args.clinTrials)

# sanity check: category and clinicalTrials file should include the same number of genes
if (args.dgidb_categ is not None) and (args.clinTrials is not None):
    if len(dictTrialInfo.keys()) != len(dictDGIDB_categ.keys()):
        print("Error! Different number of genes!")
    else:
        print("Number of dgidb genes: %s." %(len(dictDGIDB_categ.keys())) )

# finally parse existing variant table and create new file with more detailed information per gene

infile = open(args.inputFile,'r')
outfile = open(args.outFile,'w')
outfileIndependent = open(args.outFile + "_dgidbIndependent.txt",'w')

firstInputLine = infile.readline().strip()
outHeader = firstInputLine
if args.pathwayDB is not None:
    outHeader += "\tpathway"
if args.dgidb_categ is not None:
    outHeader += "\tDGIDB-drugs(Score,Type)\tDGIDB-geneName\tDGIDB-categories"
if args.clinTrials is not None:
    outHeader += "\tClinicalTrials_cancerType (phase,isRecruiting (y/n))\tClinicalTrials_notCancerType (phase,isRecruiting (y/n))"

outfile.write(outHeader + "\n")
outfileIndependent.write(outHeader + "\n")

index_geneCol = getGeneColumnIndex(firstInputLine)

matchedDGIDB = []

for line in infile:
    lineSplit = line.strip().split("\t")
    gene = lineSplit[index_geneCol]

    outfileLine = line.strip()

    # first include pathway

    if args.pathwayDB is not None:
        if gene in dictPathways.keys():
            outfileLine = outfileLine + "\t%s" %(dictPathways[gene])
        else:
            outfileLine = outfileLine + "\t."

    # dgidb and clinical trials
    # Note that clinical trials file is not available without dgidb file
    if args.dgidb_categ is None:
        outfile.write(outfileLine + "\n")
        outfileIndependent.write(outfileLine + "\n")
        continue

    includeInIndependentFile = True
    for dgidbGene in dictDGIDB_categ.keys():
        if (dgidbGene in gene) and (len(gene.strip()) == len(dgidbGene.strip())):
            includeInIndependentFile = False

            if dgidbGene not in matchedDGIDB:
                matchedDGIDB.append(dgidbGene)

            # drug and clincical trials info
            drugNames = ""
            clinicalTrialCancerType = ""
            clincialTrialNotCancerType = ""
            if dgidbGene not in dictTrialInfo.keys():
                print("Error. Gene %s not conatined in dgidb clincial trials dict!" %(dgidbGene))
                continue

            for drugInteractions in dictTrialInfo[dgidbGene].keys():
                drugNames += drugInteractions + "(" + dictTrialInfo[dgidbGene][drugInteractions][0] + "," + dictTrialInfo[dgidbGene][drugInteractions][1] + ");"
                clinTrialInfo = dictTrialInfo[dgidbGene][drugInteractions][2]

                foundCT = False
                foundNotCT = False
                numTrialsCT = 0
                numTrialsNotCT = 0
                clinicalTrialCancerType_temp = ""
                clincialTrialNotCancerType_temp_1 = "" # for the different phases. Only report trial IDs for the highest phase found for this drug
                clincialTrialNotCancerType_temp_2 = ""
                clincialTrialNotCancerType_temp_3 = ""
                if ";" in clinTrialInfo:
                    clinSplit = clinTrialInfo.strip().split(";")
                    for clinEntry in clinSplit:
                        trialInfo = clinEntry.split(",")
                        if "yes" in trialInfo[3]: 
                            clinicalTrialCancerType_temp += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
                            foundCT = True
                            numTrialsCT += 1
                        else:
                            if "3" in trialInfo[2]:
                                clincialTrialNotCancerType_temp_3 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
                            elif "2" in trialInfo[2]:
                                clincialTrialNotCancerType_temp_2 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
                            else:
                                clincialTrialNotCancerType_temp_1 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
                            foundNotCT = True
                            numTrialsNotCT += 1
                    if not foundCT:
                        clinicalTrialCancerType_temp += ".;"
                    if not foundNotCT:
                        clincialTrialNotCancerType_temp_1 += ".;"
                        clincialTrialNotCancerType_temp_2 += ".;"
                        clincialTrialNotCancerType_temp_3 += ".;"
                elif "," in clinTrialInfo:
                    trialInfo = clinTrialInfo.split(",")
                    if "yes" in trialInfo[3]:
                        clinicalTrialCancerType_temp += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
                        foundCT = True
                        numTrialsCT += 1
                    else:
                        # not cancer-type specific, check which phase
                        if "3" in trialInfo[2]:
                            clincialTrialNotCancerType_temp_3 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
                        elif "2" in trialInfo[2]:
                            clincialTrialNotCancerType_temp_2 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
                        else:
                            clincialTrialNotCancerType_temp_1 += trialInfo[0] + "(" + trialInfo[2] + "," + trialInfo[1] + ");"
                        foundNotCT = True
                        numTrialsNotCT += 1
                if not foundCT:
                    clinicalTrialCancerType_temp += ".;"
                if not foundNotCT:
                    clincialTrialNotCancerType_temp_1 += ".;"
                    clincialTrialNotCancerType_temp_2 += ".;"
                    clincialTrialNotCancerType_temp_3 += ".;"

                # now append everything to string for outfile
                clinicalTrialCancerType += drugInteractions + "(" + str(numTrialsCT) + "):" + clinicalTrialCancerType_temp
                if len(clincialTrialNotCancerType_temp_3) > 0:
                    clincialTrialNotCancerType += drugInteractions + "(" + str(numTrialsNotCT) + "):" + clincialTrialNotCancerType_temp_3
                elif len(clincialTrialNotCancerType_temp_2) > 0:
                    clincialTrialNotCancerType += drugInteractions + "(" + str(numTrialsNotCT) + "):" + clincialTrialNotCancerType_temp_2
                else:
                    clincialTrialNotCancerType += drugInteractions + "(" + str(numTrialsNotCT) + "):" + clincialTrialNotCancerType_temp_1

            outfileLine = outfileLine + "\t" + drugNames + "\t" + dictDGIDB_categ[dgidbGene][0] + "\t" + dictDGIDB_categ[dgidbGene][1]
            outfileLine = outfileLine + "\t" + clinicalTrialCancerType + "\t" + clincialTrialNotCancerType
            outfile.write(outfileLine + "\n")
            outfileIndependent.write(outfileLine + "\n")

    if includeInIndependentFile:
        outfileLine = outfileLine + "\tNA\tNA\tNA" # drug names and dgidb_categ info columns
        outfileLine = outfileLine + "\tNA\tNA" # clinical trial info
        outfileIndependent.write(outfileLine + "\n")

infile.close()
outfile.close()
outfileIndependent.close()

if args.dgidb_categ is not None:
    print("Matched %s of %s genes from dgidb." %(len(matchedDGIDB),len(dictDGIDB_categ.keys())))
