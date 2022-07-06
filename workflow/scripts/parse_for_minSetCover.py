#!/usr/bin/env python

'''
Parse the output of {sample}.dgidb.txt.CompleteTable.ClinicalTrials.txt to generate
input table for minSetCover.py. Input file: database results of drugs interacting with DE genes
Output file: Tab-separated file with three columns drug,clusters,weight

Anne Richter, April 2018
Franziska Singer, September 2018
'''

import argparse
import os.path
import re

''' 
functions
'''

def readIn_drugList(dict_prioDrugs):
	drugFile = open(args.drug_list,'r')

	drugFile.readline()  # skip header line

	for drugLine in drugFile:
		lineSplit = drugLine.strip().split('\t')
		drug = lineSplit[0].upper()
		rank = int(lineSplit[1].split()[0])

		if drug not in dict_prioDrugs.keys():
			dict_prioDrugs[drug] = rank
		else:
			print("Warning! Drug %s is contained multiple times in drug priority file!" %(drug))
	drugFile.close()
	return(dict_prioDrugs)

'''
main body
'''

parser = argparse.ArgumentParser(description='Script reformats the table {sample}.dgidb.txt.CompleteTable.ClinicalTrials.txt to the input format of minSetCover.py')
parser.add_argument('--inFiles', dest='inFiles', required=True, nargs='+', help='Input tables {sample}.dgidb.txt.CompleteTable.ClinicalTrials.txt')
parser.add_argument('--outFile', dest='outFile', required=True, help='Path and name of output table')
parser.add_argument('--colName_clinTrial', dest='colNameClinTrial', required=True, help='Identifier for the column containing clinical trial info')
parser.add_argument('--colName_DGIDB_Score', dest='colName_DGIDB_score', required=True, help='Identifier for the column containing the DGIDB score')
parser.add_argument('--drug_list', dest='drug_list', required=True, help='File with prioritized drugs')
args = parser.parse_args()


dict_prioDrugs = {}
dict_prioDrugs = readIn_drugList(dict_prioDrugs)

# dictionary with drugName as key and a list as value
# list will contain :   [0] list with lineSplit
#                       [1] list of cluster IDs
#                       [2] list with information about cancer specificity in trials
#                       [3] drugScore
dictDrugs = {}

maxScore = int(0)

# iterate over input files
for file in args.inFiles:
	infile = open(file,'r')
	fileName = os.path.basename(file)
	# cluster name is read out of file name
	# cluster name is the only number that occurrs before the first '.' of the file name
	clusterName = fileName.split(".")[-7]
	headerLine = infile.readline()
	headerSplit = headerLine.strip().split('\t')
	index_clinTrialCol = -1
	indexDGIDB_score = -1
	for headPos in range(0,len(headerSplit)):
		if args.colNameClinTrial in headerSplit[headPos]:
			index_clinTrialCol = headPos
			continue
		if args.colName_DGIDB_score == headerSplit[headPos]:
			indexDGIDB_score = headPos
			continue
	if (index_clinTrialCol == -1) or (indexDGIDB_score == -1):
		print("Error! Did not find clinical trials column using identifiers: " + args.colNameClinTrial + " and " + args.colName_DGIDB_score)
		sys.exit(1)

	for line in infile:
		lineSplit = line.strip().split('\t')
		drugName = lineSplit[1].upper()
		# string that describes clinical trials is item of lineSplit
		trialsString = lineSplit[index_clinTrialCol]
		trialsList = trialsString.split(';')
		# making the trials strings accessible
		allcancerSpecificity = []
		for trial in trialsList:
			trialProperties = trial.split(',')
			# extract if any trial considers drug cancertype-specific
			if len(trialProperties) > 1:
				cancerSpecificity = trialProperties[3]
			else:
				cancerSpecificity = 'noTrial'
			allcancerSpecificity.append(cancerSpecificity)
		# for the weighting count if drug is cancertype-specific in any trial
		yescount = allcancerSpecificity.count('yes')

		# generate score for each drug
		# trialScore is 1 if no clinical trials are available, 2 if drug is never cancertype-specific and 3 if drug is cancertype specific in at least one trial
		trialScore = 0
		if cancerSpecificity[-1] == 'noTrial':
			trialScore = float(1.0)
		elif yescount == 0:
			trialScore = float(2.0)
		elif yescount > 0:
			trialScore = float(3.0)
		# score of dgidb
		scoreScore = float(lineSplit[indexDGIDB_score])
	
		prioRank_score = 1.0  # non-prioritized drugs have the smallest weight in terms of ranking

		if drugName in dict_prioDrugs.keys():
			rank = dict_prioDrugs[drugName]
			prioRank_score = float(len(dict_prioDrugs.keys())+2 - rank) # nonPrio_drugs = 1.0; lowest ranked drug has prioRna = 2.0 etc

		drugScore = trialScore * scoreScore * prioRank_score
		print("Drug: %s, trialScore: %s, dgidbScore: %s, prioRankScore: %s, drugScore: %s" %(drugName,trialScore,scoreScore,prioRank_score,drugScore))
		# maximum score is kept
		if maxScore < drugScore:
			maxScore = drugScore
		# One drug can show many drug-gene interactions with DE genes in any cluster
		# For the weight for minSetCover the best score of each drug with any DE gene in any clutser is kept
		# if dictDrugs[drugName] already exists and latest score is better elements [0], [2] and [3] are updated. Latest clusterName is added to [1]
		if drugName in dictDrugs:
			if dictDrugs[drugName][3] < drugScore:
				dictDrugs[drugName][0] = lineSplit
				dictDrugs[drugName][1].append(clusterName)
				dictDrugs[drugName][2] = allcancerSpecificity
				dictDrugs[drugName][3] = drugScore
			else:
				# latest clusterName is added to [1] in any case
				dictDrugs[drugName][1].append(clusterName)
		elif drugName not in dictDrugs:
			# new entry with drugName in dictDrugs
			dictDrugs.setdefault(drugName, []).append(lineSplit)
			dictDrugs[drugName].append([])
			dictDrugs[drugName][1].append(clusterName)
			dictDrugs[drugName].append(allcancerSpecificity)
			dictDrugs[drugName].append(drugScore)
	infile.close()

outfile = open(args.outFile,'w')
outHeader = "drug\tclusters\tweight"
outfile.write(outHeader + "\n")

#  per drugName one line in output file
for key in sorted(dictDrugs):
	clusterCol = ",".join(set(dictDrugs[key][1]))
	# the weight for minSetCover algorithm is maxscore divided by drugScore. Minimum (best) weight is therefore 1
	drugWeight = maxScore / dictDrugs[key][3]
	outfileLine = '{0}\t{1}\t{2}'.format(key,clusterCol,drugWeight)
	outfile.write(outfileLine + "\n")

outfile.close()
