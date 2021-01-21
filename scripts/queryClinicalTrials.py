#!/usr/bin/env python

'''
Query clinical trials.gov
assumes that all melanoma studies have been downloaded and unzipped already:
Example:
wget "https://clinicaltrials.gov/search?term=melanoma&studyxml=true" -O melanomaClinicalTrials_Query_2016-03-14.zip
unzip melanomaClinicalTrials_Query_2016-03-14.zip -d melanomaClinicalTrials_Query_2016-03-14

Franziska Singer, March 2016
'''

import sys
import numpy as np
import os
import re

def checkIdentifier(thisCondition):
	acceptedTags = ["ALL", "AML", "CANCER", "TUMOR", "TUMOUR", "NEOPLASM", "CARCINOMA", "GLIOBLASTOMA", "LEUKEMIA", "LEUCEMIA", "LEUKAEMIA", "LEUCAEMIA", "LYMPHOMA", "MELANOMA", "MALIGNANCY", "MALIGNANCIES"]

	for tag in acceptedTags:
		if tag in thisCondition:
			return True
	return False

def checkMeasureIdentifier(measure):
	acceptedTags = ["OS", "OR", "ORR", "PFS", "OVERALL SURVIVAL", "OBJECTIVE RESPONSE", "TUMOR RESPONSE", "DFS", "DISEASE-FREE SURVIVAL", "OBJECTIVE RESPONSE RATE", "PROGRESSION-FREE", "BEST OVERALL RESPONSE", "CLINICAL BENEFIT RATE", "DISEASE CONTROL RATE", "PROGRESSION FREE", "SURVIVAL RATE", "EVENT FREE SURVIVAL", "EVENT-FREE SURVIVAL"]

	for tag in acceptedTags:
		if tag in measure:
			return True
	return False

if len(sys.argv) <= 1:
	print("Query clinical trial xml files given drug gene interactions.")
	print("Usage: python queryClinicalTrials.py [inputInteractionTable] [outfile] [inputClinicalTrialsFolder] [cancerType_list] [cancerType_blackList]")
	sys.exit(1)
if "-h" in sys.argv[1]:
	print("Query clinical trial xml files given drug gene interactions.")
	print("Usage: python queryClinicalTrials.py [inputInteractionTable] [outfile] [inputClinicalTrialsFolder] [cancerType_list] [cancerType_blackList]")
	sys.exit(1)


inputInteractionTable = sys.argv[1]
outfileName = sys.argv[2]
inputClinicalTrialsFolder = sys.argv[3]
cancerTypeList = sys.argv[4]  # comma separated list of condition names that are accepted, e.g. for metastatic melanoma it could be "melanoma,solid tumor"
cancerType_blackList = sys.argv[5] # comma separated list of conditions that should not be counted as cancer type specific, e.g. for small cell lung cancer "non-small cell, small cell colon cancer"

cancerTypes = [s.upper() for s in cancerTypeList.split(",")]  # .upper() to ensure comparability to trial xml files
cancerTypes_notSpec = [s.upper() for s in cancerType_blackList.split(",")]

print("Input clinicalTrials folder: %s.\nList of cancer type specific conditions: %s.\nList of black listed conditions: %s.\n" %(inputClinicalTrialsFolder,str(cancerTypes),str(cancerTypes_notSpec)))

print("Default filter criteria:\nCONDITION matches Cancer, Carcinoma, Tumor, Lymphoma, Leukemia, or Neoplasm (or subsets and abbreviations)\nSTUDY_TYPE matches interventional\nPRIMARY PURPOSE matches treatment\nPRIMARY OUTCOME matches overall survival, objective response rate, clinical benefit rate, disease control rate, best overall response, tumor response, disease-free survival, objective response, or progression-free survival/duration of response (or subsets or abbreviations)")

# first parse all studies and store in dict
dictStudies = {}  # map drug names to study ids and recruiting status
studyCount = 0
filteredOut = 0
for file in os.listdir(inputClinicalTrialsFolder):
	xmlFile = "%s%s" %(inputClinicalTrialsFolder,os.path.basename(file))
	if os.path.isfile(xmlFile) and xmlFile.endswith(".xml"):
		tempXML = open(xmlFile,'r')
		#print(xmlFile)
		studyID = ""
		recruiting = "n"
		phase = ""
		isCancerTypeSpecific = "no" # in case of multiple conditions, at least one has to match in order to define this study as cancer type specific
		studyCount += 1
		drugName = ""
		filteredStudy = False
		found_one_good_condition = False
		found_one_good_measure = False
		for line in tempXML:
			if "<study_type>" in line:   # study type must match "interventional"
				studyType = line.strip().split(">")[1].split("<")[0].upper()
				if not "INTERVENTIONAL" in studyType:
					filteredOut += 1
					filteredStudy = True
					#print(studyType)
					break
				continue
			
			if "<primary_purpose>" in line:   # study type must match "interventional"
				purpose = line.strip().split(">")[1].split("<")[0].upper()
				if not "TREATMENT" in purpose:
					filteredOut += 1
					filteredStudy = True
					#print("purpose:" + purpose)
					break
				continue
			
			if "<nct_id>" in line:   # study ID identifier
				studyID = line.strip().split(">")[1].split("<")[0]
				if "NCT" not in studyID:
					print("Warning: check identifier %s for xml file %s." %(studyID,os.path.basename(file)))
				continue
			
			if "<overall_status>" in line:   # check recruiting status, distinguish recruiting and not recruiting studies
				if ("Recruiting" in line) and ("not" not in line):
					recruiting = "y"
				continue
			
			if "<phase>" in line:  # check phase of the study
				phases =  re.findall(r'\d+', line.strip().split(">")[1].split("<")[0])  # returns a list with all numbers found
				if len(phases) == 0:
					phase = "N/A"
				else:
					phase =  "/".join(phases)
				
				continue
				
			if "<condition>" in line: # compare to white and black list with conditions, to determine if study is cancer type specific or not
				thisCondition = line.strip().split(">")[1].split("<")[0].upper()
				foundIdentifier = checkIdentifier(thisCondition)
				if foundIdentifier:
					found_one_good_condition = True

				if (thisCondition in cancerTypes) and (thisCondition not in cancerTypes_notSpec):
					# this is a cancer type-specific trial
					isCancerTypeSpecific = "yes"
				continue
				
			if "<intervention_name>" in line:   # identifier for drug name
				drugName = line.strip().split(">")[1].split("<")[0].upper()
				continue
			
			if "<measure>" in line:   # primary outcome measure must match one of a priori specified tags
				measure = line.strip().split(">")[1].split("<")[0].upper()

				foundMeasureTag = checkMeasureIdentifier(measure)

				if foundMeasureTag:
					found_one_good_measure = True
				continue
		
		#after xml parsed, insert to dict
		if not filteredStudy:
			if not found_one_good_condition: # none of the conditions matched, so filter out the study
				filteredOut += 1
				filteredStudy = True
				continue
			if not found_one_good_measure: # none of the measures matched, so filter out the study
				filteredOut += 1
				filteredStudy = True
				continue
			if drugName not in dictStudies:
				dictStudies[drugName] = []
			dictStudies[drugName].append("%s,%s,%s,%s"%(studyID,recruiting,phase,isCancerTypeSpecific))
		
		tempXML.close()

print("Parsed %s studies. Filtered out %s studies." %(studyCount,filteredOut))


# now go through given drug gene interactions
infile = open(inputInteractionTable,'r')
outfile = open(outfileName,'w')   # write extended input table to outfile

allDrugs = 0
foundDrug = 0

for line in infile:
	if line.startswith("Gene"):
		outfile.write(line.strip() + "\tClinicalTrials (ID,isRecruiting y/n,Phase,isCancerTypeSpecific yes/no)\n")
		continue
	
	lineSplit = line.strip().split("\t")
	
	drugName = lineSplit[1]
	allDrugs += 1
	if drugName in dictStudies:
		foundDrug += 1
		outfile.write(line.strip() + "\t" + ";".join(dictStudies[drugName]) + "\n")   # append the information from dictionary to the table
	else:
		outfile.write(line.strip() + "\t.\n")   # nothing to append
		
infile.close()
outfile.close()

print("Found clinical trial information for %s of %s possible drugs." %(foundDrug,allDrugs))
