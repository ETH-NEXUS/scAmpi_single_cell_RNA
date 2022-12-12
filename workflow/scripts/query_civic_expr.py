#!/usr/bin/env python

'''
Query CIViC DB for expression data
Version queries offline cache provided in civicpy module
(Cache needs to be updated by user to get latest data available)
Antoine Hanns, Dec 2022
'''
import sys
import argparse
from civicpy import civic

## Load relevant functions from CIViCutils package

sys.path.append('/cluster/work/nexus/antoine/Projects/2022_12_integrate_Civicutils/git_files/civicutils/civicutils/')
from read_and_write import readInExpr,write_match
from query import query_civic
from filtering import filter_civic
from match import match_in_civic,annotate_ct,filter_ct,process_drug_support

'''
Script
'''

### Define global variable to ensure cache file is only loaded once even if several queries are performed
global isLoad
isLoad = False


parser = argparse.ArgumentParser(description='Query CIViC to retrieve drug information for expression.')
parser.add_argument('--inputTable', dest='inputFile', required=True, help='Input table with genes and logFC values, needs to be tab separated.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Name of the output file.')
parser.add_argument('--cancerTypeList', dest='cancerTypeList', required=True, help='Comma-separated list of accepted cancer types. Partial matches will be sought.')
parser.add_argument('--blackList', dest='blackList', required=True, help='Comma-separated list of not accepted cancer types. Partial matches will be sought, and those matched will be excluded from the list of considered evidence.')
parser.add_argument('--highLevelList', dest='highLevelList', required=True, help='Comma-separated list of high level cancer types (e.g. Cancer). Only exact matches will be sought. The user should be aware that results for high level cancer types will only be retrieved when no match for cancer specific types (ie. --cancerTypeList) is found.')
parser.add_argument('--colName_gene', dest='colName_gene', required=True, help='Name of column containing gene symbols.')
parser.add_argument('--colName_logFC', dest='colName_logFC', required=True, help='Name of column containing logFC values.')
parser.add_argument('--strictExpression', dest='strictExpression', required=True, choices=['y','n'], help='y/n. Specify whether a stricter CIVIC evidence interpretation should be applied for term *EXPRESSION* (y = apply). If yes, column name "pct_nonzero" will be required and EXPRESSION will be taken into account for interpretation only for genes where pct_nonzero > 0 (ie. gene is actually expressed in the given malignant cluster).')

args = parser.parse_args()

print("\nParameters:\n inputTable: %s\n outFile: %s\n colName_gene: %s\n colName_logFC: %s\n strictExpression: %s\n" %(args.inputFile,args.outFile,args.colName_gene,args.colName_logFC,args.strictExpression))

cancerTypeList = args.cancerTypeList
blackList = args.blackList
highLevelList = args.highLevelList
cancerTypeList = [s.strip().upper() for s in cancerTypeList.split(',')]
blackList = [s.strip().upper() for s in blackList.split(',')]
highLevelList = [s.strip().upper() for s in highLevelList.split(',')]

print('\nWhite listed cancer types: {}'.format(','.join(cancerTypeList)))
print('Black listed cancer types: {}'.format(','.join(blackList)))
print('High level cancer types: {}'.format(','.join(highLevelList)))

## Sanity check for "empty" input lists (ie. in the form of [''])
## Turn them into "real" empty lists for implementation purposes
if cancerTypeList == ['']:
    cancerTypeList = []
if blackList == ['']:
    blackList = []
if highLevelList == ['']:
    highLevelList = []


# Read in file of input SNV variants
(rawData,epxrData,extraHeader) = readInExpr(args.inputFile)

# Query input genes in CIVIC
varMap = query_civic(list(epxrData.keys()), identifier_type="entrez_symbol")
# Filter undesired evidences to avoid matching later on
varMap = filter_civic(varMap, evidence_status_in = ['ACCEPTED'], var_origin_not_in = ['GERMLINE'], output_empty=False)

# Match input SNV variants in CIVIC, pick highest tier available per input gene+variant
# Tier hierarchy: 1 > 1b > 2 > 3 > 4
(matchMap,matchedIds,varMap) = match_in_civic(epxrData, dataType="EXPR", identifier_type="entrez_symbol", select_tier="highest", varMap=varMap)

# Annotate matched CIVIC evidences with cancer specificity of the associated diseases
disease_name_not_in = blackList
disease_name_in = cancerTypeList
alt_disease_names = highLevelList
annotMap = annotate_ct(varMap, disease_name_not_in, disease_name_in, alt_disease_names)

# Filter CIVIC evidences to pick only those for the highest cancer specificity available
# ct hierarchy: ct > gt > nct
annotMap = filter_ct(annotMap,select_ct="highest")

# Get custom dictionary of support from data.yml (provided within the package)
# This defines how each combination of evidence direction + clinical significance in CIVIC is classified in terms of drug support (eg. sensitivity, resistance, unknown, etc.)
supportDict = {'SUPPORTS': {'SENSITIVITYRESPONSE':'POSITIVE', 'RESISTANCE':'NEGATIVE', 'REDUCED SENSITIVITY':'NEGATIVE', 'ADVERSE RESPONSE':'NEGATIVE'}, 'DOES_NOT_SUPPORT': {'RESISTANCE':'UNKNOWN_DNS', 'SENSITIVITYRESPONSE':'UNKNOWN_DNS', 'REDUCED SENSITIVITY':'UNKNOWN_DNS', 'ADVERSE RESPONSE':'UNKNOWN_DNS'}}

# Process drug support of the matched variants using the annotated CIVIC evidences
annotMatch = process_drug_support(matchMap,annotMap,supportDict)

# Write to output
# Report the CT classification of each disease, and write column with the overall drug support of the match for each available CT class
write_match(annotMatch, annotMap, rawData, extraHeader, dataType="EXPR", outfile=args.outFile, hasSupport=True, hasCt=True, writeCt=True, writeSupport=True, writeComplete=True)
