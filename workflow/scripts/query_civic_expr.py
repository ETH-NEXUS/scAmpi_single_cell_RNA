#!/usr/bin/env python

'''
Query CIViC DB for expression data
Version queries offline cache provided in civicpy module
(Cache needs to be updated by user to get latest data available)
Lourdes Rosano, Jul 2019
'''

import sys
import argparse
from civicpy import civic
import re
import os

# Load package and import relevant functions
sys.path.append("/cluster/work/nexus/antoine/Projects/2023_07_CIViCutils_add_header_option/civicutils/civicutils/")
#import civicutils
#from civicutils.read_and_write import read_in_expr, get_dict_support, write_match
#from civicutils.query import query_civic
#from civicutils.filtering import filter_civic
#from civicutils.match import match_in_civic, annotate_ct, filter_ct, process_drug_support

from read_and_write import read_in_expr, get_dict_support, write_match
from query import query_civic
from filtering import filter_civic
from match import match_in_civic, annotate_ct, filter_ct, process_drug_support

### Define global variable to ensure cache file is only loaded once even if several queries are performed
global isLoad
isLoad = False


'''
Script
'''

parser = argparse.ArgumentParser(description='Query CIViC to retrieve therapy information for expression.')
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


print(civic.__version__)
print(civic.LOCAL_CACHE_PATH)
print(os.getenv('CIVICPY_CACHE_FILE'))

## Sanity check for "empty" input lists (ie. in the form of [''])
## Turn them into "real" empty lists for implementation purposes
if cancerTypeList == ['']:
    cancerTypeList = []
if blackList == ['']:
    blackList = []
if highLevelList == ['']:
    highLevelList = []


# Read-in file of input SNV variants
(raw_data, expr_data, extra_header) = read_in_expr(args.inputFile, Gene_name=args.colName_gene, Logfc_name=args.colName_logFC)

# Query input genes in CIViC
var_map = query_civic(list(expr_data.keys()), identifier_type="entrez_symbol")

# Filter undesired evidences to avoid matching later on
var_map = filter_civic(var_map, evidence_status_in=["ACCEPTED"], var_origin_not_in=["GERMLINE"], output_empty=False)

# Match input SNV variants in CIViC, pick highest tier available per input gene+variant
# Tier hierarchy: 1 > 1b > 2 > 3 > 4
(match_map, matched_ids, var_map) = match_in_civic(expr_data, data_type="EXPR", identifier_type="entrez_symbol", select_tier="highest", var_map=var_map)

# Annotate matched CIViC evidences with cancer specificity of the associated diseases

annot_map = annotate_ct(var_map, blackList, cancerTypeList, highLevelList)

# Filter CIViC evidences to pick only those for the highest cancer specificity available
# ct hierarchy: ct > gt > nct
annot_map = filter_ct(annot_map, select_ct="highest")

# Get custom dictionary of support from data.yml (provided within the package)
# This defines how each combination of evidence direction + clinical significance in CIViC is classified in terms of drug response (e.g. sensitivity, resistance, unknown, etc.)
support_dict = get_dict_support()

# Process consensus drug support for the matched variants using the underlying CIViC evidences annotated 
annot_match = process_drug_support(match_map, annot_map, support_dict)

# Write to output
# Do not report the CT classification of each disease, and write column with the drug responses predicted for each available CT class of every variant match
write_match(annot_match, annot_map, raw_data, extra_header, data_type="EXPR", outfile=args.outFile, has_support=True, has_ct=True, write_ct=False, write_support=True, write_complete=False)




#print('\nTotal # matches: {}'.format(nMatches))
#print('Total # no matches: {}'.format(noMatches))
#print('Total # unavailable: {}'.format(notFound))
#print('---------------------')
# TODO: only report # of genes instead of whole list?
#print('Genes with no CIViC variant data:\n {}'.format(','.join(genesNotFound)))
#print('\nUnmatched expression:\n {}'.format('\n '.join(not_matched)))
