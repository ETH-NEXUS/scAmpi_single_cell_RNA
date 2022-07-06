# collect all drugs identified for all clusters (or clusters with a specified tag)
# roche tumor profiler single cell

import sys
import argparse

# TODO: generalize accepted input files
# TODO: enable tag-based input file seclection, for instance ony clusters with tumor cells.

'''
function definitions
'''


'''
main body
'''
parser = argparse.ArgumentParser(description='Collect drugs per cluster and create a new overview table that specifies drugs per cluster.')
parser.add_argument('--inputFolder', dest='inFolder', required=True, help='Folder with input files. This should be a table per cluster with drug information, e.g. resulting from clinical trials query, XXX.dgidb.txt.CompleteTable.ClinicalTrials.txt.')
parser.add_argument('--outFile', dest='outName', required=True, help='Name of the output file.')

args = parser.parse_args()

