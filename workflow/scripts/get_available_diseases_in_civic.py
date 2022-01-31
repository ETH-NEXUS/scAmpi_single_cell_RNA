#!/usr/bin/env python

'''
Retrieve list of available disease names in CIVIC cache file
Lourdes Rosano, Feb 2021
'''

import sys
import argparse
from civicpy import civic


'''
Script
'''

parser = argparse.ArgumentParser(description='Parse CIVIC cache file and extract list of available disease names.')
parser.add_argument('--outFile', dest='outFile', required=True, help='Output file listing disease names available in CIVIC cache file.')

args = parser.parse_args()

# Load offline CIVIC cache file
civic.load_cache(on_stale='ignore')

# Retrieve all gene records available in the cache file
all_results = civic.get_all_genes()

diseases = []
# Iterate all available gene records and retrieve all available disease names
for gene_record in all_results:
    gene_variants = gene_record.variants
    for variant_record in gene_variants:
        evidence_items = variant_record.evidence_items
        for evidence_record in evidence_items:
            # Use uppercase for consistency of disease names
            disease_name = evidence_record.disease.name.strip().upper()
            if disease_name not in diseases:
                diseases.append(disease_name)

print("Total # CIVIC diseases: %s" %(len(diseases)))
# Sort alphabetically
diseases.sort()

# Write output list
outfile = open(args.outFile, 'w')
for disease in diseases:
    outfile.write(disease + "\n")
outfile.close()
