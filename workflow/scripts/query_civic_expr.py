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


### Define what combination of directions and significances define 'POSITIVE' and 'NEGATIVE' support
supportDict = {'SUPPORTS': {'SENSITIVITYRESPONSE':'POSITIVE', 'RESISTANCE':'NEGATIVE', 'REDUCED SENSITIVITY':'NEGATIVE', 'ADVERSE RESPONSE':'NEGATIVE'}, 'DOES_NOT_SUPPORT': {'RESISTANCE':'UNKNOWN_DNS', 'SENSITIVITYRESPONSE':'UNKNOWN_DNS', 'REDUCED SENSITIVITY':'UNKNOWN_DNS', 'ADVERSE RESPONSE':'UNKNOWN_DNS'}}

### Define global variable to ensure cache file is only loaded once even if several queries are performed
global isLoad
isLoad = False


'''
Functions
'''

### Functions to query CIVIC and process result into structured dictionary

# Check that a given identifier type is contained in the list of allowed values
def sanity_check_identifier_type(identifier_type):
    if identifier_type not in ["entrez_id", "entrez_symbol", "civic_id"]:
        print("\nError! '%s' is not a valid identifier_type. Please provide one of 'entrez_id', 'entrez_symbol' or 'civic_id'!" %(identifier_type))
        sys.exit(1)
    return None


# Given a list of CIVIC gene records, parse and reformat into an structured dictionary
# Assume gene records already contain all information about associated variants and corresponding evidence items
def reformat_results(results, identifier_type):

    varMap = {}
    retrieved_genes = []    # keep track of genes that could be retrieved from CIVIC
    no_variants = []        # keep track of genes retrieved from CIVIC but with no variants available
    all_variants = []       # keep track of all variants retrieved from CIVIC

    # Check that id type corresponds to one of the allowed options
    ignore = sanity_check_identifier_type(identifier_type)

    # Iterate individual gene records to retrieve associated variants and evidence information
    for gene_record in results:
        # Retrieve all ids associated to this gene (CIVIC, entrez, symbol)
        gene_civic_id = str(gene_record.id)
        gene_id = str(gene_record.entrez_id)
        # Use uppercase for consistency of the gene symbols
        gene_symbol = gene_record.name.strip().upper()
        # Retrieve variants associated to this gene (can be empty)
        gene_variants = gene_record.variants

        # Use the provided gene id (civic, entrez or symbol) to uniquely identify genes
        if identifier_type == "civic_id":
            gene_key = gene_civic_id
        if identifier_type == "entrez_id":
            gene_key = gene_id
        if identifier_type == "entrez_symbol":
            gene_key = gene_symbol

        # Keep track of gene ids which could be retrieved from CIVIC
        if gene_key not in retrieved_genes:
            retrieved_genes.append(gene_key)

        # NOTE: it seems that only genes having at least 1 variant record available are included in the offline cache (eg. gene ADORA1 has no variants and is found via API but not in the cache)
        # Skip genes that do not have any variants available in CIVIC
        if not gene_variants:
            if gene_key not in no_variants:
                no_variants.append(gene_key)
            continue

        # Iterate variant records associated to the current gene
        # Retrieve all relevant info listed for each variant
        for variant_record in gene_variants:
            # Internal variant id in CIVIC
            variant_id = str(variant_record.id)
            # Keep track of all variants retrieved from CIVIC
            if variant_id not in all_variants:
                all_variants.append(variant_id)
            # Variant name in CIVIC; use uppercase for consistency
            variant_name = variant_record.name.strip().upper()
            hgvs_expressions = variant_record.hgvs_expressions
            # Score to assess the accumulation of evidence for each variant (quantity and quality)
            civic_score = variant_record.molecular_profiles[0].molecular_profile_score
            # Sanity check for empty scores
            if civic_score is None:
                civic_score = "NULL"

            # List of evidence records available for the current variant (can be empty)
            evidence_items = variant_record.molecular_profiles[0].evidence_items

            # Iterate through the listed evidence items and store relevant information for this variant
            # Variants which do not have any clinical data associated to them will be directly skipped
            for evidence_record in evidence_items:

                # TODO: add sanity check for all relevant elements being present and non-empty (entrez_name, variant name, evidence_type, disease name). For now, assume this should never happen.

                # Use uppercase for consistency of the tags
                evidence_status = evidence_record.status.strip().upper()
                evidence_type = evidence_record.evidence_type.strip().upper()
                # TODO: extend query to be applicable to new evidence types 'ONCOGENIC' and 'FUNCTIONAL'
                # Skip for now, as they trigger errors due to empty disease names (also, they have not been fully tested)
                if evidence_type not in ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']:
                    continue
                disease = evidence_record.disease.name.strip().upper()

                evidence_level = evidence_record.evidence_level                     # just in case, should never be None
                variant_origin = evidence_record.variant_origin                     # can be None
                evidence_direction = evidence_record.evidence_direction             # can be None
                significance = evidence_record.significance       # can be None

                # Skip records that are not accepted evidence
                if (evidence_status != "ACCEPTED"):
                    continue

                # Sanity check for empty evidence direction, significance or level
                # 'NULL' is introduced to distinguish from 'N/A' tag
                if evidence_direction is None:
                    evidence_direction = "NULL"
                else:
                    evidence_direction = evidence_direction.strip().upper()
                if significance is None:
                    significance = "NULL"
                else:
                    significance = significance.strip().upper()
                if evidence_level is None:
                    evidence_level = "NULL"
                else:
                    evidence_level = evidence_level.strip().upper()

                # Combine the direction and significance of the evidence in one term
                evidence = evidence_direction + ':' + significance

                # At this point, currently evaluate evidence item has passed all the checks/filters
                # Keep track of the current evidence item under the corresponding variant and gene
                if gene_key not in varMap.keys():
                    varMap[gene_key] = {}

                # Variant name should be unique within gene (found some duplicates but all were submitted, not accepted data)
                if variant_name not in varMap[gene_key].keys():
                    varMap[gene_key][variant_name] = {}
                    varMap[gene_key][variant_name]['id'] = variant_id
                    varMap[gene_key][variant_name]['civic_score'] = civic_score

                    # Keep original HGVS annotations (empty list when nothing is available)
                    # Use uppercase to avoid mismatches due to case
                    varMap[gene_key][variant_name]['hgvs'] = [h.strip().upper() for h in hgvs_expressions]

                # FIXME: there is no sanity check for detecting possible variant name duplicates. For now, assume this should never happen.
                if evidence_type not in varMap[gene_key][variant_name].keys():
                    varMap[gene_key][variant_name][evidence_type] = {}
                if disease not in varMap[gene_key][variant_name][evidence_type].keys():
                    varMap[gene_key][variant_name][evidence_type][disease] = {}

                therapies = []
                evidence_therapies = evidence_record.therapies
                for evidence_therapy in evidence_therapies:
                    therapy_name = evidence_therapy.name.strip().upper()
                    if therapy_name not in therapies:
                        therapies.append(therapy_name)

                # When more than 1 therapy are listed for the same evidence item, 'therapy_interaction_type' is not null and defines the nature of this multiple therapy entry
                therapy_interaction = evidence_record.therapy_interaction_type
                if therapy_interaction is not None:
                    therapy_interaction = therapy_interaction.strip().upper()
                    # 'Substitutes' indicates that therapies can be considered individually
                    if therapy_interaction != "SUBSTITUTES":
                        # Remaining terms ('Sequential' and 'Combination') indicate that therapies should be considered together, so join their names into a single tag
                        # Sort therapies alphabetically to ensure that their order in the combination treatment is always the same
                        therapies.sort()
                        therapies = ["+".join(therapies)]

                if not therapies:
                    # Only non-Predictive evidences and Predictive ones without therapies will have this dummy level
                    # Introduced for consistency purposes within the varMap structure
                    therapies = ["NULL"]

                # Iterate through therapies to add evidences associated to them
                #   For non-Predictive evidences or Predictive with empty therapies, therapies=['NULL']
                #   For Predictive and interaction=None, len(therapies) = 1
                #   For Predictive and interaction='Substitutes', len(therapies)>1
                #   For Predictive and interaction!='Substitutes', len(therapies)=1 (combiantion of several using '+')
                for therapy in therapies:
                    if therapy not in varMap[gene_key][variant_name][evidence_type][disease].keys():
                        varMap[gene_key][variant_name][evidence_type][disease][therapy] = {}
                    if evidence not in varMap[gene_key][variant_name][evidence_type][disease][therapy].keys():
                        varMap[gene_key][variant_name][evidence_type][disease][therapy][evidence] = {}
                    if evidence_level not in varMap[gene_key][variant_name][evidence_type][disease][therapy][evidence].keys():
                        varMap[gene_key][variant_name][evidence_type][disease][therapy][evidence][evidence_level] = []
                    # Group all publications associated to the same level. Do not check publication status
                    ## On 25.01.2019, source structure was changed to introduce ASCO abstracts as a source type
                    # FIXME: there is no sanity check for empty ID, however assume this should never happen
                    # FIXME: check for type of source? currently, both PUBMED and ASCO would be considered
                    varMap[gene_key][variant_name][evidence_type][disease][therapy][evidence][evidence_level].append(str(evidence_record.source.citation_id).strip())

    # TODO: assertions are currently not considered, as they mostly consist of free text summarizing the available evidence

    return (varMap,retrieved_genes,no_variants,all_variants)


# Given a list of gene identifiers, query CIVIC for known variants and return a structured dictionary with the relevant results
# List of gene ids can be: CIVIC id, entrez id or gene symbol
def query_civic_genes(genes, identifier_type="entrez_symbol"):

    # Check that provided argument is a list (even if length = 1)
    if (not isinstance(genes, list)) or (not genes):
        print("\nError! Please provide a list of gene ids to query in CIVIC.")
        sys.exit(1)

    # Check that id type corresponds to one of the allowed options
    ignore = sanity_check_identifier_type(identifier_type)

    ## Load offline cache of CIVICdb
    # Ensure only loaded once
    # NOTE: it is important that on_stale='ignore' to avoid attempts to update the cache via internet access
    global isLoad
    if not isLoad:
        print("Loading offline CIVIC cache...")
        is_successful = civic.load_cache(on_stale='ignore')
        if is_successful:
            isLoad = True
        else:
            print("\nError! Local cache file could not be successfully loaded!")
            sys.exit(1)

    # Querying functionality provided by module `civicpy` can only use internal CIVIC ids
    # NOTE: when using `civic.get_genes_by_ids()` with CIVIC ids, if any of them is not contained in the cache file, then it will try to directly query the CIVICdb, causing this script to crash

    # Workaround: retrieve all gene records available in the offline cache and parse them to match to input genes
    # Avoid at all costs directly querying the CIVIC db: even for CIVIC ids, we will match to available records from the offline cache
    all_results = civic.get_all_genes()

    # Iterate individual records and retrieve only those matching the provided gene ids
    results = []
    for gene_record in all_results:
        toKeep = False

        if (identifier_type == "civic_id"):
            gene_id = str(gene_record.id)                    # expectation is a single number
            if gene_id in genes:
                toKeep = True

        if (identifier_type == "entrez_id"):
            gene_id = str(gene_record.entrez_id)             # expectation is a single number
            if gene_id in genes:
                toKeep = True

        if (identifier_type == "entrez_symbol"):
            gene_id = gene_record.name.strip()          # expectation is single string
            aliases = gene_record.aliases               # expectation is list of strings
            # Perform union of all gene symbols available for the current record
            tmp_symbols = list(set([gene_id]) | set(aliases))
            # Try to match current gene record (any alias) to the provided gene symbols
            for tmp_symbol in tmp_symbols:
                # Use uppercase to ensure consistency of gene symbols
                this_symbol = tmp_symbol.strip().upper()
                if this_symbol in genes:
                    toKeep = True
                    break

        if toKeep:
            results.append(gene_record)

    # At this point, all CIVIC results for queried genes have been retrieved in a list
    # Process gene records into a dictionary with structured format
    # gene -> variants -> evidence_items
    (dict_results,retrieved_genes,no_variants,all_variants) = reformat_results(results, identifier_type)
    return (dict_results,retrieved_genes,no_variants,all_variants)


### Other functions

## Given the header, find the input columns containing genes and logFC values (use provided arguments directly)
## When strictExpression = "y", column name "pct_nonzero" will be searched for as well (name is implicitly assumed).
def getColumnIndices(firstInputLine):
    firstLineSplit = firstInputLine.split('\t')
    index_geneCol = -1
    index_logCol = -1
    index_pct = -1
    for pos in range(0,len(firstLineSplit)):
        # Avoid mismatches due to case by always using uppercase
        if args.colName_gene.upper() == firstLineSplit[pos].upper():
            index_geneCol = pos
        if args.colName_logFC.upper() == firstLineSplit[pos].upper():
            index_logCol = pos
        # Only when input argument strictExpression = "y", search for 1 additional column
        if args.strictExpression == "y":
            # Name of required column is assumed to be "pct_nonzero" (default as it is implicitly defined within the pipeline)
            if "pct_nonzero" == firstLineSplit[pos].lower(): 
                index_pct = pos

    if (index_geneCol == -1) or (index_logCol == -1):
        print("Error! Could not match all input columns in header %s." %(firstInputLine))
        sys.exit(1)
    # Additional check for column "pct_nonzero" when strictExpression = "y"
    if args.strictExpression == "y":
        if index_pct == -1:
            print("Error! Could not find input column 'pct_nonzero' in header %s." %(firstInputLine))
            sys.exit(1)

    return (index_geneCol,index_logCol,index_pct)


## Given a single CIVIC record name, return whether it corresponds to a EXPRESSION record related to exons, and the type of expression change
## For this, attemp to match the variant name to special EXPRESSION exon cases present in CIVIC (eg. EXON 1-2 EXPRESSION, EXON 5 OVEREXPRESSION...)
def expr_is_exon_string(varName):
    exprType = ''
    isExon = False
    if re.search('^EXON [0-9-]+ EXPRESSION$', varName):
        isExon = True
        exprType = 'EXPRESSION'
    elif re.search('^EXON [0-9-]+ OVEREXPRESSION$', varName):
        isExon = True
        exprType = 'OVEREXPRESSION'
    elif re.search('^EXON [0-9-]+ UNDEREXPRESSION$', varName):
        isExon = True
        exprType = 'UNDEREXPRESSION'

    return (isExon,exprType)


# Write information about a single cancer type item into one or more structured strings
# i.e. DISEASE[|therapy1,therapy2..](direction, significance(level(PMID,..,PMID),level(..)));
## For 'predictive' evidence (writetherapy=True), keep dictionary of therapy support for the current gene/line
## Format: therapy -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
def write_evidences(item, cancer, ct, therapyMap, writetherapy=False):
    evidences = []
    # For each therapy associated with the given cancer type
    for therapy in item.keys():
        # For each evidence associated with the given therapy
        # Evidences are simplified by using the combined form 'direction:significance'
        for evidence in item[therapy].keys():
            ## For each evidence (ie combination of direction+clin_signf), count how many different evidence items support it
            ## At this stage, we find count evidence items by counting how many different combinations of level+pmids there are for the same therapy, disease and evidence
            pmids = []
            # If therapy=True, write therapy information, i.e. DISEASE|therapy(..)
            if writetherapy:
                # Always one therapy (single or combination with '+')
                outString = cancer + '|' + therapy + '('
            else:
                outString = cancer + '('
            # Split the evidence direction and significance
            direction, clin_signf = evidence.split(':')
            outString += direction + ',' + clin_signf + '('
            # There may be several levels grouped per evidence
            levels = []
            for level in item[therapy][evidence].keys():
                # There may be several publications (i.e. PMIDs) grouped per level
                levels.append(level + '(' + ','.join(item[therapy][evidence][level]) + ')')
                # Count how many different evidence items support this particular evidence item
                for z in item[therapy][evidence][level]:
                    # Distinguish cases where the same publication is used to support different and identical evidence levels (they nonetheless count as separate evidence items)
                    pmids.append(z)
#                     new_z = level + '_' + z
#                     if new_z not in pmids:
#                         pmids.append(new_z)
            outString += ','.join(levels) + '))'
            evidences.append(outString)

            ## For 'predictive' evidence, keep dictionary of therapy support for the current gene/line
            ## Format: therapy -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
            if writetherapy:
                ## Only add current therapy and ct to dictionary if this match should indeed be considered for interpretation of the therapy support
                if therapy not in therapyMap.keys():
                    therapyMap[therapy] = {}
                ## Possible values for cancer specificity: 'ct' (specific), 'gt' (general), 'nct' (non-specific)
                if ct not in therapyMap[therapy].keys():
                    therapyMap[therapy][ct] = []

                ## Evaluate therapy support based on the combination of direction + significance
                ## Each combination has an associated support: POSITIVE, NEGATIVE, UNKNOWN_DNS or UNKNOWN_BLANK
                if ('NULL' in direction) or ('N/A' in direction) or ('NULL' in clin_signf) or ('N/A' in clin_signf):
                    thisSupport = 'UNKNOWN_BLANK'
                else:
                    if direction not in supportDict.keys():
                        print("Error! Could not find direction %s in support dictionary." %(direction))
                        sys.exit(1)
                    if clin_signf not in supportDict[direction].keys():
                        print("Error! Could not find significance %s in support dictionary." %(clin_signf))
                        sys.exit(1)
                    thisSupport = supportDict[direction][clin_signf]

                ## Keep track of number of occurrences for each support type for the given therapy
                ## Here, take into account the number of supporting PMIDs associated to each evidence item
                for z in pmids:
                    therapyMap[therapy][ct].append(thisSupport)

    return (evidences,therapyMap)


# Return in a list all evidence items (in written form) for a given evidence type 
# i.e. DISEASE1[|therapy1,therapy2..](direction1,significance1(level1(PMID,..,PMID),level2(..)));
#      DISEASE1[|therapy1,therapy2..](direction2,significance2(level1(PMID,..,PMID),level2(..)));
# For each evidence type, return either:
#   - Info on white listed cancer types (eg. 'Melanoma')
#   - If previous is not available, info on high level cancer types (eg. 'Cancer', 'Solid tumor')
#   - If previous is not available, info on all available cancer types for the given variant (except those included in the black list)
## For 'predictive' evidence (writetherapy=True), keep dictionary of therapy support for the current gene/line
## Format: therapy -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
def match_cancer_and_return_evidences(vardata, evidence_type, cancerTypes, cancerTypes_notSpec, highLevelTypes, therapyMap, writetherapy=False):
    ## Evidences returned for the matched cancer specificity
    evidences = []

    ## Keep track of cancer types that did not pass the black-list filter
    blackMatched = []

    ## Keep track of cancer types that passed the black-list filter
    cleanSet = []
    ## Final list of matched cancer types (associated to either 'ct', 'gt' or 'nct')
    matched = []
    ## Keep track of which type of cancer specificy was matched in the end ('ct', 'gt' or 'nct')
    ct = ''

    # If the given evidence type is not present for the variant, no evidences will be returned by the function
    if evidence_type in vardata.keys():

        ## 1) First, remove cancer types that partially match black-listed terms (if any are provided)
        # NOTE: PARTIAL matches to the black list are allowed! eg:
        #   - including 'small' will remove 'non-small cell lung cancer' and 'lung small cell carcinoma'
        #   - including 'non-small' will remove 'non-small cell lung cancer' but not 'lung small cell carcinoma'
        if cancerTypes_notSpec:
            # Iterate available cancer types and keep track of those that partially match to black list (there can be several)
            for cancerType in vardata[evidence_type].keys():
                # To find partial matches, it is necessary to iterate through both lists (input and civic)
                for blackListed in cancerTypes_notSpec:
                    # Search for partial match of INPUT cancer type in CIVIC cancer type e.g. 'Melanoma' (input list) in 'Skin Melanoma' (CIVIC) and not opposite
                    if blackListed in cancerType:
                        if cancerType not in blackMatched:
                            blackMatched.append(cancerType)
            # Iterate available cancer types once again to retrieve those that passed the black list filter
            for cancerType in vardata[evidence_type].keys():
                # Retrieve valid cancer types only
                if cancerType in blackMatched:
                    continue
                if cancerType not in cleanSet:
                    cleanSet.append(cancerType)

        ## If no black list was provided, then all available cancer types constitute the clean set
        else:
            cleanSet = list(vardata[evidence_type].keys())

        ## 2) Now, iterate the list of "allowed" cancer types (ie. passing the black list filter) and attempt to match to white-listed terms
        # NOTE: again, PARTIAL matches to the white list are allowed! eg:
        #   - including 'melanoma' will match 'melanoma', 'skin melanoma' and 'uveal melanoma', but not 'skin cancer' (hypothetical)
        #   - including 'uveal melanoma' will only match 'uveal melanoma'
        for cleanType in cleanSet:
            # To find partial matches, it is necessary to iterate through both lists (input and civic)
            for whiteListed in cancerTypes:
                # Search for partial match of INPUT cancer type in CIVIC cancer type e.g. 'Melanoma' (input list) in 'Skin Melanoma' (CIVIC) and not opposite
                if whiteListed in cleanType:
                    # Keep track of cancer types that passed the white list filter
                    if cleanType not in matched:
                        ct = 'ct'
                        matched.append(cleanType)

        ## 3) When nothing could be matched to the white-list terms, attempt a second-best match strategy to higher-level cancer types
        ## In CIVIC, some 'general' cancer types are included as disease, eg. 'cancer' or 'solid tumor'. Hence, must be exact matches since they are DB-specific
        # NOTE: here, only PERFECT matches are allowed!
        #   - including 'cancer' will only match 'cancer' and not 'lung cancer'
        if not matched:
            for cleanType in cleanSet:
                if cleanType in highLevelTypes:
                    if cleanType not in matched:
                        ct = 'gt'
                        matched.append(cleanType)

        ## 4) If nothing could be matched to either the white-list or higher-level cancer types, return all 'allowed' (ie. not black-listed) cancer types available in CIVIC
        ## These will be considered as non-specific cancer types (ie. off label)
        if not matched:
            for cleanType in cleanSet:
                if cleanType not in matched:
                    ct = 'nct'
                    matched.append(cleanType)

        ## Now, list 'matched' contains all cancer types for which their evidences items will be returned
        ## They can correspond to either 'ct' (from white-list), 'gt' (from high-level list) or 'nct' (if nothing could be matched, return all that is available)
        for cancerType in matched:
            ## Return evidence items in already formatted strings
            ## Also, return dictionary of therapy support for the current gene/line
            (strings,therapyMap) = write_evidences(vardata[evidence_type][cancerType], cancerType, ct, therapyMap, writetherapy)
            for s in strings:
                evidences.append(s)

    return (evidences,therapyMap)



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


infile = open(args.inputFile,'r')
firstInputLine = infile.readline().strip()

## Retrieve indices for all relevant columns (ie. gene and logFC)
## If not all are found, script will exit with an error
# When strictExpression = "y", additional column "pct_nonzero" will be sought as well (name is assumed).
# When strictExpression = "n", index_pct = -1 and will be ignored.
(index_geneCol,index_logFC,index_pct) = getColumnIndices(firstInputLine)

## Already write new header into output file
outfile = open(args.outFile,'w')
outHeader = firstInputLine
outHeader += "\tCIViC_Score\tCIViC_therapy_Support\tCIViC_Predictive\tCIViC_Diagnostic\tCIViC_Prognostic\tCIViC_Predisposing"
outfile.write(outHeader + "\n")

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


### Retrieve all input genes to query CIViC for variant info
genes = []
naGenes = []
for line in infile:
    lineSplit = line.strip().split("\t")
    ## Sanity check of correct format: 1 gene per line
    if ';' in lineSplit[index_geneCol]:
        print("Warning! Multiple genes %s detected within the same line in input file." %(lineSplit[index_geneCol]))
        sys.exit(1)
    # Uppercase as a sanity check
    gene = lineSplit[index_geneCol].strip().upper()
    logfc_annot = lineSplit[index_logFC].strip().upper()
    # Check gene is not empty as a sanity check
    if gene:
        ## Skip genes with logFC='NA' (ie could not asses differential expression), since empty values '.' will be assigned in CIVIC columns anyway
        if logfc_annot == 'NA':
            if gene not in naGenes:
                naGenes.append(gene)
            continue
        if gene not in genes:
            genes.append(gene)
infile.close()

# When no genes are found, write empty files (necessary for snakemake pipeline) and exit without error
# eg. empty input file containing only header because patient had no DE genes at all
if (not genes) and (not naGenes):
    print("\nDid not find any genes in column \'{}\' of file {}.".format(args.colName_gene, args.inputFile))
    outfile.close()
    sys.exit(0)
# When no genes were found because all existing were skipped due to logFC='NA',
# write output identical to input file plus additional CIVIC columns being empty
elif (not genes) and naGenes:
    print("\nAll genes parsed have logFC=NA.")
    ## Write corresponding "empty" output file: identical to input file + empty new CIVIC columns
    infile = open(args.inputFile,'r')
    next(infile)
    for line in infile:
        lineSplit = line.strip().split("\t")
        outFileline = "\t".join(lineSplit) + '\t.\t.\t.\t.\t.\t.\n'
        outfile.write(outFileline)
    infile.close()
    outfile.close()
    sys.exit(0)


### Query CIViC for the genes of interest

print('\nTotal # genes to query: {}'.format(len(genes)))
print('\nRetrieving data from CIViC...')

(varMap,retrieved_genes,no_variants,all_variants) = query_civic_genes(genes, identifier_type="entrez_symbol")

retrieved = set(retrieved_genes)
unmatched = list(set(genes) - retrieved)
print("Found %s/%s genes in CIVIC associated to %s variants. Found %s CIVIC genes that had no variants available." %(len(retrieved_genes),len(genes),len(all_variants),len(no_variants)))
print('\nGenes with no CIVIC data: {}'.format(','.join(unmatched)))


### Iterate through input table once more and create output table
### Only new format of input table is allowed!

# Keep track of all matches and non-matches
nMatches = 0               # no. lines where a CIVIC match was found
noMatches = 0              # no. lines where a CIVIC match was not found
notFound = 0               # no. lines where gene was not in CIVIC
not_matched = []           # gene+expression where gene was found in CIVIC but not matched to expression
genesNotFound = []         # Genes not found in CIVIC

print('Matching input data to CIVIC records...')
infile = open(args.inputFile,'r')
next(infile)
for lineIndx,line in enumerate(infile):
    lineSplit = line.strip().split("\t")
    ## Output line will consist of current line + additional columns
    outfileLine = line.strip()

    ## There should only be 1 gene and 1 logFC value per line
    ## Sanity check of correct format already done when parsing genes for CIVIC query
    gene = lineSplit[index_geneCol].strip().upper()
    logfc_annot = lineSplit[index_logFC].strip().upper()
    ## When strictExpression = "y", additional column "pct_nonzero" will be evaluated as well
    ## This argument specifies whether a stricter CIVIC evidence interpretation should be applied for term *EXPRESSION*
    ## If yes, EXPRESSION will be taken into account for interpretation only for genes where pct_nonzero > 0 (ie. gene is actually expressed in the given malignant cluster)
    if args.strictExpression == "y":
        pct = float(lineSplit[index_pct].strip())

    matched = []  # matches in CIVIC
    ## If gene has logFC='NA' (ie could no asses differential expression), skip and write empty values '.' for the new output columns
    if logfc_annot != 'NA':

        ## If gene is in CIVIC, attempt match to any CIVIC expression record
        if gene in varMap.keys():
            logfc = float(logfc_annot)

            ## Expression tags will be used to match expression records in CIVIC
            expr_change = ''
            input_tags = []
            ## Generate expression tags to be able to match record in CIVIC
            ## As of 29/09/2019, only relevant records in CIVIC are: 'OVEREXPRESSION', 'EXPRESSION', 'UNDEREXPRESSION' (ignore remaining)
            if logfc > 0:
                input_tags.append('OVEREXPRESSION')
                expr_change = 'OVEREXPRESSION'
            elif logfc < 0:
                input_tags.append('UNDEREXPRESSION')
                expr_change = 'UNDEREXPRESSION'
            else:
                print("Warning! Unknown logFC value for gene %s (logFC = %s)." %(gene,logfc_annot))
                sys.exit(1)

            ## When strictExpression = "y", additional column "pct_nonzero" will be evaluated
            ## EXPRESSION will be taken into account for interpretation only for genes where pct_nonzero > 0
            if args.strictExpression == "y":
                if pct > 0:
                    input_tags.append('EXPRESSION')
            ## When strictExpression = "n", always take EXPRESSION into account
            else:
                input_tags.append('EXPRESSION')

            ## Special case for gene 'CDKN2A': synonym 'p16' is used in CIVIC records (eg. 'p16 EXPRESSION')
            ## Include these cases as well to be able to match them
            new_tags = []
            if gene == 'CDKN2A':
                for this_tag in input_tags:
                    new_tag = 'P16 ' + this_tag
                    if new_tag not in new_tags:
                        new_tags.append(new_tag)
            for new_tag in new_tags:
                if new_tag not in input_tags:
                    input_tags.append(new_tag)


            ## Possible tag combinations are either OVEREXPRESSION,EXPRESSION (when logFC>0) or UNDEREXPRESSION,EXPRESSION (when logFC<0)

            ## Iterate available expression tags for the given gene and attempt match to a CIVIC record for that gene
            ## Opposite to the SNV and CNV case, here there is no tier hierarchy, ie. gene expression is either matched in CIVIC or not
            matched = []
            civicNames = list(varMap[gene].keys())
            for expr_tag in input_tags:
                if expr_tag in civicNames:
                    if expr_tag not in matched:
                        matched.append(expr_tag)

            ## Special case for EXPRESSION records related to EXONS: manually add them as matches when necessary
            ## Specific expression change must be identical in order to consider them as matched: eg. 'OVEREXPRESSION' and 'EXON 18 OVEREXPRESSION'
            for this_string in civicNames:
                (isExon,exprType) = expr_is_exon_string(this_string)
                if isExon and exprType:
                    # TODO: Records of the type 'P16 EXPRESSION' would be missed
                    if exprType in matched:
                        if this_string not in matched:
                            matched.append(this_string)

            ## At least 1 expression tag was matched in CIVIC
            if matched:
                nMatches += 1
            ## No match of expression tag in CIVIC for the current gene
            else:
                noMatches += 1
                ## Keep track of gene+expression that had CIVIC info available for the gene but could not be matched to specific expression
                gene_expr = gene + ":" + expr_change
                if gene_expr not in not_matched:
                    not_matched.append(gene_expr)

        ## If gene is not in CIVIC, keep track of unmatched gene
        else:
            notFound += 1
            if gene not in genesNotFound:
                genesNotFound.append(gene)
            matched = [] # in this case, list of CIVIC matches will be empty


    ### Process resulting matches in CIVIC (if any) for the current gene, and prepare info to be written to output. Always one gene per line only.
    geneScores = [] # to keep track of the CIVIC score for each matched expression record
    resultMap = {}  # to keep track of all evidence items retrieved from each matched CIVIC record
    therapyMap = {}    # to keep track of evidence for CIVIC therapy predictions for the given gene/line (across all matched records)

    for match in matched:
        matchData = varMap[gene][match]
        # Store the score for each individual CIVIC record for the current gene
        geneScores.append(match + ':' + str(matchData['civic_score']))

        # Process CIVIC evidence associated to the current record, and prepare results to be written to output
        # Important to be consistent with the column order in the outFile header
        for evidence_type in ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']:
            evidence_items = []
            if evidence_type not in resultMap.keys():
                resultMap[evidence_type] = []
            if evidence_type == 'PREDICTIVE':
                writetherapy = True
            else:
                writetherapy = False
            # Returns a list of evidence items for the given CIVIC record and evidence type, already in written form
            # Only evidences coming from the matched cancer specificity are returned (either 'ct' if in white list, 'gt' in in high level list or 'nct' if unspecific)
            (evidence_items,therapyMap) = match_cancer_and_return_evidences(matchData, evidence_type, cancerTypeList, blackList, highLevelList, therapyMap, writetherapy)
            # Add resulting evidence strings for the current record under the appropriate evidence type
            for i in evidence_items:
                # Add matched CIVIC record name to each evidence item string
                resultMap[evidence_type].append(match + ':' + i)


    ### Parse and process available therapy prediction evidence (if any) to agree on CIVIC support decision
    ### One support decision per available therapy; prioritize based on ct and apply majority vote for support decision

    ## Use dictionary of therapy support for the current gene/line
    ## Format: therapy -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
    supportStrings = []
    ## If no predictive evidence is available, corresponding column in output will be empty (ie. '.')
    for therapy in therapyMap.keys():
        ## Prioritize cancer specificity: ct > gt > nct
        thisCT = ''
        if 'ct' in therapyMap[therapy].keys():
            thisCT = 'ct'
        elif 'gt' in therapyMap[therapy].keys():
            thisCT = 'gt'
        elif 'nct' in therapyMap[therapy].keys():
            thisCT = 'nct'
        else:
            print("Error! Unexpected ct case for gene %s." %(gene))
            sys.exit(1)

        ## Given the selected ct, count number of occurrences for each possible support type (if any)
        count_pos = therapyMap[therapy][thisCT].count('POSITIVE')
        count_neg = therapyMap[therapy][thisCT].count('NEGATIVE')
        count_unk = therapyMap[therapy][thisCT].count('UNKNOWN_BLANK')
        count_dns = therapyMap[therapy][thisCT].count('UNKNOWN_DNS')

        ## Pool UNKNOWN_BLANK and UNKNOWN_DNS together (as both result in unknown CIVIC support)
        count_total_unk = count_unk + count_dns
        ## Sanity check that there is at least some support
        if (count_pos == 0) and (count_neg == 0) and (count_total_unk == 0):
            print("Error! Unexpected support case for gene %s." %(gene))
            sys.exit(1)

        ## Resolve contradicting evidence (if any) by majority vote
        tempSupport = ''
        ## For this, pool UNKNOWN_BLANK and UNKNOWN_DNS together
        ## Whenever there is a tie of "confident" (pos or neg) vs "non-confident" (unk), choose the confident one
        if (count_total_unk > count_pos) and (count_total_unk > count_neg):
            tempSupport = "CIVIC_UNKNOWN"
        elif count_pos == count_neg:
            tempSupport = "CIVIC_CONFLICT"
        elif (count_pos > count_neg) and (count_pos >= count_total_unk):
            tempSupport = "CIVIC_SUPPORT"
        elif (count_neg > count_pos) and (count_neg >= count_total_unk):
            tempSupport = "CIVIC_RESISTANCE"
        else:
            print("Error! Unexpected support case for gene %s." %(gene))
            sys.exit(1)

        ## Build support string for current therapy(+gene in line)
        ## Format: therapy:CT:SUPPORT
        therapiesupport = therapy + ':' + thisCT.upper() + ':' + tempSupport
        supportStrings.append(therapiesupport)


    ### Write results for current line to output (ie single gene + logFC value per line)

    ## Report CIVIC scores matched for the current gene's expression
    if geneScores:
        outfileLine += '\t' + ';'.join(geneScores)
    else:
        outfileLine += '\t.'
    ## Report CIVIC support decision for all the predicted therapies in the line
    ## Format: therapy1:CT:SUPPORT;therapy2:CT:SUPPORT;..;therapyN:CT:SUPPORT
    if supportStrings:
        outfileLine += '\t' + ';'.join(supportStrings)
    else:
        outfileLine += '\t.'
    ## Report CIVIC evidence associated to all relevant records in the current line
    ## Important to be consistent with the column order in the outFile header
    for evidence_type in ['PREDICTIVE', 'DIAGNOSTIC', 'PROGNOSTIC', 'PREDISPOSING']:
        isEmpty = True
        if evidence_type in resultMap.keys():
            if resultMap[evidence_type]:
                isEmpty = False
                outfileLine += '\t' + ';'.join(resultMap[evidence_type])
        # Check for cases when a field is empty
        if isEmpty:
            outfileLine += '\t.'

    ## Write output line for current gene to output
    outfile.write(outfileLine + '\n')

infile.close()
outfile.close()

print('\nTotal # matches: {}'.format(nMatches))
print('Total # no matches: {}'.format(noMatches))
print('Total # unavailable: {}'.format(notFound))
print('---------------------')
# TODO: only report # of genes instead of whole list?
print('Genes with no CIViC variant data:\n {}'.format(','.join(genesNotFound)))
print('\nUnmatched expression:\n {}'.format('\n '.join(not_matched)))
