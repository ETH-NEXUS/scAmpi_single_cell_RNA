#!/usr/bin/env python

'''
API call to CIViC DB for expression data
Requires Conda environment: async_API

Lourdes Rosano, Jul 2019
'''

## Execute with environment 'async_API'

import sys
import argparse
import asyncio
import aiohttp
import re
import time
# from itertools import permutations


### Define what combination of directions and clinical significances define 'POSITIVE' and 'NEGATIVE' support
supportDict = {'SUPPORTS': {'SENSITIVITY/RESPONSE':'POSITIVE', 'RESISTANCE':'NEGATIVE', 'REDUCED SENSITIVITY':'NEGATIVE', 'ADVERSE RESPONSE':'NEGATIVE'}, 'DOES NOT SUPPORT': {'RESISTANCE':'UNKNOWN_DNS', 'SENSITIVITY/RESPONSE':'UNKNOWN_DNS', 'REDUCED SENSITIVITY':'UNKNOWN_DNS', 'ADVERSE RESPONSE':'UNKNOWN_DNS'}}


'''
Functions
'''

### Functions for asynchronous API calls
### Adapted from 'https://pawelmhm.github.io/asyncio/python/aiohttp/2016/04/22/asyncio-aiohttp.html'

async def fetch(url, session):
    async with session.get(url) as response:
        # TODO: allow for more status codes (only 4XX and 5xx raise error)
        if response.status != 200:
            print("Error! Query {} failed (Status code: {})".format(url, response.status))
            sys.exit(1)
        return await response.json()

async def bound_fetch(sem, url, session):
    # Getter function with semaphore.
    async with sem:
        return await fetch(url, session)


## Query a list of elements (ids or id batches, ie. several ids per call allowed) using urlbase
## isBatch indicates whether elements correspond to single ids or batches of ids
async def run(elements, urlbase, isBatch=False):
    tasks = []

    # Create instance of Semaphore to limit number of concurrent requests
    sem = asyncio.Semaphore(32)

    # Fetch all responses within one Client session, keep connection alive for all requests
    async with aiohttp.ClientSession(headers={"Connection": "close"}) as session:
        for e in elements:
            if isBatch:
                # Join all IDs from the same batch for the corresponding call
                idString = ','.join(e)
            else:
                idString = e
            task = asyncio.ensure_future(bound_fetch(sem, urlbase.format(idString), session))
            tasks.append(task)

        responses = await asyncio.gather(*tasks, return_exceptions=True)
        return responses


## Execute asynchronous API calls to query CIVIC (can be used for both genes and variants)
## Included into a function in order to make error-handling work
def query_civic(elements, urlbase, isBatch=False):
    old_loop = asyncio.get_event_loop()
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    future = asyncio.ensure_future(run(elements, urlbase, isBatch))
    results = loop.run_until_complete(future)
    loop.close()
    asyncio.set_event_loop(old_loop)
    return results



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
# i.e. DISEASE[|DRUG1,DRUG2..](direction, significance(level(PMID,..,PMID),level(..)));
## For 'predictive' evidence (writeDrug=True), keep dictionary of drug support for the current gene/line
## Format: drug -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
def write_evidences(item, cancer, ct, drugMap, writeDrug=False):
    evidences = []
    # For each drug associated with the given cancer type
    for drug in item.keys():
        # For each evidence associated with the given drug
        # Evidences are simplified by using the combined form 'direction:significance'
        for evidence in item[drug].keys():
            ## For each evidence (ie combination of direction+clin_signf), count how many different evidence items support it
            ## At this stage, we find count evidence items by counting how many different combinations of level+pmids there are for the same drug, disease and evidence
            pmids = []
            # If drug=True, write drug information, i.e. DISEASE|DRUG(..)
            if writeDrug:
                # Always one drug (single or combination with '+')
                outString = cancer + '|' + drug + '('
            else:
                outString = cancer + '('
            # Split the evidence direction and clinical significance
            direction, clin_signf = evidence.split(':')
            outString += direction + ',' + clin_signf + '('
            # There may be several levels grouped per evidence
            levels = []
            for level in item[drug][evidence].keys():
                # There may be several publications (i.e. PMIDs) grouped per level
                levels.append(level + '(' + ','.join(item[drug][evidence][level]) + ')')
                # Count how many different evidence items support this particular evidence item
                for z in item[drug][evidence][level]:
                    # Distinguish cases where the same publication is used to support different and identical evidence levels (they nonetheless count as separate evidence items)
                    pmids.append(z)
#                     new_z = level + '_' + z
#                     if new_z not in pmids:
#                         pmids.append(new_z)
            outString += ','.join(levels) + '))'
            evidences.append(outString)

            ## For 'predictive' evidence, keep dictionary of drug support for the current gene/line
            ## Format: drug -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
            if writeDrug:
                ## Only add current drug and ct to dictionary if this match should indeed be considered for interpretation of the drug support
                if drug not in drugMap.keys():
                    drugMap[drug] = {}
                ## Possible values for cancer specificity: 'ct' (specific), 'gt' (general), 'nct' (non-specific)
                if ct not in drugMap[drug].keys():
                    drugMap[drug][ct] = []

                ## Evaluate drug support based on the combination of direction + clinical significance
                ## Each combination has an associated support: POSITIVE, NEGATIVE, UNKNOWN_DNS or UNKNOWN_BLANK
                if ('NULL' in direction) or ('N/A' in direction) or ('NULL' in clin_signf) or ('N/A' in clin_signf):
                    thisSupport = 'UNKNOWN_BLANK'
                else:
                    if direction not in supportDict.keys():
                        print("Error! Could not find direction %s in support dictionary." %(direction))
                        sys.exit(1)
                    if clin_signf not in supportDict[direction].keys():
                        print("Error! Could not find clinical significance %s in support dictionary." %(clin_signf))
                        sys.exit(1)
                    thisSupport = supportDict[direction][clin_signf]

                ## Keep track of number of occurrences for each support type for the given drug
                ## Here, take into account the number of supporting PMIDs associated to each evidence item
                for z in pmids:
                    drugMap[drug][ct].append(thisSupport)

    return (evidences,drugMap)


# Return in a list all evidence items (in written form) for a given evidence type 
# i.e. DISEASE1[|DRUG1,DRUG2..](direction1,significance1(level1(PMID,..,PMID),level2(..)));
#      DISEASE1[|DRUG1,DRUG2..](direction2,significance2(level1(PMID,..,PMID),level2(..)));
# For each evidence type, return either:
#   - Info on white listed cancer types (eg. 'Melanoma')
#   - If previous is not available, info on high level cancer types (eg. 'Cancer', 'Solid tumor')
#   - If previous is not available, info on all available cancer types for the given variant (except those included in the black list)
## For 'predictive' evidence (writeDrug=True), keep dictionary of drug support for the current gene/line
## Format: drug -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
def match_cancer_and_return_evidences(vardata, evidence_type, cancerTypes, cancerTypes_notSpec, highLevelTypes, drugMap, writeDrug=False):
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
            ## Also, return dictionary of drug support for the current gene/line
            (strings,drugMap) = write_evidences(vardata[evidence_type][cancerType], cancerType, ct, drugMap, writeDrug)
            for s in strings:
                evidences.append(s)

    return (evidences,drugMap)



'''
Script
'''

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
outHeader += "\tCIViC_Score\tCIViC_Drug_Support\tCIViC_Predictive\tCIViC_Diagnostic\tCIViC_Prognostic\tCIViC_Predisposing"
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
skipGenes = []
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

    ## Skip genes with logFC='NA' (ie could not asses differential expression), since empty values '.' will be assigned in CIVIC columns anyway
    if logfc_annot == 'NA':
        if gene not in naGenes:
            naGenes.append(gene)
        continue
    ## Skip exception cases like 'BX255925.3' or 'AL929554.17'
    if '.' in gene:
        if gene not in skipGenes:
            skipGenes.append(gene)
        continue
    if gene not in genes:
        genes.append(gene)
infile.close()

# When no genes are found, write empty files (necessary for snakemake pipeline) and exit without error
# eg. empty input file containing only header because patient had no DE genes at all
if (not genes) and (not skipGenes) and (not naGenes):
    print("\nDid not find any genes in column \'{}\' of file {}.".format(args.colName_gene, args.inputFile))
    outfile.close()
    sys.exit(0)
# When no genes were found because all existing were skipped (either logFC='NA' or symbol containing dots '.'),
# write output identical to input file plus additional CIVIC columns being empty
elif (not genes) and (skipGenes or naGenes):
    print("\nCould not find any genes that allow query to CIVIC (either all with logFC=NA or non-allowed gene symbols).")
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


### Asynchronous batch queries to CIViC (genes)

## There seems to be a maximum request size for the gene batch query
## Generate gene batches of smaller size to be able to query successfully
batchSize = 200
geneBatches = [genes[x:x+batchSize] for x in range(0, len(genes), batchSize)]

urlGenes = 'https://civic.nexus.ethz.ch/api/genes/{}?identifier_type=entrez_symbol'
# urlGenes = 'https://civicdb.org/api/genes/{}?identifier_type=entrez_symbol'
print('\nTotal # genes to query: {}'.format(len(genes)))
if skipGenes:
    print('\nOmitted gene symbols: {}\n'.format(','.join(list(skipGenes))))
if naGenes:
    print('\nOmitted genes due to unavailable logFC: {}\n'.format(','.join(list(naGenes))))
print('Retrieving gene data from CIViC...')
## Query CIVIC for the genes of interest
start = time.time()
results = query_civic(elements=geneBatches, urlbase=urlGenes, isBatch=True)
end = time.time()
print('  Response time: {}'.format(end-start))

## Error handling: when an exception occurs for a given query, it will be returned within the results
for r in results:
    if isinstance(r,Exception):
        print("Something went wrong! Exception: %s" %(r))
        sys.exit(1)

### Gather all variants for retrieved genes
print('\nParsing results...')
allvariants = []
retrieved = []
# Not all gene records will contain variants
for lb in results:
    # Sanity check for queries containing one single gene
    # Instead of a nested list (batches and genes), will return one single list (gene)
    if isinstance(lb, dict):
        lb = [lb]
    for dg in lb:
        # Only perfect matches are allowed (even gene aliases will not give a match)
        if dg['name'] not in retrieved:
            retrieved.append(dg['name'])
        for dv in dg['variants']:
            if dv['id'] not in allvariants:
                allvariants.append(dv['id'])

retrieved = set(retrieved)
unmatched = list(set(genes) - retrieved)
# TODO: report only total # of genes instead of whole list? Specially for CNVs this list becomes huge
print('  Genes Matched: {}'.format(','.join(list(retrieved))))
print('\n  Genes Not Matched: {}'.format(','.join(unmatched)))
print('\nTotal # variants to query: {}'.format(len(allvariants)))


### Asynchronous queries to CIViC (variants)

urlVars = 'https://civic.nexus.ethz.ch/api/variants/{}'
# urlVars = 'https://civicdb.org/api/variants/{}'
print('Retrieving variant data from CIViC...')
## Query CIVIC for the variants of interest
start = time.time()
results = query_civic(elements=allvariants, urlbase=urlVars, isBatch=False)
end = time.time()
print('  Response time: {}'.format(end-start))

## Error handling: when an exception occurs for a given query, it will be returned within the results
for r in results:
    if isinstance(r,Exception):
        print("Something went wrong! Exception: %s" %(r))
        sys.exit(1)

### Parse results and create dictionary of gene-variant information
print('\nParsing results...')
varMap = {}
for dv in results:
    # Skip variants which do not have any clinical data associated to them
    if not dv['evidence_items']:
#     if not (dv['evidence_items'] or dv['assertions']):
        continue
    # Iterate through the evidence items and store relevant information
    for item in dv['evidence_items']:
        # Sanity check that all critical elements are present and non-empty
        if not (dv['entrez_name'] and dv['name'] and item['evidence_type'] and item['disease']['name']):
            continue
        # Skip records that are not accepted evidence
        if item['status'].upper() != 'ACCEPTED':
            continue

        # Gene names in CIVIC are HUGO symbols (uppercase) but do a sanity check nevertheless
        gene = dv['entrez_name'].upper()
        variant = dv['name'].upper()
        evidenceType = item['evidence_type'].upper()
        cancerType = item['disease']['name'].upper()
        # Sanity check for empty evidence direction, clinical significance or level
        # 'NULL' is introduced to distinguish from 'N/A' tag
        if item['evidence_direction'] is None:
            item['evidence_direction'] = 'NULL'
        if item['clinical_significance'] is None:
            item['clinical_significance'] = 'NULL'
        if item['evidence_level'] is None:
            item['evidence_level'] = 'NULL'
        # Combine the direction and significance of the evidence in one term
        evidence = item['evidence_direction'].upper() + ':' + item['clinical_significance'].upper()
        level = item['evidence_level'].upper()
        if gene not in varMap.keys():
            varMap[gene] = {}
        # Variant name should be unique within gene
        # (found some duplicates but all were submitted, not accepted data)
        if variant not in varMap[gene].keys():
            varMap[gene][variant] = {}
            # Internal CIViC ID
            varMap[gene][variant]['id'] = dv['id']
            # Score to assess the accumulation of evidence for each variant (quantity and quality)
            # Sanity check for empty scores
            if dv['civic_actionability_score'] is not None:
                varMap[gene][variant]['civic_score'] = dv['civic_actionability_score']
            else:
                varMap[gene][variant]['civic_score'] = 'NULL'

            # Keep original HGVS annotations (empty list when nothing is available)
            # Use uppercase to avoid mismatches due to case
            varMap[gene][variant]['hgvs'] = [h.upper() for h in dv['hgvs_expressions']]

        # TODO: there is no sanity check for detecting possible variant name duplicates
        if evidenceType not in varMap[gene][variant].keys():
            varMap[gene][variant][evidenceType] = {}
        if cancerType not in varMap[gene][variant][evidenceType].keys():
            varMap[gene][variant][evidenceType][cancerType] = {}
        if item['drugs']:
            drugs = [d['name'].upper() for d in item['drugs']]
            # When more than 1 drug are listed for the same evidence item, 'drug_interaction_type' is not null and defines the nature of this multiple drug entry 
            if item['drug_interaction_type'] is not None:
                # 'Substitutes' indicates that drugs can be considered individually
                if item['drug_interaction_type'].upper() != 'SUBSTITUTES':
                    # Remaining terms ('Sequential' and 'Combination') indicate that drugs should be considered together, so join their names into a single tag
                    # Sort drugs alphabetically to ensure that their order in the combination treatment is always the same
                    drugs.sort()
                    drugs = ['+'.join(drugs)]
#                     # Consider all possible permutations of the drug list
#                     drugs = ['+'.join(per) for per in permutations(drugs)]
#                     drugMatch = None
#                     for drug in drugs:
#                         if drug in varMap[gene][variant][evidenceType][cancerType].keys():
#                             drugMatch = drug
#                             break
#                     # If drug combination is new, get the first permutation
#                     if drugMatch is None:
#                         drugs = [drugs[0]]
#                     else:
#                         drugs = [drugMatch]
        else:
            # Only non-Predictive evidences and Predictive ones without drugs will have this dummy level
            # Introduced for consistency purposes within the varMap structure
            drugs = ['NULL']

        # Iterate through drugs to add evidences associated to them
        #   For non-Predictive evidences or Predictive with empty drugs, drugs=['NULL']
        #   For Predictive and interaction=None, len(drugs) = 1
        #   For Predictive and interaction='Substitutes', len(drugs)>1
        #   For Predictive and interaction!='Substitutes', len(drugs)=1 (combiantion of several using '+')
        for drug in drugs:
            if drug not in varMap[gene][variant][evidenceType][cancerType].keys():
                varMap[gene][variant][evidenceType][cancerType][drug] = {}
            if evidence not in varMap[gene][variant][evidenceType][cancerType][drug].keys():
                varMap[gene][variant][evidenceType][cancerType][drug][evidence] = {}
            if level not in varMap[gene][variant][evidenceType][cancerType][drug][evidence].keys():
                varMap[gene][variant][evidenceType][cancerType][drug][evidence][level] = []
            # Group all publications associated to the same level. Do not check publication status
            ## On 25.01.2019, source structure was changed to introduce ASCO abstracts as a source type
            ## TODO: sanity check for empty ID. Check for type of source?
            varMap[gene][variant][evidenceType][cancerType][drug][evidence][level].append(item['source']['citation_id'])
#             varMap[gene][variant][evidenceType][cancerType][drug][evidence][level].append(item['source']['pubmed_id'])

    # TODO: iterate through assertions and repeat above process
    # for item in dv['assertions']:


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
                    ## TODO: Records of the type 'P16 EXPRESSION' would be missed
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
    drugMap = {}    # to keep track of evidence for CIVIC drug predictions for the given gene/line (across all matched records)

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
                writeDrug = True
            else:
                writeDrug = False
            # Returns a list of evidence items for the given CIVIC record and evidence type, already in written form
            # Only evidences coming from the matched cancer specificity are returned (either 'ct' if in white list, 'gt' in in high level list or 'nct' if unspecific)
            (evidence_items,drugMap) = match_cancer_and_return_evidences(matchData, evidence_type, cancerTypeList, blackList, highLevelList, drugMap, writeDrug)
            # Add resulting evidence strings for the current record under the appropriate evidence type
            for i in evidence_items:
                # Add matched CIVIC record name to each evidence item string
                resultMap[evidence_type].append(match + ':' + i)


    ### Parse and process available drug prediction evidence (if any) to agree on CIVIC support decision
    ### One support decision per available drug; prioritize based on ct and apply majority vote for support decision

    ## Use dictionary of drug support for the current gene/line
    ## Format: drug -> ct -> [support1,support2,..,support1] (keep track of all occurrences)
    supportStrings = []
    ## If no predictive evidence is available, corresponding column in output will be empty (ie. '.')
    for drug in drugMap.keys():
        ## Prioritize cancer specificity: ct > gt > nct
        thisCT = ''
        if 'ct' in drugMap[drug].keys():
            thisCT = 'ct'
        elif 'gt' in drugMap[drug].keys():
            thisCT = 'gt'
        elif 'nct' in drugMap[drug].keys():
            thisCT = 'nct'
        else:
            print("Error! Unexpected ct case for gene %s." %(gene))
            sys.exit(1)

        ## Given the selected ct, count number of occurrences for each possible support type (if any)
        count_pos = drugMap[drug][thisCT].count('POSITIVE')
        count_neg = drugMap[drug][thisCT].count('NEGATIVE')
        count_unk = drugMap[drug][thisCT].count('UNKNOWN_BLANK')
        count_dns = drugMap[drug][thisCT].count('UNKNOWN_DNS')

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

        ## Build support string for current drug(+gene in line)
        ## Format: DRUG:CT:SUPPORT
        drugSupport = drug + ':' + thisCT.upper() + ':' + tempSupport
        supportStrings.append(drugSupport)


    ### Write results for current line to output (ie single gene + logFC value per line)

    ## Report CIVIC scores matched for the current gene's expression
    if geneScores:
        outfileLine += '\t' + ';'.join(geneScores)
    else:
        outfileLine += '\t.'
    ## Report CIVIC support decision for all the predicted drugs in the line
    ## Format: DRUG1:CT:SUPPORT;DRUG2:CT:SUPPORT;..;DRUGN:CT:SUPPORT
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
