#!/bin/bash

#################################################
#### File name: replace_oldHGNC_symbols.sh
#### Author: Anne Bertolini
#### Date created: November 2019
##################################################

# Replace old HGNC symbols in any tab delimited file with the current symbol.
# Takes as input:
# 1) the table with new and old symbols. The table should be output from the
#    HGNC multi-symbol checker. In column 1 there is the input (old) symbol and in
#    column 3 the approved symbol. Headers are "Input" and "Approved symbol".
#    The table should be filtered to contain only those lines that are supposed
#    to get changed.
# 2) the text file that will be changed so that it contains the up-to-date
#    symbols.

# read in output table of HGNC Multi-symbol checker
inChangedSymbols=${1}
# any text file that is getting updated.
inTextFile=${2}

SEARCH=( $(cat ${inChangedSymbols} | cut -f1 | grep -v Input ))
REPLACE=( $(cat ${inChangedSymbols} | cut -f3 | grep -v Approved))

# count total number of replacements
SUM_REPLACED=0

# Give out length of array
echo "Genes to be replaced with new symbols: ${#SEARCH[@]}"
echo "New symbols: ${#REPLACE[@]}"

for ((i = 0; i < ${#SEARCH[@]}; ++i)) do
	position=$(( $i + 1 ))
	echo "${SEARCH[$i]} -> ${REPLACE[$i]}"
	# count exact (word-regexp) matches
	OCCU=( $(grep -wo ${SEARCH[$i]} ${inTextFile} | wc -l ) )
	echo "$OCCU ${SEARCH[$i]}"
	SUM_REPLACED=$(($SUM_REPLACED + $OCCU))
#	echo "$SUM_REPLACED : Sum of replacements"
	# replace only exact matches of ${SEARCH[$i]}
	sed -i "s/\<${SEARCH[$i]}\>/${REPLACE[$i]}/g" ${inTextFile}
	echo ""
	wait
done
echo "Total number of replacements: $SUM_REPLACED"
