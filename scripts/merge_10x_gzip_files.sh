#!/bin/bash

# Merges gzip files which match a pattern to one big gzip file and removes the 
# original files

#pattern1=$1/'*L001_R1*'
#pattern2=$1/'*L001_R2*'
#pattern3=$1/'*L001_I1*'

#pattern1=$1/'*L002_R1*'
#pattern2=$1/'*L002_R2*'
#pattern3=$1/'*L002_I1*'

#pattern1=$1/'*L003_R1*'
#pattern2=$1/'*L003_R2*'
#pattern3=$1/'*L003_I1*'

#pattern1=$1/'*L004_R1*'
#pattern2=$1/'*L004_R2*'
#pattern3=$1/'*L004_I1*'

function merge {
        echo Processing $1
	first_file=$(echo $1 | cut -f1 -d " " )
        basename=`basename $first_file` 
        path=`dirname $first_file`
        #echo $path
        #echo $basename
	newfile=`echo $basename | sed  's/BSSE/MERGED_BSSE/'`
        echo -e "New file will be called: $newfile\n"

	cat $@ > $path/$newfile
        #rm $@
}

for i in {1..4}
do 
	pattern1=$1/*L00${i}_R1*
	pattern2=$1/*L00${i}_R2*
	pattern3=$1/*L00${i}_I1*
	
	merge "$pattern1"
	merge "$pattern2"
	merge "$pattern3"
done
