#TumorProfiler, copy all files specified in suffix file lists
# prepare result upload
# create subfolders, copy files, copy and adapt readme, create the files [sample].raw_files.txt, [sample]__Summary.txt (in vsrs), and [sample].subClone.txt (in vsrs and derived)
# create md5suns
# create file listing the missing files (missing for upload)

set -e

outDir=$1
derivedList=$2
fastqDir=$3
sampleName=$4
cellranger_dir=$5
de_dir=$6
fastqc_dir=$7
enriched_dir=$8
celltype_dir=${9}
atypical_dir=${10}
rawH5_dir=${11}
cluster_dir=${12}
mean_dir=${13}
filtered_dir=${14}
plotting_dir=${15}
plotexpr_dir=${16}
gsva_dir=${17}
sequencing_runName=${18}
noTumor_orMaligOnly=${19} # either no_tumor or malig_only, defines what missing files are created

# create sub directories
mkdir -p ${outDir}"raw/"
mkdir -p ${outDir}"derived/"

# initialize helper files
touch ${outDir}raw/${sequencing_runName}__raw_files.txt
touch ${outDir}derived/${sampleName}__Summary.txt
touch ${outDir}derived/${sampleName}__fastqc_files.txt

# read in suffix lists, directly remove all special suffices  ; note: mac ^M removal with sed? sed 's/^M//g'
suffix_derived=$(cat ${derivedList} | grep -v "subClone" | grep -v "fastqc_files.txt" | grep -v "raw.h5" | grep -v "_summary." | grep -v "__Summary.txt" | grep -v "heatmap_enrichment.png" | grep -v "association.txt" )
#echo "Derived: ${suffix_derived}"

# fill subdirectory "raw"

for fastq in ${fastqDir}*.fastq.gz
do
	fastqName=$(basename ${fastq})
	cp ${fastq} ${outDir}raw/${sequencing_runName}__${fastqName}
	echo ${sequencing_runName}__${fastqName} >> ${outDir}raw/${sequencing_runName}__raw_files.txt
done

# fill subdirectory "derived"
# check what missing files need to be listed
if [ ${noTumor_orMaligOnly} == "no_tumor" ]; then
	echo ${sampleName}__subClone_diffExp_files.txt > ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__subClone_enrichment_files.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__subClone_clinicalAnno_files.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__heatmap_enrichment.png >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__DEmalignant.heatmap_enrichment.png >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugCombination.filteredDrugs.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugCombination.allDrugs.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugToCluster.filteredDrugs.vennplot.png >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugToCluster.allDrugs.vennplot.png >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugToCluster.allDrugs.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__full_druglist_to_subclones.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drug_prediction_umap.png >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__clusters_cell_count_percent.txt >> ${outDir}derived/${sampleName}__missing_files.txt
elif [ ${noTumor_orMaligOnly} == "malig_only" ]; then
	cp ${enriched_dir}/vs_other_malignant/${sampleName}.genes_cells_filtered.corrected.atypical_removed.DEmalignant.heatmap_enrichment.png ${outDir}derived/${sampleName}__DEmalignant.heatmap_enrichment.png
	echo ${sampleName}__subClone_diffExp_files.txt > ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__subClone_enrichment_files.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__subClone_clinicalAnno_files.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__heatmap_enrichment.png >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugCombination.filteredDrugs.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugCombination.allDrugs.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugToCluster.filteredDrugs.vennplot.png >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugToCluster.allDrugs.vennplot.png >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drugToCluster.allDrugs.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__full_druglist_to_subclones.txt >> ${outDir}derived/${sampleName}__missing_files.txt
	echo ${sampleName}__drug_prediction_umap.png >> ${outDir}derived/${sampleName}__missing_files.txt
fi

# fastqc files
for fastqc in ${fastqc_dir}*.html
do
	fastqcName=$(basename ${fastqc})
	cp ${fastqc} ${outDir}derived/${sampleName}__${fastqcName}
	echo ${sampleName}__${fastqcName} >> ${outDir}derived/${sampleName}__fastqc_files.txt
done

# h5 file raw

cp ${rawH5_dir}${sampleName}.h5 ${outDir}derived/${sampleName}__raw.h5

# stats files from cellranger

cp ${cellranger_dir}${sampleName}.metrics_summary.csv ${outDir}derived/${sampleName}__metrics_summary.csv
cp ${cellranger_dir}${sampleName}.web_summary.html ${outDir}derived/${sampleName}__web_summary.html

# phenograph_celltype_association.txt overview table

cp ${atypical_dir}${sampleName}.genes_cells_filtered.corrected.atypical_removed.phenograph_celltype_association.txt ${outDir}derived/${sampleName}__phenograph_celltype_association.txt

# count all remaining copied files to check if their number matches files from suffix lists

counter_derived=0

# double for-loop, maybe there is a better way
# celltype classification
echo "celltypes"
for celltype_file in ${celltype_dir}*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${celltype_file} == *${derivedSuffix} ]]; then
			echo ${celltype_file}
			cp ${celltype_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done


echo "cluster"
# clustering
for cluster_file in ${cluster_dir}*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${cluster_file} == *${derivedSuffix} ]]; then
			echo ${cluster_file}
			cp ${cluster_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done

echo "mean"
# mean expression
for mean_file in ${mean_dir}*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${mean_file} == *${derivedSuffix} ]]; then
			echo ${mean_file}
			cp ${mean_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done

echo "qc plots of filtered matrix"
# qc plots of filtered matrix
for filtered_qc_file in ${filtered_dir}*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${filtered_qc_file} == *${derivedSuffix} ]]; then
			echo ${filtered_qc_file}
			cp ${filtered_qc_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done

echo "plots sample composition"
# plots about the sample composition
for plot_file in ${plotting_dir}*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${plot_file} = *${derivedSuffix} ]]; then
			echo ${plot_file}
			cp ${plot_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done

echo "plots QC"
# plots about the sample composition
for qcPlot_file in ${plotting_dir}/qc_visualization/*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${qcPlot_file} = *${derivedSuffix} ]]; then
			echo ${qcPlot_file}
			cp ${qcPlot_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done

echo "plots gene expression"
# plots showing the gene expression
for plotexpr_file in ${plotexpr_dir}*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${plotexpr_file} = *${derivedSuffix} ]]; then
			echo ${plotexpr_file}
			cp ${plotexpr_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done

echo "atypical removed RDS"
# overview table over phenograph clusters and their cell type composition
for rds_file in ${atypical_dir}*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${rds_file} = *${derivedSuffix} ]]; then
			echo ${rds_file}
			cp ${rds_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done

echo "gsva plots"
# plots of gsva analysis
for gsvaplot_file in ${gsva_dir}*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${gsvaplot_file} = *${derivedSuffix} ]]; then
			echo ${gsvaplot_file}
			cp ${gsvaplot_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done

echo "bam files"
# bam file generated by cellranger in the very first analysis step
# get short sample name
[[ ${sampleName} =~ (.+)_scR_.+ ]]
short_sampleName=${BASH_REMATCH[1]}

# differenciate between samples sequenced on novaseq or nextseq
# if sample was sequenced on novaseq bam file is in preprocessing directory
if [[ ${fastqDir} =~ trial_novaseq_preprocessing ]] ; then
	[[ ${fastqDir} =~ (.+)/(fastqs|openbis)/.+ ]]
	root_bamdir=$(echo ${BASH_REMATCH[1]})
	bamdir=$(echo ${root_bamdir}/analysis/cellranger_run_index_hopping/${short_sampleName}/outs/)
else
# if sample was sequenced with nextseq bam file is in directory `cellranger_run` in analysis
	bamdir=$(echo ${cellranger_dir}${short_sampleName}/outs/)
# in bamdir search for bamfiles
fi

for bam_file in ${bamdir}*
do
	for derivedSuffix in ${suffix_derived}
	do
		if [[ ${bam_file} = *${derivedSuffix} ]]; then
			echo ${bam_file}
			cp ${bam_file} ${outDir}derived/${sampleName}__${derivedSuffix}
			counter_derived=$((counter_derived + 1))
		fi
	done
done


printf "Number of derived suffices: "
echo "${suffix_derived}" | wc -l

echo "Copied derived files: " $counter_derived

# change to directory and create md5 checksums
cd ${outDir}raw/

for f_raw in *
do
	outName=${f_raw}.md5
	md5sum ${f_raw} > ${outName}
done

cd ${outDir}derived/

for f_derived in *
do
	outName=${f_derived}.md5
	md5sum ${f_derived} > ${outName}
done
