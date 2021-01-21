/cluster/project/nexus/utilities/sharedPrograms/snakemake/snakemake_4.8.0/bin/snakemake -s /cluster/work/bewi/ngs/projects/tumorProfiler/code/clinical/roche_tumorprofiler_singlecell_2018/snake_example.snake --configfile /cluster/work/bewi/ngs/projects/tumorProfiler/code/clinical/roche_tumorprofiler_singlecell_2018/config_clinical.json --rulegraph > rulegraph_TP_master.dot

module load gcc/4.8.2 graphviz/2.38

dot -Tpng rulegraph_TP_master.dot > rulegraph_TP_master.png
