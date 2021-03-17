snakemake -s snake_example.snake --configfile /config_clinical.json --rulegraph > rulegraph_TP_master.dot

module load gcc/4.8.2 graphviz/2.38

dot -Tpng rulegraph_TP_master.dot > rulegraph_TP_master.png
