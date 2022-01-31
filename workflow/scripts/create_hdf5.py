import argparse
from utils import cellranger_to_hdf5

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--genes_file", required=True, help="tsv file containing gene names")
parser.add_argument("-m", "--matrix_file", required=True, help="file containing the geneXcell matrix")
parser.add_argument("-b", "--barcodes_file", required=True, help="file containing the cell ids(barcodes)")
parser.add_argument("-o", "--output_file", required=True, help="path to the output hdf5 file")
parser.add_argument("--n_threads", help="the number of threads", type=int,
        default=1)
args = parser.parse_args()



cellranger_to_hdf5(args.genes_file, args.matrix_file, args.barcodes_file,
        args.output_file)
