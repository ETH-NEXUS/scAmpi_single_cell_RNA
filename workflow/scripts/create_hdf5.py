import argparse
import numpy as np
from scipy.io import mmread
import h5py
import os

def cellranger_to_hdf5 (genes, matrix_file, barcodes, out_path):

    # this is needed for the filesystems that do not allow file locking
    os.environ['HDF5_USE_FILE_LOCKING'] = "FALSE"
    # to read a matrix market file
    matrix = mmread(matrix_file.__str__()).astype("float32").todense().T
    matrix = np.asarray(matrix)

    gene_ids = np.genfromtxt(genes.__str__(), dtype='S16')[:,0]
    gene_names = np.genfromtxt(genes.__str__(), dtype='S16')[:, 1]
    cell_names = np.genfromtxt(barcodes.__str__(), dtype='S16')

    # removing all-zero-genes accross all cells
    detected_genes_index = ~(matrix == 0).all(axis=0)
    matrix = matrix[:,detected_genes_index]
    gene_ids = gene_ids[detected_genes_index]
    gene_names = gene_names[detected_genes_index]


    f = h5py.File(out_path.__str__(), "w")
    f.create_dataset(name = 'raw_counts', data = matrix)
    gg = f.create_group('gene_attrs')
    gg.create_dataset(name = 'gene_names', data = gene_names)
    gg.create_dataset(name = 'gene_ids', data = gene_ids)
    cg = f.create_group('cell_attrs')
    cg.create_dataset(name = 'cell_names', data = cell_names)
    cg.create_dataset(name = 'cells_on_rows', data = True)

    f.close()

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
