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

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
