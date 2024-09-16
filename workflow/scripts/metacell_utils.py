# Author: Matthias Lienhard
# Date: 2024-08-12
# Description: Shared util functions for metacell scripts

import h5py
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix


def h5_to_ad(h5file, celltype_mapping_file=None):
    "import the h5 file"
    # Open the HDF5 file and read the counts data
    with h5py.File(h5file, 'r') as f:
        # Load counts data
        counts = csr_matrix(f['raw_counts'][:])
        # Load cell names
        cell_names = f['cell_attrs/cell_names'][:]
        cell_names = [name.decode('utf-8') for name in cell_names]
        # Load gene ids
        gene_ids = f['gene_attrs/gene_ids'][:]
        gene_ids = [gid.decode('utf-8') for gid in gene_ids]
        # Load gene names
        gene_names = f['gene_attrs/gene_names'][:]
        gene_names = [name.decode('utf-8') for name in gene_names]
    # Create AnnData object with sparse matrix representation
    adata = sc.AnnData(X=counts, obs=pd.DataFrame(
        index=cell_names), var=pd.DataFrame({'gene_names': gene_names}, index=gene_ids))
    if celltype_mapping_file is not None:
        celltype_mapping = dict(pd.read_csv(
            celltype_mapping_file, sep='\t').set_index('barcodes')['celltype_final'])
        if not set(celltype_mapping.keys()).issubset(set(adata.obs_names)):
            print("WARNING: Some cell barcodes in the celltype mapping file do not match the cell barcodes in the h5 input")
        bc_found = [bc in set(celltype_mapping.keys())
                    for bc in adata.obs_names]
        print(
            f'found celltype assignments for {sum(bc_found)}/{len(bc_found)} cells')
        # Map the cell types to the AnnData object
        adata.obs['celltype'] = adata.obs_names.map(celltype_mapping)
    return adata
