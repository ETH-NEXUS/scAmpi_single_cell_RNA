# Author: Matthias Lienhard
# Date: 2024-08-12
# Description: Running Metacells2 metacells identification

import h5py
import matplotlib.pyplot as plt  # For plotting
import metacells as mc           # The Metacells package
import numpy as np               # For array/matrix operations
import pandas as pd              # For data frames
import os                        # For filesystem operations
import seaborn as sns             # For plotting
from scipy.sparse import csr_matrix, vstack
from scanpy import AnnData
# from math import hypot           # For plotting
import argparse


# minimal set of genes to exclude
EXCLUDED_GENE_NAMES = [
    "XIST", "MALAT1",   # Sex-specific genes.
    "NEAT1"             # Non-coding.
]
EXCLUDED_GENE_PATTERNS = ["MT-.*"]  # Mytochondrial.
#minimal set of genes to maskd
LATERAL_NAMES = [
    "AURKA", "MCM3", "MCM4", "MCM7", "MKI67", "PCNA", "RRM2", "SMC4", "TPX2",  # Cell-cycle
    "FOS", "HSP90AB1", "TXN",                                                  # Stress
]
LATERAL_PATTERNS = ["RP[LS].*"] # Ribosomal


def main(h5file, out_dir, out_prefix,celltype_mapping_file,target_metacell_umis,target_metacell_size, exclude_names,exclude_pattern, lateral_names, lateral_pattern ):
    # read the data
    sample=out_prefix.rstrip('_.')
    cells=h5_to_ad(h5file, celltype_mapping_file=celltype_mapping_file)
    # prepare for metacells2
    mc.ut.top_level(cells)
    mc.ut.set_name(cells, sample)
    print(f"Imported dataset with {cells.n_obs} cells, {cells.n_vars} genes")
    # filtering
    mc.pl.exclude_genes(
        cells,
        excluded_gene_names=exclude_names, 
        excluded_gene_patterns=exclude_pattern,
        random_seed=42)
    # masking
    mc.pl.mark_lateral_genes(cells, lateral_gene_names=lateral_names, lateral_gene_patterns=lateral_pattern)
    # assigning cells to metacells
    mc.pl.divide_and_conquer_pipeline(
        cells, 
        target_metacell_umis=target_metacell_umis, 
        target_metacell_size=target_metacell_size,
        random_seed=42)


    metacells= mc.pl.collect_metacells(cells, name=f"{sample}.one-pass.metacells", random_seed=42)
    print(f'found {metacells.shape[0]} metacells')
    assign_df=get_assignment(cells)
    
    assign_df.to_csv(f'{out_dir}/{out_prefix}_metacells2_assignment.tsv', sep='\t')
    assign_df['UMI_counts'] =  cells.X.sum(axis=1).A1
    assign_df['fractionMT'] = cells.obs['fractionMT'].reindex(assign_df['cell_barcode']).values
    metacell_fractionMT = (assign_df.groupby('metacell1_id')[['fractionMT', 'UMI_counts']]
                            .apply(lambda x: np.average(x['fractionMT'], weights=x['UMI_counts'])))
    metacells.obs['fractionMT']=metacell_fractionMT


    # create metecells adata object
    h5_out=f'{out_dir}/{out_prefix}_metacells2.h5'
    #write to file
    export_metacells( metacells, h5_out)
    #h5_out=f'{out_dir}/{out_prefix}_seacells_soft.h5'
    #export_seacells(SEACell_ad, gene_dict, h5_out)

def aggregate_metacell_counts(assign_df, cells):
    barcode_df=dict(assign_df['metacell1_id'])
    cell_metacell_ids=np.array([barcode_df[bc] for bc in cells.obs.index])
    metacell_ids=sorted(set(cell_metacell_ids))
    metacell_ids.remove('Outliers')
    aggregated_counts=[csr_matrix(cells.X[cell_metacell_ids==mid].sum(0)) for mid in metacell_ids]
    aggregated_counts = vstack(aggregated_counts, format='csr')
    metacell_obs = pd.DataFrame( {'total_umis':aggregated_counts.sum(1).A1},index=metacell_ids)
    metacell_var=cells.var['gene_ids']
    aggregated_metacells=AnnData(X=aggregated_counts, obs=metacell_obs, var=metacell_var)
    return aggregated_metacells

def h5_to_ad(h5file, celltype_mapping_file=None):
    "import the h5 file"
    # Open the HDF5 file and read the counts data
    with h5py.File(h5file, 'r') as f:
        # Load counts data
        counts = csr_matrix(f['raw_counts'][:]).astype('float32')
        # Load cell names
        cell_names = f['cell_attrs/cell_names'][:]
        cell_names = [name.decode('utf-8') for name in cell_names]
        # Load gene ids
        gene_ids = f['gene_attrs/gene_ids'][:]
        gene_ids = [gid.decode('utf-8') for gid in gene_ids]
        # Load gene names
        gene_names = f['gene_attrs/gene_names'][:]
        gene_names = [name.decode('utf-8') for name in gene_names]
        # Load fractionMT
        fractionMT = f['cell_attrs/fractionMT'][:]
    # Create AnnData object with sparse matrix representation
    adata = AnnData(X=counts, obs=pd.DataFrame(index=cell_names), var=pd.DataFrame({'gene_ids':gene_ids},index=gene_names))
    adata.obs['fractionMT'] = fractionMT
    if celltype_mapping_file is not None:
        celltype_mapping = dict(pd.read_csv(celltype_mapping_file, sep='\t').set_index('barcodes')['celltype_final'])
        if not set(celltype_mapping.keys()).issubset(set(adata.obs_names)):
            print("WARNING: Some cell barcodes in the celltype mapping file do not match the cell barcodes in the h5 input")
        bc_found=[bc in set(celltype_mapping.keys()) for bc in adata.obs_names]
        print(f'found celltype assignments for {sum(bc_found)}/{len(bc_found)} cells')
        # Map the cell types to the AnnData object
        adata.obs['celltype'] = adata.obs_names.map(celltype_mapping)
    return adata

def get_assignment(cells):
    df=cells.obs[['metacell_name', 'metacell']]
    # Rename the columns
    df.rename(columns={'metacell_name': 'metacell1_id'}, inplace=True)
    # Add cell barcode column
    df['cell_barcode'] = df.index
    # Reorder columns to have 'cell_barcode' as the first column
    df = df[['cell_barcode', 'metacell1_id', 'metacell']]
    return(df)


def export_metacells(metacells, h5_out=None):
    "export metacells agregated counts to h5 file"
    if h5_out is not None:
        with h5py.File(h5_out, 'w') as f:
            # Save counts
            f.create_dataset('raw_counts', data=metacells.layers['total_umis'], dtype=int)
            # Save cell names (metacell names)
            metacell_names = np.array(metacells.obs_names, dtype='S')
            f.create_dataset('cell_attrs/cell_names', data=metacell_names)
            # Save gene ids
            gene_names = np.array(metacells.var_names, dtype='S')
            f.create_dataset('gene_attrs/gene_names', data=gene_names)
            # Save gene names
            gene_ids = metacells.var.gene_ids.values.astype('S')
            f.create_dataset('gene_attrs/gene_ids', data=gene_ids)
            # save the fractionMT
            fraction_mt=metacells.obs['fractionMT'].values
            f.create_dataset('cell_attrs/fractionMT', data=fraction_mt)
            # Indicate that cells are on rows
            f.create_dataset('cell_attrs/cells_on_rows', data=np.array(True))
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Running Metacells2 metacells identification")
    parser.add_argument("-i", metavar='in_file.h5', required=True, help="Input file with raw counts")
    parser.add_argument("-o", metavar='out_dir', required=True, help="Output directory")
    parser.add_argument("-p", metavar="out_prefix", required=True, help="Prefix for the output files")
    parser.add_argument("-c", metavar="celltype_mapping", help="Optional: File with barcode celltype mappings. Used for celltype purity of metacells.",default=None)
    parser.add_argument("--target_metacell_umis", type=int, default=100000, help="Target number of umis per metacell")
    parser.add_argument("--target_metacell_size", type=int, default=30, help="Target number of cells per metacell")
    parser.add_argument("--exclude_names", type=str, default=','.join(EXCLUDED_GENE_NAMES), help=", seperated list of genes to exclude (sex, non-coding, mt, ...)")
    parser.add_argument("--exclude_pattern", type=str, default=EXCLUDED_GENE_PATTERNS[0], help=", seperated list of pattern to exclude (sex, non-coding, mt, ...)")
    parser.add_argument("--lateral_names", type=str, default=','.join(EXCLUDED_GENE_NAMES), help=", seperated list of genes to mask as lateral (cell cycle, stress response, ribosomal, ...)")
    parser.add_argument("--lateral_pattern", type=str, default=EXCLUDED_GENE_PATTERNS[0], help=", seperated list of pattern to mask as lateral (cell cycle, stress response, ribosomal ...)")
    
    args = parser.parse_args()
    # parse the tags
    h5file=args.i
    out_dir=args.o
    out_prefix=args.p
    celltype_mapping_file=args.c
    target_umis = int(args.target_metacell_umis)
    target_size = int(args.target_metacell_size)
    exclude_names = args.exclude_names.split(',')
    exclude_pattern = args.exclude_pattern.split(',')
    lateral_names = args.lateral_names.split(',')
    lateral_pattern = args.lateral_pattern.split(',')
    print(f"input file: {h5file}")
    print(f"output files: {out_dir}/{out_prefix}*\n")
    if celltype_mapping_file is not None:
        print(f"celltype mapping: {celltype_mapping_file}")
    else:
        print("no celltype mapping provided")
    main(h5file, out_dir, out_prefix, celltype_mapping_file,target_umis,target_size, exclude_names,exclude_pattern, lateral_names, lateral_pattern)
    



