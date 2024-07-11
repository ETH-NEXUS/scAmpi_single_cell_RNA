import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
import SEACells
import matplotlib
matplotlib.use('Agg')  # Use Agg backend for rendering plots

import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import warnings

# parameters of the script:
# -i: h5 Input file with raw counts
# -o: Output directory
# -p: Prefix for the output files
# Optional:
# -c: File with barcode celltype mappings. Used for celltype purity and umap QC plots. [default: None]
# --n_cells: number of metacells [default: 90]
# --n_pc: pca dimensions [default: 50]
# --n_waypoint_eigs: Number of eigenvalues to consider when initializing metacells [default: 10]
# --n_topvar_genes: number of top variable genes [default: 1500]
# --n_neighbors: # neighbors considered for the neighborhood graph [default: 15]
# --iter: number of iterations to fit, provided as min,max [default: 10,50]


class SEACellsWarning(RuntimeWarning):
    # seacells raises them, so we need to catch and transform it to warning
    pass


# high level function of metacell identification with SEACells:
def main(h5file, out_dir, out_prefix, celltype_mapping_file, n_SEACells, n_pc, n_waypoint_eigs, n_topvar_genes, n_neighbors, min_iter, max_iter):
    # import the data
    adata=h5_to_ad(h5file, celltype_mapping_file=celltype_mapping_file)
    # select topvar genes, compute neighborhood graph and so on
    adata=prepare_seacells(adata, n_topvar_genes, n_pc, n_neighbors)
    # fit the actual model - this takes a minute or two
    model=fit_seacells_model(adata=adata, n_SEACells=n_SEACells, n_waypoint_eigs=n_waypoint_eigs, min_iter=min_iter, max_iter=max_iter)
    # do some QC plots:
    # Check for convergence 
    model.plot_convergence(save_as=f'{out_dir}/{out_prefix}_seacells_model_convergence.png') #

    if celltype_mapping_file is not None:
        # Umap plots with metacells on top
        sc.pl.scatter(adata, basis='umap', frameon=False, color='celltype')
        plt.savefig(f'{out_dir}/{out_prefix}_seacells_celltype_umap.png')
        # seacell purity barplot
        SEACell_purity = SEACells.evaluate.compute_celltype_purity(adata, 'celltype')
        plt.figure(figsize=(4,4))
        sns.boxplot(data=SEACell_purity, y='celltype_purity')
        plt.title('Celltype Purity')
        sns.despine()
        plt.savefig(f'{out_dir}/{out_prefix}_seacells_celltype_purity.png')
        plt.close()
    SEACells.plot.plot_2D(adata, key='X_umap', show=False, colour_metacells=False, save_as=f'{out_dir}/{out_prefix}_seacells_umap.png')
    SEACells.plot.plot_2D(adata, key='X_umap', show=False, colour_metacells=True, save_as=f'{out_dir}/{out_prefix}_seacells_umap_color.png')
    SEACells.plot.plot_SEACell_sizes(adata, bins=5, show=False, save_as=f'{out_dir}/{out_prefix}_seacells_sizes.png')
    # export the assignments
    hard_df=get_assignment(model, hard=True)
    hard_df.to_csv(f'{out_dir}/{out_prefix}_seacells_hard_assignment.tsv', sep='\t')
    soft_df=get_assignment(model, hard=False)
    soft_df.to_csv(f'{out_dir}/{out_prefix}_seacells_assignment.tsv', sep='\t')
    soft_df['UMI_counts'] =  adata.obs['n_counts']

    soft_df['fractionMT'] = adata.obs['fractionMT'].reindex(soft_df.index).values
    metacell_fractionMT = (soft_df.groupby('metacell1_id')[['fractionMT', 'UMI_counts']]
                            .apply(lambda x: np.average(x['fractionMT'], weights=x['UMI_counts'])))
    #  aggregate output matrix: raw counts for each metacell
    SEACell_ad = SEACells.core.summarize_by_SEACell(adata, SEACells_label='SEACell', summarize_layer='raw')
    SEACell_ad.obs['fractionMT']=metacell_fractionMT
    SEACell_ad=prepare_seacells(SEACell_ad,  n_topvar_genes, n_pc, n_neighbors)

    sc.pl.scatter(SEACell_ad, basis='umap', frameon=False, color='n_counts')
    plt.savefig(f'{out_dir}/{out_prefix}_seacells_umap_metacells_ncounts.png')

    # alternatively, "soft" (weighted) assignment
    SEACell_soft_ad = SEACells.core.summarize_by_soft_SEACell(adata, model.A_, summarize_layer='raw', minimum_weight=0.05)
    SEACell_soft_ad.obs.index=SEACell_ad.obs.index
    SEACell_soft_ad.obs['fractionMT']=metacell_fractionMT
    SEACell_soft_ad=prepare_seacells(SEACell_soft_ad,  n_topvar_genes, n_pc, n_neighbors)
    sc.pl.scatter(SEACell_soft_ad, basis='umap', frameon=False, color='n_counts')
    plt.savefig(f'{out_dir}/{out_prefix}_seacells_umap_metacells_soft_ncounts.png')

    # export to h5 files
    gene_dict=dict(adata.var["gene_names"]) # gene names are lost in the process and have to be added here
    h5_out=f'{out_dir}/{out_prefix}_seacells.h5'
    export_seacells( SEACell_ad, gene_dict, h5_out)
    # h5_out=f'{out_dir}/{out_prefix}_seacells_soft.h5'
    # export_seacells(SEACell_soft_ad, gene_dict, h5_out)

def h5_to_ad(h5file, celltype_mapping_file=None):
    "import the h5 file"
    # Open the HDF5 file and read the counts data
    with h5py.File(h5file, 'r') as f:
        # Load counts data
        counts = csr_matrix(f['raw_counts'][:], dtype=np.float32)
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
    adata = sc.AnnData(X=counts, obs=pd.DataFrame(index=cell_names), var=pd.DataFrame({'gene_names':gene_names},index=gene_ids))
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

def prepare_seacells(adata, n_topvar_genes, n_pc, n_neighbors):
    "basic scanpy preprocessing for umap"
    # copy the data to .raw
    raw_ad = sc.AnnData(adata.X, dtype=np.float32)
    raw_ad.obs_names, raw_ad.var_names = adata.obs_names, adata.var_names
    adata.raw = raw_ad 
    # Normalize cells, log transform and compute highly variable genes
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_topvar_genes)
    sc.tl.pca(adata, n_comps=n_pc, use_highly_variable=True)
    # Compute the neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pc)  # Use the number of PCs computed earlier
    # Compute UMAP embedding
    sc.tl.umap(adata)
    return adata

def fit_seacells_model(adata, n_SEACells, n_waypoint_eigs, min_iter=15, max_iter=50, build_kernel_on="X_pca"):
    "initialize and fit the model"        
    model = SEACells.core.SEACells(adata, 
                    build_kernel_on=build_kernel_on, 
                    n_SEACells=n_SEACells, 
                    n_waypoint_eigs=n_waypoint_eigs,
                    convergence_epsilon = 1e-5)
    model.construct_kernel_matrix()
    model.initialize_archetypes()
    try:
        model.fit(min_iter=min_iter, max_iter=max_iter)
    except RuntimeWarning as w:
        warnings.warn(str(w), SEACellsWarning)
    print(f'Ran for {len(model.RSS_iters)} iterations')
    return model


def get_assignment(model, hard=True):
    archetype_labels = model.get_hard_archetypes()
    if hard:
        df=model.get_hard_assignments()        
        arch_bc=[archetype_labels[i] for i in model.A_.argmax(0)]
        df.insert(0,"archetype_barcode",arch_bc)
    else:
        labels, weights=model.get_soft_assignments()
        table = []
        metacell_dict={bc:f'SEACell-{i}' for i,bc in enumerate(archetype_labels)}
        for cell_i,(cell_bc, row) in enumerate(labels.iterrows()):
            row_data = [cell_bc]
            for i in range(5):
                metacell_bc = row[i]
                metacell_id = metacell_dict[metacell_bc]
                weight = weights[cell_i][i]
                row_data.extend([metacell_id, metacell_bc, weight])
            table.append(row_data)
        # Convert the table to a DataFrame
        df = pd.DataFrame(table, columns=['cell_barcode'] + 
                                 [f'metacell{i+1}_{v}' for i in range(5)  for v in ('id', 'bc', 'weight') ] ).set_index('cell_barcode')
    return df


def export_seacells(seacells_obj,gene_dict, h5_out=None):
    "export metacells agregated counts to h5 file"
    if h5_out is not None:
        with h5py.File(h5_out, 'w') as f:
            # Save counts
            f.create_dataset('raw_counts', data=seacells_obj.layers['raw'].toarray())
            # Save cell names (metacell names)
            metacell_names = np.array(seacells_obj.obs_names, dtype='S')
            f.create_dataset('cell_attrs/cell_names', data=metacell_names)
            # Save gene names
            gene_ids = np.array(seacells_obj.var_names, dtype='S')
            f.create_dataset('gene_attrs/gene_ids', data=gene_ids)
            # Save gene names
            gene_names = np.array([gene_dict[gid] for gid in seacells_obj.var_names], dtype='S')
            f.create_dataset('gene_attrs/gene_names', data=gene_names)
            # save the fractionMT
            fraction_mt=seacells_obj.obs['fractionMT'].values
            f.create_dataset('cell_attrs/fractionMT', data=fraction_mt)
            # Indicate that cells are on rows
            f.create_dataset('cell_attrs/cells_on_rows', data=np.array(True))
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Running SEACells metacells identification")
    parser.add_argument("-i", metavar='in_file.h5', required=True, help="Input file with raw counts")
    parser.add_argument("-o", metavar='out_dir', required=True, help="Output directory")
    parser.add_argument("-p", metavar="out_prefix", required=True, help="Prefix for the output files")
    parser.add_argument("-c", metavar="celltype_mapping", help="Optional: File with barcode celltype mappings. Used for celltype purity of metacells.",default=None)
    parser.add_argument("--n_cells", type=int, default=90, help="Number of metacells [default: 90]")
    parser.add_argument("--n_pc", type=int, default=50, help="PCA dimensions [default: 50]")
    parser.add_argument("--n_waypoint_eigs", type=int, default=10, help="Number of eigenvalues to consider when initializing metacells [default: 10]")
    parser.add_argument("--n_topvar_genes", type=int, default=1500, help="Number of top variable genes [default: 1500]")
    parser.add_argument("--n_neighbors", type=int, default=15, help="Number of neighbors considered for the neighborhood graph [default: 15]")
    parser.add_argument("--iter", type=str, default="10,50", help="Number of iterations to fit, provided as min,max [default: 10,50]")

    args = parser.parse_args()
    # parse the tags
    h5file=args.i
    out_dir=args.o
    out_prefix=args.p
    celltype_mapping_file=args.c
    n_SEACells = int(args.n_cells)
    n_pc = int(args.n_pc)
    n_waypoint_eigs = int(args.n_waypoint_eigs)
    n_topvar_genes = int(args.n_topvar_genes)
    n_neighbors = int(args.n_neighbors)
    min_iter, max_iter = tuple(map(int, args.iter.split(',')))
    print(f"input file: {h5file}")
    print(f"output files: {out_dir}/{out_prefix}*\n")
    if celltype_mapping_file is not None:
        print(f"celltype mapping: {celltype_mapping_file}")
    else:
        print("no celltype mapping provided")
    main(h5file, out_dir, out_prefix, celltype_mapping_file, n_SEACells, n_pc, n_waypoint_eigs, n_topvar_genes, n_neighbors, min_iter, max_iter)
    



