import pandas as pd
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
matplotlib.use('Agg')  # Use Agg backend for rendering plots
import argparse
import warnings



# test call
# sample="MAHACEB-T"
# celltypes=f"results/celltyping/{sample}.cts_final.txt"
# metacelltypes=f"results/celltyping/{sample}_seacells.cts_final.txt"
# metacell_assignment=f"results/metacells/{sample}.genes_cells_filtered_seacells_assignment.tsv"
# main(celltypes, metacelltypes, metacell_assignment, "metacell_tests/SEACells/cmp", sample)

def main(celltype_file, metacelltype_file, assignment_file, out_dir, out_prefix, reportfile):
    metacell_counts=get_metacell_counts(celltype_file, metacelltype_file, assignment_file)
    metacell_counts.to_csv(f"{out_dir}/{out_prefix}_metacell_celltype_counts.tsv", sep='\t', index=False)
    with open(reportfile, 'w') as f:
        f.write(get_report(metacell_counts))
    # Plotting
    plt.figure(figsize=(10, 6))
    sns.histplot(metacell_counts['purity'], kde=True, bins=30)
    # Adding titles and labels
    plt.title('Distribution of Metacell Purity')
    plt.xlabel('Purity')
    plt.ylabel('Density')
    plt.savefig(f'{out_dir}/{out_prefix}_seacells_celltype_hist.png')

    sns.kdeplot(metacell_counts['purity'], shade=True)
    # Adding titles and labels
    plt.title('Distribution of Metacell Purity')
    plt.xlabel('Purity')
    plt.ylabel('Density')
    plt.savefig(f'{out_dir}/{out_prefix}_seacells_celltype_dens.png')
  
def get_report(metacell_counts):
    # print some statistics for the logfile
    n=len(metacell_counts) # number of metacells
    # check metacells with >20% unknown/uncertain
    un_cells=metacell_counts.query("composition.str.contains('un')")
    report = f'{len(un_cells)}/{n} metacells contain > 20% cells with unknown/uncertain celltypes:\n\n\n'
    report += un_cells.to_markdown(tablefmt='grid')
    # check very mixed cells
    mixed=metacell_counts.query("purity<.5 and dominant == 'mixed'")
    report += f'\n\n\n{len(mixed)}/{n} metacells have mixed celltypes with purity < 50%:\n\n\n'
    report += mixed.to_markdown(tablefmt='grid')
    # check for disagreements
    n_pure_all=len(metacell_counts.query("purity>.8"))
    pure_disagree=metacell_counts.query("purity>.8 and metacell_type != dominant")
    groups=["B.cell", "T.cell", "Macrophage", "Melanoma"]
    disagree_idx=[idx for idx,row in pure_disagree.iterrows() if not any(gr in row.metacell_type and gr in row.dominant for gr in groups)]
    report += f'\n\n\n{len(pure_disagree)}/{n_pure_all} cells with pure celltype disagree with aggregated celltype call, '
    report += f'of which {len(disagree_idx)} are from different groups:\n\n\n'
    report += pure_disagree.loc[disagree_idx].to_markdown(tablefmt='grid')
    return report


def get_metacell_counts(celltype_file, metacelltype_file, assignment_file):
    """ 
    Get a metacell count table

    From the celltype calls of individual cells and metacells, calculate a table with the following infos per metacell:
    * number of cells
    * called celltype for aggregated cells
    * dominant celltype for individual cells
    * celltype composition in individual cells
    * purity (e.g. fraction of predominant celltype)
    * number of cells for individual celltypes
    """
    # Read the data
    celltypes_df = pd.read_csv(celltype_file, sep='\t')
    metacelltypes_df = pd.read_csv(metacelltype_file, sep='\t')
    metacelltypes_dict=dict(metacelltypes_df.set_index('barcodes')['celltype_final'])
    metacell_assignment_df = pd.read_csv(assignment_file, sep='\t')
    # Merge the cell type info with the metacell assignment info
    merged_df = pd.merge(metacell_assignment_df[['cell_barcode', 'metacell1_id']], celltypes_df, left_on='cell_barcode', right_on='barcodes', how='left')
    merged_df = merged_df[['cell_barcode', 'metacell1_id', 'celltype_final']]
    # Rename columns for clarity
    merged_df.columns = ['cell_barcode', 'metacell_id', 'celltype']
    # Calculate the cell counts for each metacell
    metacell_counts = merged_df.groupby('metacell_id')['celltype'].value_counts().unstack().fillna(0).astype(int)
    # Calculate the purity for each metacell
    metacell_purity = merged_df.groupby('metacell_id')['celltype'].value_counts(normalize=True).unstack().fillna(0)
    #dominant = metacell_purity.idxmax(axis=1)
    dominant, composition = [], []
    # define "dominant" celltype and composition per metacell
    # composition contains all celltypes assigned to >20% of the cells
    # if a metacell is composed, according to definition above, of more than one celltype
    # the dominant celltype is set to "mixed". 
    # If all celltypes of a mixed metacell contain one of the group substrings, 
    # they are called group_mixed (e.g. "T.cell_mixed")
    groups=["B.cell", "T.cell", "Macrophage", "Melanoma"]
    for cid, row in metacell_purity.iterrows():
        # sort by contribution
        fractions=row.sort_values(ascending=False)
        # consider all celltypes assigned to > 20% of the cells
        relevant=[celltype for celltype, frac in fractions.items() if frac>.2]
        # format to string
        composition_string = ", ".join([f"{ct} ({row[ct]:.1%})" for ct in relevant])
        # in case there is just one celltype
        dominant_string=relevant[0]
        # else its mixed
        if len(relevant)>1:
            dominant_string='mixed'
            # check for each group
            for gr in groups:
                # whether all celltypes contain the groupname as substring
                if all(gr in ct for ct in relevant):
                    # in this case it gets group_mixed
                    dominant_string=f'{gr}_mixed'                
        dominant.append(dominant_string)
        composition.append(composition_string)

    # add the gathered information to the cell count table
    n_cells=metacell_counts.sum(1)
    metacell_type=[metacelltypes_dict.get(cid, 'no_call') for cid in metacell_counts.index]
    metacell_counts.insert(0,'purity', metacell_purity.max(axis=1))
    metacell_counts.insert(0,'dominant', dominant)
    metacell_counts.insert(1,'composition', composition)
    metacell_counts.insert(0,'metacell_type', metacell_type)
    metacell_counts.insert(0,'n_cells', n_cells)

    return metacell_counts
    


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Running SEACells metacells identification")
    parser.add_argument("-i", metavar='individual_celltypes.cts_final.txt', required=True, help="celltypes called for individual cells from celltyping.R (rule celltyping)")
    parser.add_argument("-m", metavar='metacelltypes.cts_final.txt', required=True, help="celltypes called for aggregated metacells from celltyping.R (rule celltyping)")
    parser.add_argument("-a", metavar="metacell_assignment.tsv", required=True, help="table with cell to metacell assignment from run_seacells.py")
    parser.add_argument("-o", metavar="output_dir", required=True, help="output_directory")
    parser.add_argument("-p", metavar="prefix", required=True, help="prefix used for output files")
    parser.add_argument("-r", metavar="path/to/reportfile.rst",help="report is written to this file", required=False, default=None)

    args = parser.parse_args()
    # parse the tags
    celltype_file=args.i
    metacelltype_file=args.m
    assignment_file=args.a
    out_dir=args.o
    out_prefix=args.p
    reportfile=args.r
    main(celltype_file, metacelltype_file, assignment_file, out_dir, out_prefix, reportfile)

