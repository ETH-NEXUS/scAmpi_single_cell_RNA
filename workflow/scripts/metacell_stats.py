# Author: Matthias Lienhard
# Date: 2024-08-22
# Description: Aggregate statistics on metacells.
#              So far, only a summary plot on the
#              celltype purity is implemented.


import matplotlib.pyplot as plt
import os
import argparse
import seaborn as sns
import pandas as pd
import matplotlib
matplotlib.use('Agg')


def main(input_files, output_prefix):
    plot_purity(input_files, output_prefix, True)
    plot_purity(input_files, output_prefix, False)
    plot_cells_per_metacell(input_files, output_prefix)
    # potentially add more plots and tables


def plot_purity(ct_mapping_files, out_prefix, subcelltypes=False):
    if subcelltypes:
        out_filename = f'{out_prefix}subcelltype_purity.png'
        col_name = "sub_purity"
        title = 'Celltype Purity of Metacells Across Samples'
    else:
        out_filename = f'{out_prefix}celltype_purity.png'
        col_name = "purity"
        title = 'Subcelltype Purity of Metacells Across Samples'
    purity, sample_order = get_stats(ct_mapping_files, col_name)
    plt.figure(figsize=(14, 8))
    # sns.violinplot(x='Sample', y='purity', data=purity, order=sample_order)
    sns.boxplot(x='Sample', y=col_name, data=purity, order=sample_order)
    plt.xticks(rotation=90)  # Rotate sample labels for better readability
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_filename)

def get_stats(ct_mapping_files, col_name):
    vals = {}
    for fn in ct_mapping_files:
        sample = os.path.basename(fn).split('_')[0]
        ct_mapping = pd.read_csv(fn, sep='\t')
        vals[sample] = ct_mapping[col_name]
    df = pd.DataFrame(vals)
    stats = pd.DataFrame(dict(median=df.median(), mean=df.mean()))
    sample_order = stats.sort_values(
        by=['median', 'mean'], ascending=False).index
    df = df.melt(var_name='Sample', value_name=col_name)
    return df, sample_order


def plot_cells_per_metacell(ct_mapping_files, out_prefix):
    n_cells, sample_order = get_stats(ct_mapping_files, "n_cells")
    plt.figure(figsize=(14, 8))
    # sns.violinplot(x='Sample', y='purity', data=purity, order=sample_order)
    sns.boxplot(x='Sample', y='n_cells', data=n_cells, order=sample_order)
    plt.xticks(rotation=90)  # Rotate sample labels for better readability
    plt.title("Number of Cells per Metacell")
    plt.tight_layout()
    logger.info(f'writing image to ')
    plt.savefig(f'{out_prefix}cells_per_metacell.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create boxplots of metacell \
            celltype purity across samples")
    parser.add_argument("output_prefix", metavar='path/to/output/dir/prefix_',
                        help="Path to the output directory.")
    parser.add_argument('input_files', type=str, nargs='+',
                        help="Paths to the input CSV files\
                              containing the purity information.")
    args = parser.parse_args()
    main(args.input_files, args.output_prefix)
