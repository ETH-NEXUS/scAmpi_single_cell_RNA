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


def main(input_files, output_file):
    plot_purity(input_files, output_file)
    # potentially add more plots and tables


def plot_purity(ct_mapping_files, out_file):
    purity = {}
    for fn in ct_mapping_files:
        sample = os.path.basename(fn).split('_')[0]
        ct_mapping = pd.read_csv(fn, sep='\t')
        purity[sample] = ct_mapping.purity
    purity = pd.DataFrame(purity)
    stats = pd.DataFrame(dict(median=purity.median(), mean=purity.mean()))
    sample_order = stats.sort_values(
        by=['median', 'mean'], ascending=False).index
    purity = purity.melt(var_name='Sample', value_name='purity')
    plt.figure(figsize=(14, 8))
    # sns.violinplot(x='Sample', y='purity', data=purity, order=sample_order)
    sns.boxplot(x='Sample', y='purity', data=purity, order=sample_order)
    plt.xticks(rotation=90)  # Rotate sample labels for better readability
    plt.title('celltype purity of metacells across samples')
    plt.tight_layout()
    plt.savefig(out_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create boxplots of metacell \
            celltype purity across samples")
    parser.add_argument("output_file", metavar='out.png',
                        help="Paths to the output png file.")
    parser.add_argument('input_files', type=str, nargs='+',
                        help="Paths to the input CSV files\
                              containing the purity information.")
    args = parser.parse_args()
    main(args.input_files, args.output_file)
