import pandas as pd
import os
from qmplot import manhattanplot, qqplot
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="parser")

# Mandatory argument
parser.add_argument("--output_dir", type=str, default="/cbica/home/wenju/Dataset/UKBB/UKBB_genetic_preprocess/S4_GWAS/All_ancestry",
                           help="Path to store the clustering results")
parser.add_argument("--output_result", type=str, default="/cbica/home/wenju/Dataset/UKBB/UKBB_genetic_preprocess/S4_GWAS/All_ancestry",
                           help="Path to store the clustering results")

def plot_qmplot(output_dir, output_result):
    df = pd.read_csv(output_result, sep='\t')
    ## plot the manhattan plot
    df.rename({'CHR': '#CHROM'}, axis=1, inplace=True)
    df.dropna(inplace=True)
    df.replace([np.inf, -np.inf], 1e-200, inplace=True)
    df['P'][df['P'] < 1e-200] = 1e-200
    # clean data
    manhattanplot(data=df,
                       xticklabel_kws={"rotation": "vertical"}, sign_marker_p=None, genomewideline=None,
                       suggestiveline=5e-8)
    plt.savefig(os.path.join(output_dir, 'manhattan_qmplot.png'))
    qqplot(data=df["P"], marker="o",
           xlabel=r"Expected $-log_{10}{(P)}$",
           ylabel=r"Observed $-log_{10}{(P)}$")
    plt.savefig(os.path.join(output_dir, 'QQ_plot.png'))
    print("STOP ...")

# output_dir="/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/pQTL/CSF2RB/fastGWA"
# output_result="/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/pQTL/CSF2RB/fastGWA/protein_pheno_normalized_residualized.fastGWA"
# plot_qmplot(output_dir, output_result)
def main(options):
    plot_qmplot(options.output_dir, options.output_result)
if __name__ == "__main__":
    commandline = parser.parse_known_args()
    options = commandline[0]
    if commandline[1]:
        raise Exception("unknown arguments: %s" % parser.parse_known_args()[1])
    main(options)


