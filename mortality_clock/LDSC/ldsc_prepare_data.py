import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="parser")

# Mandatory argument
parser.add_argument("--source_file", type=str, default="/cbica/home/wenju/Dataset/UKBB/UKBB_genetic_preprocess/S4_GWAS/All_ancestry",
                           help="Path to store the clustering results")

def prepare_4_phegwas(source_file):
    """
    This is to plot the number of significant independent SNP-Phenotype assocations for MIST and AAL.
    Args:
        mica_output_dir:

    Returns:

    """
    df = pd.read_csv(source_file, sep='\t')
    df['P'] = df['P'].astype(float)
    df['BETA'] = df['BETA'].astype(float)
    df['N'] = df['N'].astype(float)
    df['SE'] = df['SE'].astype(float)
    df = df[['SNP', 'A1', 'A2', 'N', 'P', 'BETA']]
    df.to_csv(source_file.split('.zip')[0] + '.ldsc.tsv', index=False, sep=' ')

def main(options):
    prepare_4_phegwas(options.source_file)
if __name__ == "__main__":
    commandline = parser.parse_known_args()
    options = commandline[0]
    if commandline[1]:
        raise Exception("unknown arguments: %s" % parser.parse_known_args()[1])
    main(options)