import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="parser")

# Mandatory argument
parser.add_argument("--source_file", type=str, default="/cbica/home/wenju/Dataset/UKBB/UKBB_genetic_preprocess/S4_GWAS/All_ancestry",
                           help="Path to store the clustering results")
parser.add_argument("--target_file", type=str, default="/cbica/home/wenju/Dataset/UKBB/UKBB_genetic_preprocess/S4_GWAS/All_ancestry",
                           help="Path to store the clustering results")
def convert(source_file, target_file):
    df = pd.read_csv(source_file, sep='\t')

    df.rename({'AF1': 'freq'}, axis=1, inplace=True)
    df.rename({'BETA': 'b'}, axis=1, inplace=True)
    df.rename({'SE': 'se'}, axis=1, inplace=True)
    df.rename({'P': 'p'}, axis=1, inplace=True)
    df['p'] = df['p'].astype(float)
    df['b'] = df['b'].astype(float)
    df['freq'] = df['freq'].astype(float)
    df['se'] = df['se'].astype(float)
    df['N'] = df['N'].astype(int)
    df = df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']]
    df.to_csv(target_file, index=False, sep=' ')
    print('Here...')
def main(options):
    convert(options.source_file, options.target_file)
if __name__ == "__main__":
    commandline = parser.parse_known_args()
    options = commandline[0]
    if commandline[1]:
        raise Exception("unknown arguments: %s" % parser.parse_known_args()[1])
    main(options)
