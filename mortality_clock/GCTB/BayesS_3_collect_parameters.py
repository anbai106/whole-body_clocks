import numpy as np
import pandas as pd
import argparse
import os
from scipy.stats import ttest_ind, ttest_ind_from_stats
def collect_gctb(output_dir, output_dir_result):
    """
    Plot the Fig. 3B for Ioanna's paper
    """

    df_final = pd.DataFrame(columns=['BAG', 'h2_mean', 'h2_se', 'S', 'S_se', 'Pi', 'Pi_se'])
    ### bag
    bag_list = ['brain', "adipose", "heart", "kidney", "liver", "pancreas", "spleen"]
    for i in range(len(bag_list)):
        pheno = bag_list[i]
        file = os.path.join(output_dir, pheno, pheno + '.parRes')
        if os.path.exists(file):
            print("bag is: %s" % pheno)
            df = pd.read_csv(file, delim_whitespace=True)
            h2_mean = df.iloc[[7]]['Mean'].values[0]
            h2_se = df.iloc[[7]]['SD'].values[0]
            S = df.iloc[[3]]['Mean'].values[0]
            S_se = df.iloc[[3]]['SD'].values[0]
            Pi = df.iloc[[0]]['Mean'].values[0]
            Pi_se = df.iloc[[0]]['SD'].values[0]
            dict = {'BAG': pheno, 'h2_mean': h2_mean, 'h2_se': h2_se, 'S': S, 'S_se': S_se, 'Pi': Pi,
                    'Pi_se': Pi_se}
            df_final = df_final.append(dict, ignore_index=True)
        else:
            print("bag cannot converge for SBayesS: %s" % pheno)
    df_final.to_csv(os.path.join(output_dir_result, 'GCTB_SBayesS_parameters.tsv'), index=False, sep='\t', encoding='utf-8')

    print("STOP here...")

output_dir = '/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/GCTB'
output_dir_result = "/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/Result"
collect_gctb(output_dir, output_dir_result)
