import pandas as pd
import os, random
import numpy as np
from sklearn.model_selection import ShuffleSplit
from scipy.stats import ttest_ind, chi2_contingency

def chi_square_contigency_table_independence(df_1, df_2):
    """Chi square test for categorical variables for independence
    :param df_1: dataframe to contain the subytpes info
    :param df_2: dataframe to contain the categorical info

    :return: chi-square test statistic and p-value
    ref: https://en.wikipedia.org/wiki/Chi-squared_test
    """
    contingency = pd.crosstab(df_1, df_2)

    # Chi-square test of independence.
    c, p, dof, expected = chi2_contingency(contingency.values)

    return c, p

def prepare_data(output_dir):
    """
    Split-sample analysis
    Returns:
    """
    for organ in ["Metabolic", "Immune"]:
        print("Split-sample analysis for %s" % organ)
        df_organ = pd.read_csv(os.path.join(output_dir, organ, 'BAG_pheno.txt'), sep='\t')
        df_cov = pd.read_csv(os.path.join(output_dir, organ, 'BAG_cov.txt'), sep='\t')
        df_keep = pd.read_csv(os.path.join(output_dir, organ, 'BAG_keep_for_fastgwa.txt'), sep='\t', header=None)

        if df_organ.shape[0] != df_keep.shape[0] or df_organ.shape[0] != df_cov.shape[0]:
            raise Exception("Someting is wrong here...")
        else:
            print("The FID and IID for the 3 dataframes should be identical")

        #### split sample by age- and sex-matched two splits
        output_dir_female = os.path.join(output_dir, organ, 'female')
        if not os.path.exists(output_dir_female):
            os.makedirs(output_dir_female)
        output_dir_male = os.path.join(output_dir, organ, 'male')
        if not os.path.exists(output_dir_male):
            os.makedirs(output_dir_male)

        if os.path.exists(os.path.join(output_dir_female, 'BAG_pheno.txt')) and os.path.exists(os.path.join(output_dir_male, 'BAG_pheno.txt')):
            continue
        else:
            list_female = list(df_cov.loc[df_cov['sex_f31_0_0'].isin([0])]['FID'])
            list_male = list(df_cov.loc[df_cov['sex_f31_0_0'].isin([1])]['FID'])

            df_organ_female = df_organ.loc[df_organ['FID'].isin(list_female)]
            df_organ_female.to_csv(os.path.join(output_dir_female, 'BAG_pheno.txt'), index=False, sep='\t')
            df_organ_male = df_organ.loc[df_organ['FID'].isin(list_male)]
            df_organ_male.to_csv(os.path.join(output_dir_male, 'BAG_pheno.txt'), index=False, sep='\t')

            df_cov_female = df_cov.loc[df_cov['FID'].isin(list_female)]
            df_cov_female.to_csv(os.path.join(output_dir_female, 'BAG_cov.txt'), index=False, sep='\t')
            df_cov_male = df_cov.loc[df_cov['FID'].isin(list_male)]
            df_cov_male.to_csv(os.path.join(output_dir_male, 'BAG_cov.txt'), index=False, sep='\t')

            df_keep_female = df_keep.loc[df_keep[0].isin(list_female)]
            df_keep_female.to_csv(os.path.join(output_dir_female, 'BAG_keep_for_fastgwa.txt'), index=False,
                                  sep='\t',
                                  header=False)
            df_keep_male = df_keep.loc[df_keep[0].isin(list_male)]
            df_keep_male.to_csv(os.path.join(output_dir_male, 'BAG_keep_for_fastgwa.txt'), index=False,
                                sep='\t',
                                header=False)
            
    print("Stop here... ")

output_dir = '/Users/hao/cubic-home/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/data/'
prepare_data(output_dir)