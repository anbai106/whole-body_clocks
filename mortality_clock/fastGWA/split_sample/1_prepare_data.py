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
    for organ in ['Endocrine', 'Digestive', 'Hepatic', 'Immune', "Metabolic"]:
        print("Split-sample analysis for %s" % organ)
        df_organ = pd.read_csv(os.path.join(output_dir, organ, 'BAG_pheno.txt'), sep='\t')
        df_cov = pd.read_csv(os.path.join(output_dir, organ, 'BAG_cov.txt'), sep='\t')
        df_keep = pd.read_csv(os.path.join(output_dir, organ, 'BAG_keep_for_fastgwa.txt'), sep='\t', header=None)

        if df_organ.shape[0] != df_keep.shape[0] or df_organ.shape[0] != df_cov.shape[0]:
            raise Exception("Someting is wrong here...")
        else:
            print("The FID and IID for the 3 dataframes should be identical")

        #### split sample by age- and sex-matched two splits
        output_dir_split1 = os.path.join(output_dir, organ, 'split1')
        if not os.path.exists(output_dir_split1):
            os.makedirs(output_dir_split1)
        output_dir_split2 = os.path.join(output_dir, organ, 'split2')
        if not os.path.exists(output_dir_split2):
            os.makedirs(output_dir_split2)
        flag_selection = True
        n_try = 0
        y = np.array(df_cov['age_when_attended_assessment_centre_f21003_0_0'])
        sex = df_cov['sex_f31_0_0']
        age = df_cov['age_when_attended_assessment_centre_f21003_0_0']
        tval_chi2_threshold = 0.2
        pval_threshold_ttest = 0.8

        while flag_selection:
            splits = ShuffleSplit(n_splits=1, test_size=0.5)
            splits_indices = list(splits.split(np.zeros(len(y)), y))[0]
            split_1_indice = splits_indices[0].tolist()
            split_2_indice = splits_indices[1].tolist()
            if len(split_1_indice) > len(split_2_indice):
                ### remove on index of the lagger split
                random_element = random.choice(split_1_indice)
                split_1_indice.remove(random_element)
            elif len(split_1_indice) < len(split_2_indice):
                random_element = random.choice(split_2_indice)
                split_2_indice.remove(random_element)
            else:
                pass
            ## check if the age and sex distribution differ
            t_age, p_age = ttest_ind(age[split_1_indice].to_numpy(), age[split_2_indice].to_numpy())
            T_sex, _ = chi_square_contigency_table_independence(sex[split_1_indice].to_numpy(), sex[split_2_indice].to_numpy())

            if T_sex < tval_chi2_threshold and p_age > pval_threshold_ttest:
                flag_selection = False

                df_organ_split_1 = df_organ.iloc[split_1_indice]
                df_organ_split_1.to_csv(os.path.join(output_dir_split1, 'BAG_pheno.txt'), index=False, sep='\t')
                df_organ_split_2 = df_organ.iloc[split_2_indice]
                df_organ_split_2.to_csv(os.path.join(output_dir_split2, 'BAG_pheno.txt'), index=False, sep='\t')

                df_cov_split_1 = df_cov.iloc[split_1_indice]
                df_cov_split_1.to_csv(os.path.join(output_dir_split1, 'BAG_cov.txt'), index=False, sep='\t')
                df_cov_split_2 = df_cov.iloc[split_2_indice]
                df_cov_split_2.to_csv(os.path.join(output_dir_split2, 'BAG_cov.txt'), index=False, sep='\t')

                df_keep_split_1 = df_keep.iloc[split_1_indice]
                df_keep_split_1.to_csv(os.path.join(output_dir_split1, 'BAG_keep_for_fastgwa.txt'), index=False, sep='\t',
                                header=False)
                df_keep_split_2 = df_keep.iloc[split_2_indice]
                df_keep_split_2.to_csv(os.path.join(output_dir_split2, 'BAG_keep_for_fastgwa.txt'), index=False, sep='\t',
                                header=False)
            print("Number of tries %d" % n_try)
            n_try += 1
    print("Stop here... ")

output_dir = '/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/data/'
prepare_data(output_dir)