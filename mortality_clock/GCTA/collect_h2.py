import pandas as pd
import glob, os
import numpy as np
import rpy2.robjects as ro


def barplot_h2(output_dir, r_script_path, output_dir_result):
    """
    Plot the Fig. 3B for Ioanna's paper
    """

    mribag_list = ['brain', "adipose", "heart", "kidney", "liver", "pancreas", "spleen"]

    ### genotype array
    df_final = pd.DataFrame(columns=['mribag', 'Heritability', 'SE', 'Pvalue'])
    for mribag in mribag_list:
        file = os.path.join(output_dir, mribag, mribag + '.hsq')
        df = pd.read_csv(file, delimiter='\t', header=None)
        coefficient = np.float64(df.iloc[4][1])
        coefficient_std = np.float64(df.iloc[4][2])
        pvalue = np.float64(df.iloc[9][1])
        lrt = np.float64(df.iloc[7][1])
        if df.iloc[9][1] == '0.0000e+00':
            r = ro.r
            r.source(r_script_path)
            from rpy2.robjects import globalenv
            globalenv['lrt'] = float(lrt)
            pvalue = r.compute_p(float(lrt))[0]
        dict = {'mribag': mribag, 'Heritability': coefficient, 'SE': coefficient_std, 'Pvalue': pvalue}
        df_final = df_final.append(dict, ignore_index=True)

    df_final.to_csv(os.path.join(output_dir_result, 'GCTA_h2_results_genotype_array.tsv'), index=False, sep='\t', encoding='utf-8')

    print("STOP here...")

    ### use r to plot the barplot

output_dir = "/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/h2_gcta"
output_dir_result = "/Users/hao/cubic-home/Reproducibile_paper/AbdoImaging/Result"
r_script_path = '/Users/hao/cubic-home/Project/AbdoImaging/heritability_gcta/convert_p_in_R.R'
barplot_h2(output_dir, r_script_path, output_dir_result)