import pandas as pd
import os

def collect_data(output_dir_result):
    df_final = pd.DataFrame(columns=['BAG', 'DE', 'hazard_ratio', 'CI_lower_bound', 'CI_upper_bound', 'p_value', 'N_case', 'N_noncase'])
    df_pheno_list = pd.read_csv('/Users/hao/cubic-home/Reproducibile_paper/FemaleAgingClock/survival_analysis/data/included_ICD_disease_femalebag.tsv', sep='\t')
    pheno_list = df_pheno_list['value'].to_list()
    for pheno in pheno_list:
        tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/FemaleAgingClock/survival_analysis/output/full_cov', 'cox_hr_' + pheno + '_female.tsv')
        if not os.path.exists(tsv):
            print("The case is too few after merging the date and diagnosis data from UMelbourne and UPenn...: %s" % pheno)
        else:
            df = pd.read_csv(tsv, sep='\t')
            df=df.rename(columns={'var':'BAG'})
            df['DE'] = pheno
            df['Sex'] = 'Female'
            df = df[['BAG', 'DE', 'hazard_ratio', 'CI_lower_bound', 'CI_upper_bound', 'p_value', 'N_case', 'N_noncase']]
            df_final = df_final.append(df, ignore_index=True)

    #### remove the n_case < 200 to ensure power
    df_final = df_final.dropna(subset=df_final.columns[2:5], how='all')
    df_final = df_final[df_final["BAG"] != "Female_BAG"]
    df_final_size = df_final.loc[(df_final['N_case'] >= 10)]
    n_de = df_final_size['DE'].unique().size ### n_de=886 DEs
    n_bag = df_final_size['BAG'].unique().size ### n_bag=20 BAGs
    df_final_sig = df_final_size.loc[(df_final_size['p_value'] <= 0.05/38)]
    df_final_sig['p_value_bag_de'] = 0.05/20/n_de
    df_final_sig.to_csv(os.path.join(output_dir_result, 'Cox_HR_ICD10_BAG_sig_female_full_cov.tsv'), index=False, sep='\t', encoding='utf-8') ###29 unique diseases
    df_final_size.to_csv(os.path.join(output_dir_result, 'Cox_HR_ICD10_BAG_female_full_cov.tsv'), index=False, sep='\t', encoding='utf-8') ###29 unique diseases
    df_final_sig.to_excel(os.path.join(output_dir_result, 'Cox_HR_ICD10_BAG_sig_female_full_cov.xlsx'), index=False)
    print("Here...")

output_dir_result = '/Users/hao/cubic-home/Reproducibile_paper/FemaleAgingClock/Result'
collect_data(output_dir_result)