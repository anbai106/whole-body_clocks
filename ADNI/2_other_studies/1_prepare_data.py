import pandas as pd
#### here is to check what data columns we should use to develop the AD L'EPOCH
istaging_pickle_file = '/Users/hao/cubic-projects/ISTAGING/Pipelines/ISTAGING_Data_Consolidation_2020/v2.0/istaging.pkl.gz'
data = pd.read_pickle(istaging_pickle_file)
df_external = data.loc[data['Study'].isin(['ADNI_DOD', 'AIBL', 'BLSA', 'OASIS', 'PreventAD'])]
# df_adni_test = df_adni.iloc[:1]
# df_adni_test.to_csv("~/test.tsv", sep="\t", index=False)
df_external.to_csv("/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/adni_lepoch/external_istaging_adni_dod_aibl_blsa_oasis_prevent_ad.tsv", sep="\t", index=False)
print(df_external.columns.to_list())

print('Stop...')