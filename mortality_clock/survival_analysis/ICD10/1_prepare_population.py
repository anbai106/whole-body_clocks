import os.path
import pandas as pd
import numpy as np

# ---- Paths ----
### 7 mri mortality clock
for organ in ['brain', 'adipose', 'heart', 'kidney', 'liver', 'pancreas', 'spleen']:
    organ_big = organ.capitalize()
    tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock', organ + "_mri_mortality_clock", organ + "_mri_mortality_clock_predictions.tsv")
    df_data = pd.read_csv(tsv, sep='\t')
    df_data.rename({organ + '_mri_mortality_clock_acceleration_z': organ_big + '_mri'}, axis=1, inplace=True)
    df_data = df_data[['participant_id', organ_big + '_mri']]
    if organ == 'brain':
        df_clock = df_data
    else:
        df_clock = pd.merge(df_clock, df_data, on="participant_id", how="outer")

### 11 proteomics mortality clock
for organ in ['Reproductive_female', 'Pulmonary', 'Heart', 'Brain', 'Eye', 'Hepatic', 'Renal', 'Reproductive_male', 'Endocrine', 'Immune', 'Skin']:
    organ_small = organ.lower()
    tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock', organ + "_proteomics_mortality_clock", organ_small + "_proteomics_mortality_clock_predictions.tsv")
    df_data = pd.read_csv(tsv, sep='\t')
    df_data.rename({organ_small + '_proteomics_mortality_clock_acceleration_z': organ + '_proteomics'}, axis=1, inplace=True)
    df_data = df_data[['participant_id', organ + '_proteomics']]
    df_clock = pd.merge(df_clock, df_data, on="participant_id", how="outer")

### 5 metabolomics mortality clock
for organ in ['Endocrine', 'Digestive', 'Hepatic', 'Immune']:
    organ_small = organ.lower()
    tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock', organ + "_metabolomics_mortality_clock", organ_small + "_metabolomics_mortality_clock_predictions.tsv")
    df_data = pd.read_csv(tsv, sep='\t')
    df_data.rename({organ_small + '_metabolomics_mortality_clock_acceleration_z': organ + '_metabolomics'}, axis=1, inplace=True)
    df_data = df_data[['participant_id', organ + '_metabolomics']]
    df_clock = pd.merge(df_clock, df_data, on="participant_id", how="outer")
bag_population_list = df_clock['participant_id']

## ICD-10 code
icd10_tsv = '/cbica/home/wenju/Reproducibile_paper/BrainEye/data/UKBB_fullsample_ICD10.csv'
df_icd_diagnosis = pd.read_csv(icd10_tsv, sep=",")
df_icd_diagnosis.rename({'eid': 'participant_id'}, axis=1, inplace=True)
### rename column
col_names_diagnosis = df_icd_diagnosis.columns.to_list()
col_names_diagnosis_new = []
for e in col_names_diagnosis:
    if e.startswith('diagnoses_icd10_f'):
        e_new = e.replace("diagnoses_icd10_f41270_", "")
    else:
        e_new = e
    col_names_diagnosis_new.append(e_new)
df_icd_diagnosis.columns = col_names_diagnosis_new

### ICD-10 code date of first in-patient diagnosis
df_ukb_death = pd.read_excel('/cbica/home/wenju/Dataset/UKBB_UMelbourne/Death_related_var_from_Ye.xlsx')
df_id_match = pd.read_csv('/cbica/home/wenju/Dataset/UKBB_UMelbourne/UKB_UMelbourne_vs_Penn_match_key.csv')
df_ukb_death.rename({'eid': 'participant_id_umel'}, axis=1, inplace=True)
df_id_match.rename({'id': 'participant_id_umel'}, axis=1, inplace=True)
df_id_match.rename({'id_upenn': 'participant_id'}, axis=1, inplace=True)
df_icd = df_id_match.merge(df_ukb_death, on='participant_id_umel')
df_icd_participant_id = df_icd[['participant_id']]
df_icd_diagnosis_date = df_icd[df_icd.columns[df_icd.columns.str.startswith('41280-')]]
df_icd_diagnosis_date = pd.concat([df_icd_participant_id, df_icd_diagnosis_date], axis=1)
### rename column
col_names_date = df_icd_diagnosis_date.columns.to_list()
col_names_date_new = []
for e in col_names_date:
    if e.startswith('41280-'):
        e_new = e.replace("41280-", "")
        e_new = e_new.replace(".", "_")
    else:
        e_new = e
    col_names_date_new.append(e_new)
df_icd_diagnosis_date.columns = col_names_date_new
df_icd_diagnosis_date.dropna(axis=1, how='all', inplace=True)

### merge the overlap population between the diagnosis and date
date_pop = df_icd_diagnosis_date['participant_id']. to_list()
diagnosis_pop = df_icd_diagnosis['participant_id']. to_list()
pop_overlap = [value for value in date_pop if value in diagnosis_pop]
df_icd_diagnosis_date = df_icd_diagnosis_date.loc[df_icd_diagnosis_date['participant_id'].isin(pop_overlap)]
df_icd_diagnosis = df_icd_diagnosis.loc[df_icd_diagnosis['participant_id'].isin(pop_overlap)]

########################################################################################################################
### bag
########################################################################################################################
print("########################################################################################################################")
### merge with the brain population
df_icd_brain = df_icd_diagnosis.loc[df_icd_diagnosis['participant_id'].isin(bag_population_list)]
df_icd_data = df_icd_brain.iloc[:, 1:]
cn_boolean = df_icd_data.isnull().all(1)
df_icd_cn = df_icd_brain[cn_boolean.values]
df_icd_pt = df_icd_brain[~cn_boolean.values]
df_icd_pt_pop = df_icd_pt[['participant_id']]
df_icd_pt_data = df_icd_pt.iloc[:, 1:]
df_icd_cn.replace(np.nan, -1, inplace=True)
df_icd_cn = df_icd_cn[['participant_id', '0_0']]

### unique PT diagnosis
# Melt the dataframe to have all values in a single column
melted_df = df_icd_pt_data.melt()
# Get the value counts
counts_pt = melted_df['value'].value_counts().reset_index()
counts_pt.columns = ['value', 'count']
counts_pt = counts_pt.sort_values(['count'], ascending=False)
counts_pt_larger_50 = counts_pt.loc[(counts_pt['count'] >= 50)]
counts_pt_larger_50.to_csv('/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/included_ICD_mortality_clock.tsv', index=False,
    sep='\t', encoding='utf-8')
pt_list = counts_pt_larger_50['value'].to_list()
print('Here...')