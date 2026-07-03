import pandas as pd
import os

### read the metabolite population
df_metabolite_long = pd.read_csv(
    '/Users/hao/cubic-home/Reproducibile_paper/UKBB_Proteomics/imputation_2_0/UKBB_preoteomics_real_data_all_ancestry_with_related_indi_4_autocomplete_final_imputation_2_0.tsv',
    sep='\t')

df_metabolite_long = pd.read_csv('/Users/hao/cubic-home/Dataset/UKBB/UKBB_NMR_metabolomics/QC_ukbnmr/nmr_biomarker_data.csv')
df_metabolite_long.rename({'eid': 'participant_id'}, axis=1, inplace=True)
df_metabolite_long = df_metabolite_long.loc[df_metabolite_long['visit_index'].isin([1])]

df_metabolite_long = df_metabolite_long.rename(columns={"ID": "participant_id"})
metabolite_list_long = df_metabolite_long.columns[2:].to_list()

# Convert the master list to a set once outside the loop for O(1) optimal speed
metabolite_set_long = set(metabolite_list_long)

output_dir = '/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/metabolomics/data'
organ_list = 'Reproductive_female', 'Pulmonary', 'Heart', 'Brain', 'Eye', 'Hepatic', 'Renal', 'Reproductive_male', 'Endocrine', 'Immune', 'Skin'

for organ in organ_list:
    df_metabolite_bl = pd.read_csv(
        '/Users/hao/cubic-home/Reproducibile_paper/UKBB_Proteomics/MLNI/data/' + organ + '/training/training_4589.tsv',
        sep='\t')
    metabolite_list_bl = df_metabolite_bl.columns[3:].to_list()

    # REVISED: Checks if every individual baseline metabolite exists in the master dataset
    if set(metabolite_list_bl).issubset(metabolite_set_long):
        column_list = ['participant_id'] + metabolite_list_bl
        df_organ = df_metabolite_long[column_list]
        df_organ.to_csv(os.path.join(output_dir, 'proteimics_' + organ + '_2_0.tsv'), index=False, sep='\t')
    else:
        print('Organ: %s' % organ)
        print(
            'UK biobank Olink at instance 2 and 3 only have ~1400 metabolites, which does not cover the entire organ-enriched metabolites... Skip...')

### It seems that the instance 2 only covers pulmonary ProtBAG

print('Stop here...')
