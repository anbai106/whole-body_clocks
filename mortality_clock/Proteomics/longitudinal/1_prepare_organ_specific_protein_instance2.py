import pandas as pd
import os

### read the protein population
df_protein_long = pd.read_csv(
    '/Users/hao/cubic-home/Reproducibile_paper/UKBB_Proteomics/imputation_2_0/UKBB_preoteomics_real_data_all_ancestry_with_related_indi_4_autocomplete_final_imputation_2_0.tsv',
    sep='\t')
df_protein_long = df_protein_long.rename(columns={"ID": "participant_id"})
protein_list_long = df_protein_long.columns[2:].to_list()

# Convert the master list to a set once outside the loop for O(1) optimal speed
protein_set_long = set(protein_list_long)

output_dir = '/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/proteomics/data'
organ_list = ['Reproductive_female', 'Pulmonary', 'Heart', 'Brain', 'Eye', 'Hepatic', 'Renal', 'Reproductive_male', 'Endocrine', 'Immune', 'Skin']

for organ in organ_list:
    df_protein_bl = pd.read_csv(
        '/Users/hao/cubic-home/Reproducibile_paper/UKBB_Proteomics/MLNI/data/' + organ + '/training/training_4589.tsv',
        sep='\t')
    protein_list_bl = df_protein_bl.columns[3:].to_list()

    # REVISED: Checks if every individual baseline protein exists in the master dataset
    if set(protein_list_bl).issubset(protein_set_long):
        column_list = ['participant_id'] + protein_list_bl
        df_organ = df_protein_long[column_list]
        df_organ.to_csv(os.path.join(output_dir, 'proteimics_' + organ + '_2_0.tsv'), index=False, sep='\t')
    else:
        print('Organ: %s' % organ)
        print(
            'UK biobank Olink at instance 2 and 3 only have ~1400 proteins, which does not cover the entire organ-enriched proteins... Skip...')

### It seems that the instance 2 only covers pulmonary ProtBAG

print('Stop here...')
