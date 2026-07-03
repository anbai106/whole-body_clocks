import pandas as pd
import os

### read the metabolite population
df_metabolite_long = pd.read_csv('/Users/hao/cubic-home/Dataset/UKBB/UKBB_NMR_metabolomics/QC_ukbnmr/nmr_biomarker_data.csv')
df_metabolite_long.rename({'eid': 'participant_id'}, axis=1, inplace=True)
df_metabolite_long = df_metabolite_long.loc[df_metabolite_long['visit_index'].isin([1])]

metabolite_list_long = df_metabolite_long.columns[2:].to_list()

# Convert the master list to a set once outside the loop for O(1) optimal speed
metabolite_set_long = set(metabolite_list_long)

output_dir = '/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/mortality_clock/longitudinal/metabolomics/data'
organ_list = ['Endocrine', 'Digestive', 'Hepatic', 'Immune']

for organ in organ_list:
    df_metabolite_bl = pd.read_csv(
        '/Users/hao/cubic-home/Reproducibile_paper/UKBB_metabolomics/MLNI/data/' + organ + '/PT/patient_pop_non_derived.tsv',
        sep='\t')
    metabolite_list_bl = df_metabolite_bl.columns[3:].to_list()

    # REVISED: Checks if every individual baseline metabolite exists in the master dataset
    if set(metabolite_list_bl).issubset(metabolite_set_long):
        column_list = ['participant_id'] + metabolite_list_bl
        # Select the columns
        df_organ = df_metabolite_long[column_list]
        # Get dimensions before dropping
        rows_before = df_organ.shape[0]
        # Drop rows where ANY value is NaN
        df_organ = df_organ.dropna()
        # Get dimensions after dropping
        rows_after = df_organ.shape[0]
        rows_dropped = rows_before - rows_after
        # Print the final dimensions and summary
        print(f"Dimensions before dropping NaNs: {rows_before} participants/rows")
        print(f"Dimensions after dropping NaNs:  {rows_after} participants/rows")
        print(f"Total participants/rows dropped:  {rows_dropped}")
        df_organ.to_csv(os.path.join(output_dir, 'metabolomics' + organ + '_1_0.tsv'), index=False, sep='\t')
    else:
        raise Exception('Organ: %s has issues...' % organ)

print('Stop here...')
