import pandas as pd
import os

def prepare_data(df_data, cov_tsv, output_dir, fastgwa_fam):
    """
    This is a function to include unrelated individuals after running King software, and merge it with the available imaging data
    Note: the genetic PCAs are here: /cbica/projects/ISTAGING/Pipelines/ClinicalDataConsolidation_201911/Data/External_Data/UKBiobank/ukb_200622.csv
    fields: PC0_{#}
    Returns:
        NOTE: Forget the KING unrelated files, but let's not redo things...
    """
    ## fam file
    df_fam = pd.read_csv(fastgwa_fam, delimiter=' ', header=None)[[0, 1]]
    df_fam.columns = ['FID', 'participant_id']
    print(df_fam['FID'].equals(df_fam['participant_id']))
    del df_fam['FID']

    df_cov = pd.read_csv(cov_tsv)
    df_cov.rename({'eid': 'participant_id'}, axis=1, inplace=True)
    df_cov = df_cov.loc[df_cov['genetic_ethnic_grouping_f22006_0_0'].isin([1])]
    cov_list = ['sex_f31_0_0', 'age_when_attended_assessment_centre_f21003_0_0', 'genetic_principal_components_f22009_0_1', 'genetic_principal_components_f22009_0_2',
                'genetic_principal_components_f22009_0_3', 'genetic_principal_components_f22009_0_4', 'genetic_principal_components_f22009_0_5', 'genetic_principal_components_f22009_0_6',
                'genetic_principal_components_f22009_0_7', 'genetic_principal_components_f22009_0_8', 'genetic_principal_components_f22009_0_9', 'genetic_principal_components_f22009_0_10',
                'genetic_principal_components_f22009_0_11', 'genetic_principal_components_f22009_0_12', 'genetic_principal_components_f22009_0_13', 'genetic_principal_components_f22009_0_14',
                'genetic_principal_components_f22009_0_15', 'genetic_principal_components_f22009_0_16', 'genetic_principal_components_f22009_0_17', 'genetic_principal_components_f22009_0_18',
                'genetic_principal_components_f22009_0_19', 'genetic_principal_components_f22009_0_20', 'genetic_principal_components_f22009_0_21', 'genetic_principal_components_f22009_0_22',
                'genetic_principal_components_f22009_0_23', 'genetic_principal_components_f22009_0_24', 'genetic_principal_components_f22009_0_25', 'genetic_principal_components_f22009_0_26',
                'genetic_principal_components_f22009_0_27', 'genetic_principal_components_f22009_0_28', 'genetic_principal_components_f22009_0_29', 'genetic_principal_components_f22009_0_30',
                'genetic_principal_components_f22009_0_31', 'genetic_principal_components_f22009_0_32', 'genetic_principal_components_f22009_0_33', 'genetic_principal_components_f22009_0_34',
                'genetic_principal_components_f22009_0_35', 'genetic_principal_components_f22009_0_36', 'genetic_principal_components_f22009_0_37', 'genetic_principal_components_f22009_0_38',
                'genetic_principal_components_f22009_0_39', 'genetic_principal_components_f22009_0_40', 'weight_f21002_0_0', 'standing_height_f50_0_0',
                'waist_circumference_f48_0_0', 'body_mass_index_bmi_f23104_0_0',  'diastolic_blood_pressure_automated_reading_f4079_0_0',
                'systolic_blood_pressure_automated_reading_f4080_0_0']
    df_cov = df_cov[['participant_id'] + cov_list]

    df = df_data.merge(df_fam, how='inner', left_on='participant_id', right_on='participant_id')
    df = df.merge(df_cov, how='inner', left_on='participant_id', right_on='participant_id')

    phenotype_list = ['EPOCH']

    ### save the full participant population for creating the fastgwa binary file
    df_fastgwa_binary = df[['participant_id']]
    df_fastgwa_binary.rename({'participant_id': 'FID'}, axis=1, inplace=True)
    df_fastgwa_binary.insert(1, "IID", list(df_fastgwa_binary['FID']))
    df_fastgwa_binary.to_csv(os.path.join(output_dir, 'EPOCH_keep_for_fastgwa_binary_file.txt'), index=False, sep='\t',
                    header=False)

    ####################################################################################################################
    ## Fun part, let's create the genetic data for each organ age gap.
    for EPOCH in phenotype_list:
        print("Generating genotype data for fastgwa to run GWAS for EPOCH: %s" % EPOCH)
        ### pheno.txt
        df_EPOCH = df[['participant_id', EPOCH] + cov_list]
        df_EPOCH.dropna(inplace=True)

        print("After removing NAN values, the sample size is: %d" % df_EPOCH.shape[0])
        df_pheno = df_EPOCH[['participant_id', EPOCH]]
        df_pheno.rename({'participant_id': 'FID'}, axis=1, inplace=True)
        df_pheno.insert(1, "IID", list(df_pheno['FID']))
        df_pheno.to_csv(os.path.join(output_dir, 'EPOCH_pheno.txt'), index=False, sep='\t')

        ### cov.txt
        df_cov = df_EPOCH[['participant_id'] + cov_list]
        df_cov.rename({'participant_id': 'FID'}, axis=1, inplace=True)
        df_cov.insert(1, "IID", list(df_cov['FID']))
        df_cov.to_csv(os.path.join(output_dir, 'EPOCH_cov.txt'), index=False, sep='\t')
        ### save the unrelated individuls for fastgwa --keep argument
        df_fastgwa = df_cov[['FID', 'IID']]
        df_fastgwa.to_csv(os.path.join(output_dir, 'EPOCH_keep_for_fastgwa.txt'), index=False, sep='\t',
                        header=False)
    print("Stop here... ")

def create_folder_if_not_exists(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

disease = 'dementia'

tsv = '/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/all_disease_lepoch_incremental_value_scale_qc/stable_significant_disease_clock_acceleration_z_wide.tsv'
df_data = pd.read_csv(tsv, sep='\t')

### 7 proteomics  dementia clock
for organ in ['Brain', 'Hepatic', 'Endocrine', 'Heart', 'Immune', 'Reproductive_female', 'Reproductive_male']:
    organ_small = organ.lower()
    df_data.rename({disease + '_' + organ_small + '_proteomics_' + disease + '_clock_acceleration_z': 'EPOCH'}, axis=1, inplace=True)
    cov_tsv = '/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv'
    output_dir = os.path.join("/cbica/home/wenju/Reproducibile_paper/WholeBodyClock", organ + '_proteomics_' + disease + '_clock', 'fastGWA', 'data')
    create_folder_if_not_exists(output_dir)
    fastgwa_fam = "/cbica/home/wenju/Dataset/UKBB/UKBB_genetic_preprocess_all/S3_apply_all/chr_all_AllUKBBPeople.fam"
    prepare_data(df_data, cov_tsv, output_dir, fastgwa_fam)

### 4 metabolomics  dementia clock
for organ in ['Endocrine', 'Digestive', 'Hepatic', 'Metabolic']:
    organ_small = organ.lower()
    df_data.rename({disease + '_' + organ_small + '_metabolomics_' + disease + '_clock_acceleration_z': 'EPOCH'}, axis=1, inplace=True)
    cov_tsv = '/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv'
    output_dir = os.path.join("/cbica/home/wenju/Reproducibile_paper/WholeBodyClock", organ + '_metabolomics_' + disease + '_clock', 'fastGWA', 'data')
    create_folder_if_not_exists(output_dir)
    fastgwa_fam = "/cbica/home/wenju/Dataset/UKBB/UKBB_genetic_preprocess_all/S3_apply_all/chr_all_AllUKBBPeople.fam"
    prepare_data(df_data, cov_tsv, output_dir, fastgwa_fam)