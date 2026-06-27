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

    phenotype_list = ['BAG']

    ### save the full participant population for creating the fastgwa binary file
    df_fastgwa_binary = df[['participant_id']]
    df_fastgwa_binary.rename({'participant_id': 'FID'}, axis=1, inplace=True)
    df_fastgwa_binary.insert(1, "IID", list(df_fastgwa_binary['FID']))
    df_fastgwa_binary.to_csv(os.path.join(output_dir, 'BAG_keep_for_fastgwa_binary_file.txt'), index=False, sep='\t',
                    header=False)

    ####################################################################################################################
    ## Fun part, let's create the genetic data for each organ age gap.
    for BAG in phenotype_list:
        print("Generating genotype data for fastgwa to run GWAS for BAG: %s" % BAG)
        ### pheno.txt
        df_BAG = df[['participant_id', BAG] + cov_list]
        df_BAG.dropna(inplace=True)

        print("After removing NAN values, the sample size is: %d" % df_BAG.shape[0])
        df_pheno = df_BAG[['participant_id', BAG]]
        df_pheno.rename({'participant_id': 'FID'}, axis=1, inplace=True)
        df_pheno.insert(1, "IID", list(df_pheno['FID']))
        df_pheno.to_csv(os.path.join(output_dir, 'BAG_pheno.txt'), index=False, sep='\t')

        ### cov.txt
        df_cov = df_BAG[['participant_id'] + cov_list]
        df_cov.rename({'participant_id': 'FID'}, axis=1, inplace=True)
        df_cov.insert(1, "IID", list(df_cov['FID']))
        df_cov.to_csv(os.path.join(output_dir, 'BAG_cov.txt'), index=False, sep='\t')
        ### save the unrelated individuls for fastgwa --keep argument
        df_fastgwa = df_cov[['FID', 'IID']]
        df_fastgwa.to_csv(os.path.join(output_dir, 'BAG_keep_for_fastgwa.txt'), index=False, sep='\t',
                        header=False)
    print("Stop here... ")

def create_folder_if_not_exists(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

for sex in ['female', 'male']:
    for organ in ["Metabolic"]:
        sex_capital = sex.capitalize()
        df_data = pd.read_csv(
            '/Users/hao/cubic-home/Reproducibile_paper/FemaleAgingClock/MLNI/output_dir/R1/sex_pool_vs_sex_specific_metabolic_MetBAG_4_downstream.tsv',
            sep='\t')
        df_data = df_data.loc[df_data['Approach'].isin(['Sex-pooled'])]
        df_data = df_data.loc[df_data['Sex'].isin([sex_capital])]
        df_data['BAG'] = df_data['predicted_age'] - df_data['diagnosis']
        df_data = df_data[['participant_id', 'diagnosis', 'predicted_age', 'BAG']]
        df_data['participant_id'] = df_data['participant_id'].astype(int)
        cov_tsv = '/Users/hao/cubic-home/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv'
        output_dir = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/data/R1_sex_pooled', organ + '_non_derived_' + sex)
        create_folder_if_not_exists(output_dir)
        fastgwa_fam = "/Users/hao/cubic-home/Dataset/UKBB/UKBB_genetic_preprocess_all/S3_apply_all/chr_all_AllUKBBPeople.fam"
        prepare_data(df_data, cov_tsv, output_dir, fastgwa_fam)