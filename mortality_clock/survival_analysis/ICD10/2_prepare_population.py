import os.path
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")
import argparse
parser = argparse.ArgumentParser(description="parser")

# Mandatory argument
parser.add_argument("--disease", type=str, default="/cbica/home/wenju/Dataset/UKBB/UKBB_genetic_preprocess/S4_GWAS/All_ancestry",
                           help="Path to store the clustering results")
def prepare_data(disease):

    output_tsv = '/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/SA/data/' + disease + '_diagnosis_clock.tsv'
    if not os.path.exists(output_tsv):
        ### 7 mri mortality clock
        for organ in ['brain', 'adipose', 'heart', 'kidney', 'liver', 'pancreas', 'spleen']:
            organ_big = organ.capitalize()
            tsv = os.path.join('/cbica/home/wenju/Reproducibile_paper/WholeBodyClock',
                               organ + "_mri_mortality_clock", organ + "_mri_mortality_clock_predictions.tsv")
            df_data = pd.read_csv(tsv, sep='\t')
            df_data.rename({organ + '_mri_mortality_clock_acceleration_z': organ_big + '_mri'}, axis=1, inplace=True)
            df_data = df_data[['participant_id', organ_big + '_mri']]
            if organ == 'brain':
                df_clock = df_data
            else:
                df_clock = pd.merge(df_clock, df_data, on="participant_id", how="outer")

        ### 11 proteomics mortality clock
        for organ in ['Reproductive_female', 'Pulmonary', 'Heart', 'Brain', 'Eye', 'Hepatic', 'Renal',
                      'Reproductive_male', 'Endocrine', 'Immune', 'Skin']:
            organ_small = organ.lower()
            tsv = os.path.join('/cbica/home/wenju/Reproducibile_paper/WholeBodyClock',
                               organ + "_proteomics_mortality_clock",
                               organ_small + "_proteomics_mortality_clock_predictions.tsv")
            df_data = pd.read_csv(tsv, sep='\t')
            df_data.rename({organ_small + '_proteomics_mortality_clock_acceleration_z': organ + '_proteomics'}, axis=1,
                           inplace=True)
            df_data = df_data[['participant_id', organ + '_proteomics']]
            df_clock = pd.merge(df_clock, df_data, on="participant_id", how="outer")

        ### 5 metabolomics mortality clock
        for organ in ['Endocrine', 'Digestive', 'Hepatic', 'Immune']:
            organ_small = organ.lower()
            tsv = os.path.join('/cbica/home/wenju/Reproducibile_paper/WholeBodyClock',
                               organ + "_metabolomics_mortality_clock",
                               organ_small + "_metabolomics_mortality_clock_predictions.tsv")
            df_data = pd.read_csv(tsv, sep='\t')
            df_data.rename({organ_small + '_metabolomics_mortality_clock_acceleration_z': organ + '_metabolomics'},
                           axis=1, inplace=True)
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
        date_pop = df_icd_diagnosis_date['participant_id'].to_list()
        diagnosis_pop = df_icd_diagnosis['participant_id'].to_list()
        pop_overlap = [value for value in date_pop if value in diagnosis_pop]
        df_icd_diagnosis_date = df_icd_diagnosis_date.loc[df_icd_diagnosis_date['participant_id'].isin(pop_overlap)]
        df_icd_diagnosis = df_icd_diagnosis.loc[df_icd_diagnosis['participant_id'].isin(pop_overlap)]

        ########################################################################################################################
        ### MRIBAG
        ########################################################################################################################
        print(
            "########################################################################################################################")
        ### merge with the brain population
        df_icd_brain = df_icd_diagnosis.loc[df_icd_diagnosis['participant_id'].isin(bag_population_list)]
        df_icd_data = df_icd_brain.iloc[:, 1:]

        ### relax this criteron to define CN
        has_disease = df_icd_data.iloc[:, 1:].apply(lambda row: row.isin([disease]).any(), axis=1)
        cn_boolean = ~has_disease  # True (CN) if none of the disease codes are present

        df_icd_cn = df_icd_brain[cn_boolean.values]
        df_icd_pt = df_icd_brain[~cn_boolean.values]
        df_icd_pt_pop = df_icd_pt[['participant_id']]
        df_icd_cn.replace(np.nan, -1, inplace=True)
        df_icd_cn = df_icd_cn[['participant_id', '0_0']]

        filter_pt = (df_icd_pt.apply(lambda r: r.astype('string').str.match(disease).any(), axis=1))
        df_pt = df_icd_pt_pop.copy()
        df_cn = df_icd_cn.copy()
        df_cn.rename({'0_0': disease}, axis=1, inplace=True)
        df_cn[disease] = -1
        df_pt[disease] = filter_pt
        df_pt[disease] = df_pt[disease].map({False: np.nan, True: 1})
        df_disease = pd.concat([df_cn, df_pt], ignore_index=True)
        df_disease.dropna(inplace=True)
        df_disease = df_disease.merge(df_clock, how='inner', left_on='participant_id', right_on='participant_id')
        print("There are %d CN vs. %d PT for ICD-10 disease: %s" % (
            df_disease[disease].value_counts()[-1], df_disease[disease].value_counts()[1], disease))
        #### date
        population_list = df_disease['participant_id'].to_list()
        case_list = []
        date_list = []
        for ind in population_list:
            print("check diagnosis date for this participant: %s" % ind)
            row = df_disease.loc[df_disease['participant_id'].isin([ind])]
            binary_diagnosis = row[disease].values[0]
            if binary_diagnosis == 1.0:
                row_date = df_icd_diagnosis_date.loc[df_icd_diagnosis_date['participant_id'].isin([ind])]
                if row_date.shape[0] == 1:
                    row_diagnosis = df_icd_diagnosis.loc[df_icd_diagnosis['participant_id'].isin([ind])]
                    if row_diagnosis.shape[0] == 1:
                        row_diagnosis_series = row_diagnosis.squeeze(axis=0)
                        de_index_list = row_diagnosis_series[row_diagnosis_series == disease].index
                        if len(de_index_list) > 0:
                            colum_name = row_diagnosis.apply(lambda row: row[row == disease].index, axis=1).values[0][0]
                            date = row_date[colum_name].values[0]
                            case_list.append(1)
                            date_list.append(date)
                        else:
                            df_disease = df_disease.drop(df_disease[df_disease['participant_id'] == ind].index)
                    else:
                        df_disease = df_disease.drop(df_disease[df_disease['participant_id'] == ind].index)
                else:
                    df_disease = df_disease.drop(df_disease[df_disease['participant_id'] == ind].index)
            else:
                row_date = df_icd_diagnosis_date.loc[df_icd_diagnosis_date['participant_id'].isin([ind])]
                if row_date.shape[0] == 1:
                    case_list.append(0)
                    date_list.append(np.nan)
                else:
                    df_disease = df_disease.drop(df_disease[df_disease['participant_id'] == ind].index)
        df_disease['case'] = case_list
        df_disease['date'] = date_list
        df_disease.to_csv(output_tsv, index=False, sep='\t', encoding='utf-8')

def main(options):
    prepare_data(options.disease)

if __name__ == "__main__":
    commandline = parser.parse_known_args()
    options = commandline[0]
    if commandline[1]:
        raise Exception("unknown arguments: %s" % parser.parse_known_args()[1])
    main(options)

print('Here...')