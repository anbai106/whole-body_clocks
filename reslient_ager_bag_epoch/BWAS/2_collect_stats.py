import numpy as np
import pandas as pd
import os

########################################################################################################################
### 119 Brain GM MUSE
########################################################################################################################
idp_list_gm = [
    "MUSE_Volume_23", "MUSE_Volume_30", "MUSE_Volume_31", "MUSE_Volume_32", "MUSE_Volume_36", "MUSE_Volume_37",
    "MUSE_Volume_38", "MUSE_Volume_39", "MUSE_Volume_47", "MUSE_Volume_48", "MUSE_Volume_55", "MUSE_Volume_56",
    "MUSE_Volume_57", "MUSE_Volume_58", "MUSE_Volume_59", "MUSE_Volume_60", "MUSE_Volume_71", "MUSE_Volume_72",
    "MUSE_Volume_73", "MUSE_Volume_75", "MUSE_Volume_76", "MUSE_Volume_100", "MUSE_Volume_101", "MUSE_Volume_102",
    "MUSE_Volume_103", "MUSE_Volume_104", "MUSE_Volume_105", "MUSE_Volume_106", "MUSE_Volume_107", "MUSE_Volume_108",
    "MUSE_Volume_109", "MUSE_Volume_112", "MUSE_Volume_113", "MUSE_Volume_114", "MUSE_Volume_115", "MUSE_Volume_116",
    "MUSE_Volume_117", "MUSE_Volume_118", "MUSE_Volume_119", "MUSE_Volume_120", "MUSE_Volume_121", "MUSE_Volume_122",
    "MUSE_Volume_123", "MUSE_Volume_124", "MUSE_Volume_125", "MUSE_Volume_128", "MUSE_Volume_129", "MUSE_Volume_132",
    "MUSE_Volume_133", "MUSE_Volume_134", "MUSE_Volume_135", "MUSE_Volume_136", "MUSE_Volume_137", "MUSE_Volume_138",
    "MUSE_Volume_139", "MUSE_Volume_140", "MUSE_Volume_141", "MUSE_Volume_142", "MUSE_Volume_143", "MUSE_Volume_144",
    "MUSE_Volume_145", "MUSE_Volume_146", "MUSE_Volume_147", "MUSE_Volume_148", "MUSE_Volume_149", "MUSE_Volume_150",
    "MUSE_Volume_151", "MUSE_Volume_152", "MUSE_Volume_153", "MUSE_Volume_154", "MUSE_Volume_155", "MUSE_Volume_156",
    "MUSE_Volume_157", "MUSE_Volume_160", "MUSE_Volume_161", "MUSE_Volume_162", "MUSE_Volume_163", "MUSE_Volume_164",
    "MUSE_Volume_165", "MUSE_Volume_166", "MUSE_Volume_167", "MUSE_Volume_168", "MUSE_Volume_169", "MUSE_Volume_170",
    "MUSE_Volume_171", "MUSE_Volume_172", "MUSE_Volume_173", "MUSE_Volume_174", "MUSE_Volume_175", "MUSE_Volume_176",
    "MUSE_Volume_177", "MUSE_Volume_178", "MUSE_Volume_179", "MUSE_Volume_180", "MUSE_Volume_181", "MUSE_Volume_182",
    "MUSE_Volume_183", "MUSE_Volume_184", "MUSE_Volume_185", "MUSE_Volume_186", "MUSE_Volume_187", "MUSE_Volume_190",
    "MUSE_Volume_191", "MUSE_Volume_192", "MUSE_Volume_193", "MUSE_Volume_194", "MUSE_Volume_195", "MUSE_Volume_196",
    "MUSE_Volume_197", "MUSE_Volume_198", "MUSE_Volume_199", "MUSE_Volume_200", "MUSE_Volume_201", "MUSE_Volume_202",
    "MUSE_Volume_203", "MUSE_Volume_204", "MUSE_Volume_205", "MUSE_Volume_206", "MUSE_Volume_207"
]
df_brain_gm = pd.DataFrame(columns=['BAG', 'IDP', 'Log_Odds', 'SE', 'Z', 'P_value', 'OR', 'OR_CI_Lower', 'OR_CI_Upper'])
i = 0
for idp in idp_list_gm:
    tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/BWAS/Brain/T1', 'NSS_logistic_results_' + f'{idp}' + '.tsv')
    i += 1
    if os.path.exists(tsv):
        print(f'Collect: {tsv}')
        df = pd.read_csv(tsv, sep='\t')
        
        if i == 1:
            df_brain_gm = df
        else:
            df_brain_gm = pd.concat([df_brain_gm, df], ignore_index=True)
    else:
        print(f'Model does not converge: {tsv}')

df_brain_gm['Organ'] = 'Brain-GM'
df_brain_gm.to_csv('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_Brain_GM.tsv', index=False, sep='\t')

########################################################################################################################
### 192 WM DTI
########################################################################################################################
idp_list_wm = [
    "mean_fa_in_middle_cerebellar_peduncle_on_fa_skeleton_f25056_2_0",
    "mean_fa_in_pontine_crossing_tract_on_fa_skeleton_f25057_2_0",
    "mean_fa_in_genu_of_corpus_callosum_on_fa_skeleton_f25058_2_0",
    "mean_fa_in_body_of_corpus_callosum_on_fa_skeleton_f25059_2_0",
    "mean_fa_in_splenium_of_corpus_callosum_on_fa_skeleton_f25060_2_0",
    "mean_fa_in_fornix_on_fa_skeleton_f25061_2_0",
    "mean_fa_in_corticospinal_tract_on_fa_skeleton_right_f25062_2_0",
    "mean_fa_in_corticospinal_tract_on_fa_skeleton_left_f25063_2_0",
    "mean_fa_in_medial_lemniscus_on_fa_skeleton_right_f25064_2_0",
    "mean_fa_in_medial_lemniscus_on_fa_skeleton_left_f25065_2_0",
    "mean_fa_in_inferior_cerebellar_peduncle_on_fa_skeleton_right_f25066_2_0",
    "mean_fa_in_inferior_cerebellar_peduncle_on_fa_skeleton_left_f25067_2_0",
    "mean_fa_in_superior_cerebellar_peduncle_on_fa_skeleton_right_f25068_2_0",
    "mean_fa_in_superior_cerebellar_peduncle_on_fa_skeleton_left_f25069_2_0",
    "mean_fa_in_cerebral_peduncle_on_fa_skeleton_right_f25070_2_0",
    "mean_fa_in_cerebral_peduncle_on_fa_skeleton_left_f25071_2_0",
    "mean_fa_in_anterior_limb_of_internal_capsule_on_fa_skeleton_right_f25072_2_0",
    "mean_fa_in_anterior_limb_of_internal_capsule_on_fa_skeleton_left_f25073_2_0",
    "mean_fa_in_posterior_limb_of_internal_capsule_on_fa_skeleton_right_f25074_2_0",
    "mean_fa_in_posterior_limb_of_internal_capsule_on_fa_skeleton_left_f25075_2_0",
    "mean_fa_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_right_f25076_2_0",
    "mean_fa_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_left_f25077_2_0",
    "mean_fa_in_anterior_corona_radiata_on_fa_skeleton_right_f25078_2_0",
    "mean_fa_in_anterior_corona_radiata_on_fa_skeleton_left_f25079_2_0",
    "mean_fa_in_superior_corona_radiata_on_fa_skeleton_right_f25080_2_0",
    "mean_fa_in_superior_corona_radiata_on_fa_skeleton_left_f25081_2_0",
    "mean_fa_in_posterior_corona_radiata_on_fa_skeleton_right_f25082_2_0",
    "mean_fa_in_posterior_corona_radiata_on_fa_skeleton_left_f25083_2_0",
    "mean_fa_in_posterior_thalamic_radiation_on_fa_skeleton_right_f25084_2_0",
    "mean_fa_in_posterior_thalamic_radiation_on_fa_skeleton_left_f25085_2_0",
    "mean_fa_in_sagittal_stratum_on_fa_skeleton_right_f25086_2_0",
    "mean_fa_in_sagittal_stratum_on_fa_skeleton_left_f25087_2_0",
    "mean_fa_in_external_capsule_on_fa_skeleton_right_f25088_2_0",
    "mean_fa_in_external_capsule_on_fa_skeleton_left_f25089_2_0",
    "mean_fa_in_cingulum_cingulate_gyrus_on_fa_skeleton_right_f25090_2_0",
    "mean_fa_in_cingulum_cingulate_gyrus_on_fa_skeleton_left_f25091_2_0",
    "mean_fa_in_cingulum_hippocampus_on_fa_skeleton_right_f25092_2_0",
    "mean_fa_in_cingulum_hippocampus_on_fa_skeleton_left_f25093_2_0",
    "mean_fa_in_fornix_cresstria_terminalis_on_fa_skeleton_right_f25094_2_0",
    "mean_fa_in_fornix_cresstria_terminalis_on_fa_skeleton_left_f25095_2_0",
    "mean_fa_in_superior_longitudinal_fasciculus_on_fa_skeleton_right_f25096_2_0",
    "mean_fa_in_superior_longitudinal_fasciculus_on_fa_skeleton_left_f25097_2_0",
    "mean_fa_in_superior_frontooccipital_fasciculus_on_fa_skeleton_right_f25098_2_0",
    "mean_fa_in_superior_frontooccipital_fasciculus_on_fa_skeleton_left_f25099_2_0",
    "mean_fa_in_uncinate_fasciculus_on_fa_skeleton_right_f25100_2_0",
    "mean_fa_in_uncinate_fasciculus_on_fa_skeleton_left_f25101_2_0",
    "mean_fa_in_tapetum_on_fa_skeleton_right_f25102_2_0",
    "mean_fa_in_tapetum_on_fa_skeleton_left_f25103_2_0",
    "mean_md_in_middle_cerebellar_peduncle_on_fa_skeleton_f25104_2_0",
    "mean_md_in_pontine_crossing_tract_on_fa_skeleton_f25105_2_0",
    "mean_md_in_genu_of_corpus_callosum_on_fa_skeleton_f25106_2_0",
    "mean_md_in_body_of_corpus_callosum_on_fa_skeleton_f25107_2_0",
    "mean_md_in_splenium_of_corpus_callosum_on_fa_skeleton_f25108_2_0",
    "mean_md_in_fornix_on_fa_skeleton_f25109_2_0",
    "mean_md_in_corticospinal_tract_on_fa_skeleton_right_f25110_2_0",
    "mean_md_in_corticospinal_tract_on_fa_skeleton_left_f25111_2_0",
    "mean_md_in_medial_lemniscus_on_fa_skeleton_right_f25112_2_0",
    "mean_md_in_medial_lemniscus_on_fa_skeleton_left_f25113_2_0",
    "mean_md_in_inferior_cerebellar_peduncle_on_fa_skeleton_right_f25114_2_0",
    "mean_md_in_inferior_cerebellar_peduncle_on_fa_skeleton_left_f25115_2_0",
    "mean_md_in_superior_cerebellar_peduncle_on_fa_skeleton_right_f25116_2_0",
    "mean_md_in_superior_cerebellar_peduncle_on_fa_skeleton_left_f25117_2_0",
    "mean_md_in_cerebral_peduncle_on_fa_skeleton_right_f25118_2_0",
    "mean_md_in_cerebral_peduncle_on_fa_skeleton_left_f25119_2_0",
    "mean_md_in_anterior_limb_of_internal_capsule_on_fa_skeleton_right_f25120_2_0",
    "mean_md_in_anterior_limb_of_internal_capsule_on_fa_skeleton_left_f25121_2_0",
    "mean_md_in_posterior_limb_of_internal_capsule_on_fa_skeleton_right_f25122_2_0",
    "mean_md_in_posterior_limb_of_internal_capsule_on_fa_skeleton_left_f25123_2_0",
    "mean_md_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_right_f25124_2_0",
    "mean_md_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_left_f25125_2_0",
    "mean_md_in_anterior_corona_radiata_on_fa_skeleton_right_f25126_2_0",
    "mean_md_in_anterior_corona_radiata_on_fa_skeleton_left_f25127_2_0",
    "mean_md_in_superior_corona_radiata_on_fa_skeleton_right_f25128_2_0",
    "mean_md_in_superior_corona_radiata_on_fa_skeleton_left_f25129_2_0",
    "mean_md_in_posterior_corona_radiata_on_fa_skeleton_right_f25130_2_0",
    "mean_md_in_posterior_corona_radiata_on_fa_skeleton_left_f25131_2_0",
    "mean_md_in_posterior_thalamic_radiation_on_fa_skeleton_right_f25132_2_0",
    "mean_md_in_posterior_thalamic_radiation_on_fa_skeleton_left_f25133_2_0",
    "mean_md_in_sagittal_stratum_on_fa_skeleton_right_f25134_2_0",
    "mean_md_in_sagittal_stratum_on_fa_skeleton_left_f25135_2_0",
    "mean_md_in_external_capsule_on_fa_skeleton_right_f25136_2_0",
    "mean_md_in_external_capsule_on_fa_skeleton_left_f25137_2_0",
    "mean_md_in_cingulum_cingulate_gyrus_on_fa_skeleton_right_f25138_2_0",
    "mean_md_in_cingulum_cingulate_gyrus_on_fa_skeleton_left_f25139_2_0",
    "mean_md_in_cingulum_hippocampus_on_fa_skeleton_right_f25140_2_0",
    "mean_md_in_cingulum_hippocampus_on_fa_skeleton_left_f25141_2_0",
    "mean_md_in_fornix_cresstria_terminalis_on_fa_skeleton_right_f25142_2_0",
    "mean_md_in_fornix_cresstria_terminalis_on_fa_skeleton_left_f25143_2_0",
    "mean_md_in_superior_longitudinal_fasciculus_on_fa_skeleton_right_f25144_2_0",
    "mean_md_in_superior_longitudinal_fasciculus_on_fa_skeleton_left_f25145_2_0",
    "mean_md_in_superior_frontooccipital_fasciculus_on_fa_skeleton_right_f25146_2_0",
    "mean_md_in_superior_frontooccipital_fasciculus_on_fa_skeleton_left_f25147_2_0",
    "mean_md_in_uncinate_fasciculus_on_fa_skeleton_right_f25148_2_0",
    "mean_md_in_uncinate_fasciculus_on_fa_skeleton_left_f25149_2_0",
    "mean_md_in_tapetum_on_fa_skeleton_right_f25150_2_0",
    "mean_md_in_tapetum_on_fa_skeleton_left_f25151_2_0",
    "mean_icvf_in_middle_cerebellar_peduncle_on_fa_skeleton_f25344_2_0",
    "mean_icvf_in_pontine_crossing_tract_on_fa_skeleton_f25345_2_0",
    "mean_icvf_in_genu_of_corpus_callosum_on_fa_skeleton_f25346_2_0",
    "mean_icvf_in_body_of_corpus_callosum_on_fa_skeleton_f25347_2_0",
    "mean_icvf_in_splenium_of_corpus_callosum_on_fa_skeleton_f25348_2_0",
    "mean_icvf_in_fornix_on_fa_skeleton_f25349_2_0",
    "mean_icvf_in_corticospinal_tract_on_fa_skeleton_right_f25350_2_0",
    "mean_icvf_in_corticospinal_tract_on_fa_skeleton_left_f25351_2_0",
    "mean_icvf_in_medial_lemniscus_on_fa_skeleton_right_f25352_2_0",
    "mean_icvf_in_medial_lemniscus_on_fa_skeleton_left_f25353_2_0",
    "mean_icvf_in_inferior_cerebellar_peduncle_on_fa_skeleton_right_f25354_2_0",
    "mean_icvf_in_inferior_cerebellar_peduncle_on_fa_skeleton_left_f25355_2_0",
    "mean_icvf_in_superior_cerebellar_peduncle_on_fa_skeleton_right_f25356_2_0",
    "mean_icvf_in_superior_cerebellar_peduncle_on_fa_skeleton_left_f25357_2_0",
    "mean_icvf_in_cerebral_peduncle_on_fa_skeleton_right_f25358_2_0",
    "mean_icvf_in_cerebral_peduncle_on_fa_skeleton_left_f25359_2_0",
    "mean_icvf_in_anterior_limb_of_internal_capsule_on_fa_skeleton_right_f25360_2_0",
    "mean_icvf_in_anterior_limb_of_internal_capsule_on_fa_skeleton_left_f25361_2_0",
    "mean_icvf_in_posterior_limb_of_internal_capsule_on_fa_skeleton_right_f25362_2_0",
    "mean_icvf_in_posterior_limb_of_internal_capsule_on_fa_skeleton_left_f25363_2_0",
    "mean_icvf_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_right_f25364_2_0",
    "mean_icvf_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_left_f25365_2_0",
    "mean_icvf_in_anterior_corona_radiata_on_fa_skeleton_right_f25366_2_0",
    "mean_icvf_in_anterior_corona_radiata_on_fa_skeleton_left_f25367_2_0",
    "mean_icvf_in_superior_corona_radiata_on_fa_skeleton_right_f25368_2_0",
    "mean_icvf_in_superior_corona_radiata_on_fa_skeleton_left_f25369_2_0",
    "mean_icvf_in_posterior_corona_radiata_on_fa_skeleton_right_f25370_2_0",
    "mean_icvf_in_posterior_corona_radiata_on_fa_skeleton_left_f25371_2_0",
    "mean_icvf_in_posterior_thalamic_radiation_on_fa_skeleton_right_f25372_2_0",
    "mean_icvf_in_posterior_thalamic_radiation_on_fa_skeleton_left_f25373_2_0",
    "mean_icvf_in_sagittal_stratum_on_fa_skeleton_right_f25374_2_0",
    "mean_icvf_in_sagittal_stratum_on_fa_skeleton_left_f25375_2_0",
    "mean_icvf_in_external_capsule_on_fa_skeleton_right_f25376_2_0",
    "mean_icvf_in_external_capsule_on_fa_skeleton_left_f25377_2_0",
    "mean_icvf_in_cingulum_cingulate_gyrus_on_fa_skeleton_right_f25378_2_0",
    "mean_icvf_in_cingulum_cingulate_gyrus_on_fa_skeleton_left_f25379_2_0",
    "mean_icvf_in_cingulum_hippocampus_on_fa_skeleton_right_f25380_2_0",
    "mean_icvf_in_cingulum_hippocampus_on_fa_skeleton_left_f25381_2_0",
    "mean_icvf_in_fornix_cresstria_terminalis_on_fa_skeleton_right_f25382_2_0",
    "mean_icvf_in_fornix_cresstria_terminalis_on_fa_skeleton_left_f25383_2_0",
    "mean_icvf_in_superior_longitudinal_fasciculus_on_fa_skeleton_right_f25384_2_0",
    "mean_icvf_in_superior_longitudinal_fasciculus_on_fa_skeleton_left_f25385_2_0",
    "mean_icvf_in_superior_frontooccipital_fasciculus_on_fa_skeleton_right_f25386_2_0",
    "mean_icvf_in_superior_frontooccipital_fasciculus_on_fa_skeleton_left_f25387_2_0",
    "mean_icvf_in_uncinate_fasciculus_on_fa_skeleton_right_f25388_2_0",
    "mean_icvf_in_uncinate_fasciculus_on_fa_skeleton_left_f25389_2_0",
    "mean_icvf_in_tapetum_on_fa_skeleton_right_f25390_2_0",
    "mean_icvf_in_tapetum_on_fa_skeleton_left_f25391_2_0",
    "mean_od_in_middle_cerebellar_peduncle_on_fa_skeleton_f25392_2_0",
    "mean_od_in_pontine_crossing_tract_on_fa_skeleton_f25393_2_0",
    "mean_od_in_genu_of_corpus_callosum_on_fa_skeleton_f25394_2_0",
    "mean_od_in_body_of_corpus_callosum_on_fa_skeleton_f25395_2_0",
    "mean_od_in_splenium_of_corpus_callosum_on_fa_skeleton_f25396_2_0",
    "mean_od_in_fornix_on_fa_skeleton_f25397_2_0",
    "mean_od_in_corticospinal_tract_on_fa_skeleton_right_f25398_2_0",
    "mean_od_in_corticospinal_tract_on_fa_skeleton_left_f25399_2_0",
    "mean_od_in_medial_lemniscus_on_fa_skeleton_right_f25400_2_0",
    "mean_od_in_medial_lemniscus_on_fa_skeleton_left_f25401_2_0",
    "mean_od_in_inferior_cerebellar_peduncle_on_fa_skeleton_right_f25402_2_0",
    "mean_od_in_inferior_cerebellar_peduncle_on_fa_skeleton_left_f25403_2_0",
    "mean_od_in_superior_cerebellar_peduncle_on_fa_skeleton_right_f25404_2_0",
    "mean_od_in_superior_cerebellar_peduncle_on_fa_skeleton_left_f25405_2_0",
    "mean_od_in_cerebral_peduncle_on_fa_skeleton_right_f25406_2_0",
    "mean_od_in_cerebral_peduncle_on_fa_skeleton_left_f25407_2_0",
    "mean_od_in_anterior_limb_of_internal_capsule_on_fa_skeleton_right_f25408_2_0",
    "mean_od_in_anterior_limb_of_internal_capsule_on_fa_skeleton_left_f25409_2_0",
    "mean_od_in_posterior_limb_of_internal_capsule_on_fa_skeleton_right_f25410_2_0",
    "mean_od_in_posterior_limb_of_internal_capsule_on_fa_skeleton_left_f25411_2_0",
    "mean_od_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_right_f25412_2_0",
    "mean_od_in_retrolenticular_part_of_internal_capsule_on_fa_skeleton_left_f25413_2_0",
    "mean_od_in_anterior_corona_radiata_on_fa_skeleton_right_f25414_2_0",
    "mean_od_in_anterior_corona_radiata_on_fa_skeleton_left_f25415_2_0",
    "mean_od_in_superior_corona_radiata_on_fa_skeleton_right_f25416_2_0",
    "mean_od_in_superior_corona_radiata_on_fa_skeleton_left_f25417_2_0",
    "mean_od_in_posterior_corona_radiata_on_fa_skeleton_right_f25418_2_0",
    "mean_od_in_posterior_corona_radiata_on_fa_skeleton_left_f25419_2_0",
    "mean_od_in_posterior_thalamic_radiation_on_fa_skeleton_right_f25420_2_0",
    "mean_od_in_posterior_thalamic_radiation_on_fa_skeleton_left_f25421_2_0",
    "mean_od_in_sagittal_stratum_on_fa_skeleton_right_f25422_2_0",
    "mean_od_in_sagittal_stratum_on_fa_skeleton_left_f25423_2_0",
    "mean_od_in_external_capsule_on_fa_skeleton_right_f25424_2_0",
    "mean_od_in_external_capsule_on_fa_skeleton_left_f25425_2_0",
    "mean_od_in_cingulum_cingulate_gyrus_on_fa_skeleton_right_f25426_2_0",
    "mean_od_in_cingulum_cingulate_gyrus_on_fa_skeleton_left_f25427_2_0",
    "mean_od_in_cingulum_hippocampus_on_fa_skeleton_right_f25428_2_0",
    "mean_od_in_cingulum_hippocampus_on_fa_skeleton_left_f25429_2_0",
    "mean_od_in_fornix_cresstria_terminalis_on_fa_skeleton_right_f25430_2_0",
    "mean_od_in_fornix_cresstria_terminalis_on_fa_skeleton_left_f25431_2_0",
    "mean_od_in_superior_longitudinal_fasciculus_on_fa_skeleton_right_f25432_2_0",
    "mean_od_in_superior_longitudinal_fasciculus_on_fa_skeleton_left_f25433_2_0",
    "mean_od_in_superior_frontooccipital_fasciculus_on_fa_skeleton_right_f25434_2_0",
    "mean_od_in_superior_frontooccipital_fasciculus_on_fa_skeleton_left_f25435_2_0",
    "mean_od_in_uncinate_fasciculus_on_fa_skeleton_right_f25436_2_0",
    "mean_od_in_uncinate_fasciculus_on_fa_skeleton_left_f25437_2_0",
    "mean_od_in_tapetum_on_fa_skeleton_right_f25438_2_0",
    "mean_od_in_tapetum_on_fa_skeleton_left_f25439_2_0"
]
df_brain_wm = pd.DataFrame(columns=['BAG', 'IDP', 'Log_Odds', 'SE', 'Z', 'P_value', 'OR', 'OR_CI_Lower', 'OR_CI_Upper'])
i = 0
for idp in idp_list_wm:
    tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/BWAS/Brain/DTI', 'NSS_logistic_results_' + f'{idp}' + '.tsv')
    i += 1
    if os.path.exists(tsv):
        print(f'Collect: {tsv}')
        df = pd.read_csv(tsv, sep='\t')
        
        if i == 1:
            df_brain_wm = df
        else:
            df_brain_wm = pd.concat([df_brain_wm, df], ignore_index=True)
    else:
        print(f'Model does not converge: {tsv}')
df_brain_wm['Organ'] = 'Brain-WM'
df_brain_wm.to_csv('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_Brain_WM.tsv', index=False, sep='\t')

########################################################################################################################
### 210 FC
########################################################################################################################
idp_list_fc = [
    "f_1", "f_2", "f_3", "f_4", "f_5", "f_6", "f_7", "f_8", "f_9", "f_10",
    "f_11", "f_12", "f_13", "f_14", "f_15", "f_16", "f_17", "f_18", "f_19", "f_20",
    "f_21", "f_22", "f_23", "f_24", "f_25", "f_26", "f_27", "f_28", "f_29", "f_30",
    "f_31", "f_32", "f_33", "f_34", "f_35", "f_36", "f_37", "f_38", "f_39", "f_40",
    "f_41", "f_42", "f_43", "f_44", "f_45", "f_46", "f_47", "f_48", "f_49", "f_50",
    "f_51", "f_52", "f_53", "f_54", "f_55", "f_56", "f_57", "f_58", "f_59", "f_60",
    "f_61", "f_62", "f_63", "f_64", "f_65", "f_66", "f_67", "f_68", "f_69", "f_70",
    "f_71", "f_72", "f_73", "f_74", "f_75", "f_76", "f_77", "f_78", "f_79", "f_80",
    "f_81", "f_82", "f_83", "f_84", "f_85", "f_86", "f_87", "f_88", "f_89", "f_90",
    "f_91", "f_92", "f_93", "f_94", "f_95", "f_96", "f_97", "f_98", "f_99", "f_100",
    "f_101", "f_102", "f_103", "f_104", "f_105", "f_106", "f_107", "f_108", "f_109", "f_110",
    "f_111", "f_112", "f_113", "f_114", "f_115", "f_116", "f_117", "f_118", "f_119", "f_120",
    "f_121", "f_122", "f_123", "f_124", "f_125", "f_126", "f_127", "f_128", "f_129", "f_130",
    "f_131", "f_132", "f_133", "f_134", "f_135", "f_136", "f_137", "f_138", "f_139", "f_140",
    "f_141", "f_142", "f_143", "f_144", "f_145", "f_146", "f_147", "f_148", "f_149", "f_150",
    "f_151", "f_152", "f_153", "f_154", "f_155", "f_156", "f_157", "f_158", "f_159", "f_160",
    "f_161", "f_162", "f_163", "f_164", "f_165", "f_166", "f_167", "f_168", "f_169", "f_170",
    "f_171", "f_172", "f_173", "f_174", "f_175", "f_176", "f_177", "f_178", "f_179", "f_180",
    "f_181", "f_182", "f_183", "f_184", "f_185", "f_186", "f_187", "f_188", "f_189", "f_190",
    "f_191", "f_192", "f_193", "f_194", "f_195", "f_196", "f_197", "f_198", "f_199", "f_200",
    "f_201", "f_202", "f_203", "f_204", "f_205", "f_206", "f_207", "f_208", "f_209", "f_210"
]
df_brain_fc = pd.DataFrame(columns=['BAG', 'IDP', 'Log_Odds', 'SE', 'Z', 'P_value', 'OR', 'OR_CI_Lower', 'OR_CI_Upper'])
i = 0
for idp in idp_list_fc:
    tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/BWAS/Brain/FC', 'NSS_logistic_results_' + f'{idp}' + '.tsv')
    i += 1
    if os.path.exists(tsv):
        print(f'Collect: {tsv}')
        df = pd.read_csv(tsv, sep='\t')
        
        if i == 1:
            df_brain_fc = df
        else:
            df_brain_fc = pd.concat([df_brain_fc, df], ignore_index=True)
    else:
        print(f'Model does not converge: {tsv}')
df_brain_fc['Organ'] = 'Brain-FC'
df_brain_fc.to_csv('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_Brain_FC.tsv', index=False, sep='\t')

########################################################################################################################
### 82 Heart
########################################################################################################################
idp_list_heart = [
    "lv_end_diastolic_volume_f24100_2_0",
    "lv_end_systolic_volume_f24101_2_0",
    "lv_stroke_volume_f24102_2_0",
    "lv_ejection_fraction_f24103_2_0",
    "lv_cardiac_output_f24104_2_0",
    "lv_myocardial_mass_f24105_2_0",
    "rv_end_diastolic_volume_f24106_2_0",
    "rv_end_systolic_volume_f24107_2_0",
    "rv_stroke_volume_f24108_2_0",
    "rv_ejection_fraction_f24109_2_0",
    "la_maximum_volume_f24110_2_0",
    "la_minimum_volume_f24111_2_0",
    "la_stroke_volume_f24112_2_0",
    "la_ejection_fraction_f24113_2_0",
    "ra_maximum_volume_f24114_2_0",
    "ra_minimum_volume_f24115_2_0",
    "ra_stroke_volume_f24116_2_0",
    "ra_ejection_fraction_f24117_2_0",
    "ascending_aorta_maximum_area_f24118_2_0",
    "ascending_aorta_minimum_area_f24119_2_0",
    "ascending_aorta_distensibility_f24120_2_0",
    "descending_aorta_maximum_area_f24121_2_0",
    "descending_aorta_minimum_area_f24122_2_0",
    "descending_aorta_distensibility_f24123_2_0",
    "lv_mean_myocardial_wall_thickness_aha_1_f24124_2_0",
    "lv_mean_myocardial_wall_thickness_aha_2_f24125_2_0",
    "lv_mean_myocardial_wall_thickness_aha_3_f24126_2_0",
    "lv_mean_myocardial_wall_thickness_aha_4_f24127_2_0",
    "lv_mean_myocardial_wall_thickness_aha_5_f24128_2_0",
    "lv_mean_myocardial_wall_thickness_aha_6_f24129_2_0",
    "lv_mean_myocardial_wall_thickness_aha_7_f24130_2_0",
    "lv_mean_myocardial_wall_thickness_aha_8_f24131_2_0",
    "lv_mean_myocardial_wall_thickness_aha_9_f24132_2_0",
    "lv_mean_myocardial_wall_thickness_aha_10_f24133_2_0",
    "lv_mean_myocardial_wall_thickness_aha_11_f24134_2_0",
    "lv_mean_myocardial_wall_thickness_aha_12_f24135_2_0",
    "lv_mean_myocardial_wall_thickness_aha_13_f24136_2_0",
    "lv_mean_myocardial_wall_thickness_aha_14_f24137_2_0",
    "lv_mean_myocardial_wall_thickness_aha_15_f24138_2_0",
    "lv_mean_myocardial_wall_thickness_aha_16_f24139_2_0",
    "lv_mean_myocardial_wall_thickness_global_f24140_2_0",
    "lv_circumferential_strain_aha_1_f24141_2_0",
    "lv_circumferential_strain_aha_2_f24142_2_0",
    "lv_circumferential_strain_aha_3_f24143_2_0",
    "lv_circumferential_strain_aha_4_f24144_2_0",
    "lv_circumferential_strain_aha_5_f24145_2_0",
    "lv_circumferential_strain_aha_6_f24146_2_0",
    "lv_circumferential_strain_aha_7_f24147_2_0",
    "lv_circumferential_strain_aha_8_f24148_2_0",
    "lv_circumferential_strain_aha_9_f24149_2_0",
    "lv_circumferential_strain_aha_10_f24150_2_0",
    "lv_circumferential_strain_aha_11_f24151_2_0",
    "lv_circumferential_strain_aha_12_f24152_2_0",
    "lv_circumferential_strain_aha_13_f24153_2_0",
    "lv_circumferential_strain_aha_14_f24154_2_0",
    "lv_circumferential_strain_aha_15_f24155_2_0",
    "lv_circumferential_strain_aha_16_f24156_2_0",
    "lv_circumferential_strain_global_f24157_2_0",
    "lv_radial_strain_aha_1_f24158_2_0",
    "lv_radial_strain_aha_2_f24159_2_0",
    "lv_radial_strain_aha_3_f24160_2_0",
    "lv_radial_strain_aha_4_f24161_2_0",
    "lv_radial_strain_aha_5_f24162_2_0",
    "lv_radial_strain_aha_6_f24163_2_0",
    "lv_radial_strain_aha_7_f24164_2_0",
    "lv_radial_strain_aha_8_f24165_2_0",
    "lv_radial_strain_aha_9_f24166_2_0",
    "lv_radial_strain_aha_10_f24167_2_0",
    "lv_radial_strain_aha_11_f24168_2_0",
    "lv_radial_strain_aha_12_f24169_2_0",
    "lv_radial_strain_aha_13_f24170_2_0",
    "lv_radial_strain_aha_14_f24171_2_0",
    "lv_radial_strain_aha_15_f24172_2_0",
    "lv_radial_strain_aha_16_f24173_2_0",
    "lv_radial_strain_global_f24174_2_0",
    "lv_longitudinal_strain_segment_1_f24175_2_0",
    "lv_longitudinal_strain_segment_2_f24176_2_0",
    "lv_longitudinal_strain_segment_3_f24177_2_0",
    "lv_longitudinal_strain_segment_4_f24178_2_0",
    "lv_longitudinal_strain_segment_5_f24179_2_0",
    "lv_longitudinal_strain_segment_6_f24180_2_0",
    "lv_longitudinal_strain_global_f24181_2_0"
]

df_heart = pd.DataFrame(columns=['BAG', 'IDP', 'Log_Odds', 'SE', 'Z', 'P_value', 'OR', 'OR_CI_Lower', 'OR_CI_Upper'])
i = 0
for idp in idp_list_heart:
    tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/BWAS/Heart', 'NSS_logistic_results_' + f'{idp}' + '.tsv')
    i += 1
    if os.path.exists(tsv):
        print(f'Collect: {tsv}')
        df = pd.read_csv(tsv, sep='\t')
        
        if i == 1:
            df_heart = df
        else:
            df_heart = pd.concat([df_heart, df], ignore_index=True)
    else:
        print(f'Model does not converge: {tsv}')
df_heart['Organ'] = 'Heart'
df_heart.to_csv('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_Heart.tsv', index=False, sep='\t')

########################################################################################################################
### 88 Eye
########################################################################################################################
idp_list_eye = [
    "overall_macular_thickness_left_f27800_0_0",
    "macular_thickness_at_the_central_subfield_left_f27802_0_0",
    "macular_thickness_at_the_inner_inferior_subfield_left_f27804_0_0",
    "macular_thickness_at_the_inner_nasal_subfield_left_f27806_0_0",
    "macular_thickness_at_the_inner_superior_subfield_left_f27808_0_0",
    "macular_thickness_at_the_inner_temporal_subfield_left_f27810_0_0",
    "macular_thickness_at_the_outer_inferior_subfield_left_f27812_0_0",
    "macular_thickness_at_the_outer_nasal_subfield_left_f27814_0_0",
    "macular_thickness_at_the_outer_superior_subfield_left_f27816_0_0",
    "macular_thickness_at_the_outer_temporal_subfield_left_f27818_0_0",
    "average_retinal_nerve_fibre_layer_thickness_left_f28500_0_0",
    "average_inner_nuclear_layer_thickness_left_f28502_0_0",
    "average_ganglion_cellinner_plexiform_layer_thickness_left_f28504_0_0",
    "inlelm_thickness_of_the_central_subfield_left_f28506_0_0",
    "inlelm_thickness_of_the_inner_subfield_left_f28508_0_0",
    "inlelm_thickness_of_the_outer_subfield_left_f28510_0_0",
    "average_inlelm_thickness_left_f28512_0_0",
    "elmisos_thickness_of_central_subfield_left_f28514_0_0",
    "elmisos_thickness_of_inner_subfield_left_f28516_0_0",
    "elmisos_thickness_of_outer_subfield_left_f28518_0_0",
    "average_elmisos_thickness_left_f28520_0_0",
    "isosrpe_thickness_of_central_subfield_left_f28522_0_0",
    "isosrpe_thickness_of_inner_subfield_left_f28524_0_0",
    "isosrpe_thickness_of_outer_subfield_left_f28526_0_0",
    "average_isosrpe_thickness_left_f28528_0_0",
    "inlrpe_thickness_of_central_subfield_left_f28530_0_0",
    "inlrpe_thickness_of_inner_subfield_left_f28532_0_0",
    "inlrpe_thickness_of_outer_subfield_left_f28534_0_0",
    "average_inlrpe_thickness_left_f28536_0_0",
    "overall_average_retinal_pigment_epithelium_thickness_left_f27822_0_0",
    "retinal_pigment_epithelium_thickness_at_central_subfield_left_f27824_0_0",
    "retinal_pigment_epithelium_thickness_at_inner_inferior_subfield_left_f27826_0_0",
    "retinal_pigment_epithelium_thickness_at_inner_nasal_subfield_left_f27828_0_0",
    "retinal_pigment_epithelium_thickness_at_inner_superior_subfield_left_f27830_0_0",
    "retinal_pigment_epithelium_thickness_at_inner_temporal_subfield_left_f27832_0_0",
    "retinal_pigment_epithelium_thickness_at_outer_inferior_subfield_left_f27834_0_0",
    "retinal_pigment_epithelium_thickness_at_outer_nasal_subfield_left_f27836_0_0",
    "retinal_pigment_epithelium_thickness_at_outer_superior_subfield_left_f27838_0_0",
    "retinal_pigment_epithelium_thickness_at_outer_temporal_subfield_left_f27840_0_0",
    "disc_diameter_after_inverse_rank_normal_transformation_left_f27851_0_0",
    "mean_of_vertical_disc_diameter_left_f27853_0_0",
    "vertical_cup_to_disc_ratio_vcdr_regressed_and_transformed_left_f27855_0_0",
    "vertical_cup_to_disc_ratio_vcdr_left_f27857_0_0",
    "total_macular_volume_left_f27820_0_0",
    "overall_macular_thickness_right_f27801_0_0",
    "macular_thickness_at_the_central_subfield_right_f27803_0_0",
    "macular_thickness_at_the_inner_inferior_subfield_right_f27805_0_0",
    "macular_thickness_at_the_inner_nasal_subfield_right_f27807_0_0",
    "macular_thickness_at_the_inner_superior_subfield_right_f27809_0_0",
    "macular_thickness_at_the_inner_temporal_subfield_right_f27811_0_0",
    "macular_thickness_at_the_outer_inferior_subfield_right_f27813_0_0",
    "macular_thickness_at_the_outer_nasal_subfield_right_f27815_0_0",
    "macular_thickness_at_the_outer_superior_subfield_right_f27817_0_0",
    "macular_thickness_at_the_outer_temporal_subfield_right_f27819_0_0",
    "average_retinal_nerve_fibre_layer_thickness_right_f28501_0_0",
    "average_inner_nuclear_layer_thickness_right_f28503_0_0",
    "average_ganglion_cellinner_plexiform_layer_thickness_right_f28505_0_0",
    "inlelm_thickness_of_the_central_subfield_right_f28507_0_0",
    "inlelm_thickness_of_the_inner_subfield_right_f28509_0_0",
    "inlelm_thickness_of_the_outer_subfield_right_f28511_0_0",
    "average_inlelm_thickness_right_f28513_0_0",
    "elmisos_thickness_of_central_subfield_right_f28515_0_0",
    "elmisos_thickness_of_inner_subfield_right_f28517_0_0",
    "elmisos_thickness_of_outer_subfield_right_f28519_0_0",
    "average_elmisos_thickness_right_f28521_0_0",
    "isosrpe_thickness_of_central_subfield_right_f28523_0_0",
    "isosrpe_thickness_of_inner_subfield_right_f28525_0_0",
    "isosrpe_thickness_of_outer_subfield_right_f28527_0_0",
    "average_isosrpe_thickness_right_f28529_0_0",
    "inlrpe_thickness_of_central_subfield_right_f28531_0_0",
    "inlrpe_thickness_of_inner_subfield_right_f28533_0_0",
    "inlrpe_thickness_of_outer_subfield_right_f28535_0_0",
    "average_inlrpe_thickness_right_f28537_0_0",
    "overall_average_retinal_pigment_epithelium_thickness_right_f27823_0_0",
    "retinal_pigment_epithelium_thickness_at_central_subfield_right_f27825_0_0",
    "retinal_pigment_epithelium_thickness_at_inner_inferior_subfield_right_f27827_0_0",
    "retinal_pigment_epithelium_thickness_at_inner_nasal_subfield_right_f27829_0_0",
    "retinal_pigment_epithelium_thickness_at_inner_superior_subfield_right_f27831_0_0",
    "retinal_pigment_epithelium_thickness_at_inner_temporal_subfield_right_f27833_0_0",
    "retinal_pigment_epithelium_thickness_at_outer_inferior_subfield_right_f27835_0_0",
    "retinal_pigment_epithelium_thickness_at_outer_nasal_subfield_right_f27837_0_0",
    "retinal_pigment_epithelium_thickness_at_outer_superior_subfield_right_f27839_0_0",
    "retinal_pigment_epithelium_thickness_at_outer_temporal_subfield_right_f27841_0_0",
    "disc_diameter_after_inverse_rank_normal_transformation_right_f27852_0_0",
    "mean_of_vertical_disc_diameter_right_f27854_0_0",
    "vertical_cup_to_disc_ratio_vcdr_regressed_and_transformed_right_f27856_0_0",
    "vertical_cup_to_disc_ratio_vcdr_right_f27858_0_0",
    "total_macular_volume_right_f27821_0_0"
]
df_eye = pd.DataFrame(columns=['BAG', 'IDP', 'Log_Odds', 'SE', 'Z', 'P_value', 'OR', 'OR_CI_Lower', 'OR_CI_Upper'])
i = 0
for idp in idp_list_eye:
    tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/BWAS/Eye', 'NSS_logistic_results_' + f'{idp}' + '.tsv')
    i += 1
    if os.path.exists(tsv):
        print(f'Collect: {tsv}')
        df = pd.read_csv(tsv, sep='\t')
        
        if i == 1:
            df_eye = df
        else:
            df_eye = pd.concat([df_eye, df], ignore_index=True)
    else:
        print(f'Model does not converge: {tsv}')
df_eye['Organ'] = 'Eye'
df_eye.to_csv('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_Eye.tsv', index=False, sep='\t')

########################################################################################################################
### 29 Abdominal
########################################################################################################################
idp_list_abdominal = [
    "Visceral_fat_volume_21085_2_0",
    "Pancreas_PDFF_fat_fraction_21090_2_0",
    "Anterior_thigh_fat_free_muscle_volume_right_22403_2_0",
    "Posterior_thigh_fat_free_muscle_volume_right_22404_2_0",
    "Anterior_thigh_fat_free_muscle_volume_left_22405_2_0",
    "Posterior_thigh_fat_free_muscle_volume_left_22406_2_0",
    "Abdominal_subcutaneous_adipose_tissue_volume_ASAT_22408_2_0",
    "Total_trunk_fat_volume_22410_2_0",
    "Total_abdominal_adipose_tissue_index_22432_2_0",
    "Abdominal_fat_ratio_22434_2_0",
    "Muscle_fat_infiltration_22435_2_0",
    "Posterior_thigh_muscle_fat_infiltration_MFI_left_23355_2_0",
    "Posterior_thigh_muscle_fat_infiltration_MFI_right_23356_2_0",
    "Anterior_thigh_muscle_fat_infiltration_MFI_left_24353_2_0",
    "Anterior_thigh_muscle_fat_infiltration_MFI_right_24354_2_0",
    "Proton_density_fat_fraction_PDFF_40061_2_0",
    "Left_kidney_volume_21081_2_0",
    "Kidney_parenchyma_right_21162_2_0",
    "Kidney_distance_21163_2_0",
    "Liver_volume_21080_2_0",
    "Liver_PDFF_fat_fraction_21088_2_0",
    "Liver_iron_21089_2_0",
    "Liver_iron_corrected_T1_ct1_40062_2_0",
    "Pancreas_volume_21087_2_0",
    "Pancreas_PDFF_fat_fraction_21090_2_0_pancreas",
    "Pancreas_iron_21091_2_0",
    "Spleen_volume_21083_2_0",
    "Spleen_iron_IDEAL_21170_2_0",
    "Spleen_iron_protocol_normalised_21173_2_0",
]
df_abdominal = pd.DataFrame(columns=['BAG', 'IDP', 'Log_Odds', 'SE', 'Z', 'P_value', 'OR', 'OR_CI_Lower', 'OR_CI_Upper'])
i = 0
for idp in idp_list_abdominal:
    tsv = os.path.join('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/BWAS/Abdominal', 'NSS_logistic_results_' + f'{idp}' + '.tsv')
    i += 1
    if os.path.exists(tsv):
        print(f'Collect: {tsv}')
        df = pd.read_csv(tsv, sep='\t')
        
        if i == 1:
            df_abdominal = df
        else:
            df_abdominal = pd.concat([df_abdominal, df], ignore_index=True)
    else:
        print(f'Model does not converge: {tsv}')
idp_dict_abdominal = {
    "Visceral_fat_volume_21085_2_0": "Adipose",
    "Pancreas_PDFF_fat_fraction_21090_2_0": "Pancreas",
    "Anterior_thigh_fat_free_muscle_volume_right_22403_2_0": "Adipose",
    "Posterior_thigh_fat_free_muscle_volume_right_22404_2_0": "Adipose",
    "Posterior_thigh_fat_free_muscle_volume_left_22406_2_0": "Adipose",
    "Abdominal_subcutaneous_adipose_tissue_volume_ASAT_22408_2_0": "Adipose",
    "Total_trunk_fat_volume_22410_2_0": "Adipose",
    "Total_abdominal_adipose_tissue_index_22432_2_0": "Adipose",
    "Abdominal_fat_ratio_22434_2_0": "Adipose",
    "Muscle_fat_infiltration_22435_2_0": "Adipose",
    "Posterior_thigh_muscle_fat_infiltration_MFI_left_23355_2_0": "Adipose",
    "Posterior_thigh_muscle_fat_infiltration_MFI_right_23356_2_0": "Adipose",
    "Anterior_thigh_muscle_fat_infiltration_MFI_left_24353_2_0": "Adipose",
    "Anterior_thigh_muscle_fat_infiltration_MFI_right_24354_2_0": "Adipose",
    "Proton_density_fat_fraction_PDFF_40061_2_0": "Adipose",
    "Left_kidney_volume_21081_2_0": "Kidney",
    "Kidney_parenchyma_right_21162_2_0": "Kidney",
    "Kidney_distance_21163_2_0": "Kidney",
    "Liver_volume_21080_2_0": "Liver",
    "Liver_PDFF_fat_fraction_21088_2_0": "Adipose",
    "Liver_iron_21089_2_0": "Liver",
    "Liver_iron_corrected_T1_ct1_40062_2_0": "Liver",
    "Pancreas_volume_21087_2_0": "Pancreas",
    "Pancreas_PDFF_fat_fraction_21090_2_0_pancreas": "Pancreas",
    "Pancreas_iron_21091_2_0": "Pancreas",
    "Spleen_volume_21083_2_0": "Spleen",
    "Spleen_iron_IDEAL_21170_2_0": "Spleen",
    "Spleen_iron_protocol_normalised_21173_2_0": "Spleen"
}
# Add Organ column to df_abdominal
df_abdominal["Organ"] = df_abdominal["IDP"].map(idp_dict_abdominal)
df_abdominal.to_csv('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_Abdominal.tsv', index=False, sep='\t')

#### concate the results across the different organs
df_final = pd.concat([df_brain_gm, df_brain_wm, df_brain_fc,
                      df_heart,
                      df_eye,
                      df_abdominal
                      ], ignore_index=True)

### remove OR>10 becasue the signals are not true
df_final = df_final[df_final["OR"] <= 10]
df_final = df_final[~df_final["BAG"].str.endswith("_MRIBAG", na=False)].copy()
df_final.to_csv('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result.tsv', index=False, sep='\t')
df_final_sig = df_final[df_final['P_value'] < 0.05/len(df_final['IDP'].unique())]
df_final_sig.to_csv('/Users/hao/cubic-home/Reproducibile_paper/SleepAging/NSS/Result/BWAS_result_bon.tsv', index=False, sep='\t')

print('stop here...')