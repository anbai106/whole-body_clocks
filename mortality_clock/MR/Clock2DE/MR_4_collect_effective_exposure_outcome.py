import pandas as pd
import os

output_dir = '/cbica/home/wenju/Reproducibile_paper/AbdoImaging/MR/Clock2DE'
output_dir_results = "/cbica/home/wenju/Reproducibile_paper/AbdoImaging/Result"
list_organ = ['brain', "adipose", "heart", "kidney", "liver", "pancreas", "spleen"]
n = 0
finngen_pheno = pd.read_csv('/cbica/projects/MULTI/processed/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv', sep='\t')['phenocode'].to_list()

for organ in list_organ:
    ## PGC
    for pgc in ['AD', 'ADHD', 'BIP', 'SCZ']:
        tsv = os.path.join(output_dir, organ, "MR_" + organ + "_2_" + pgc + "_OR.tsv")
        if os.path.exists(tsv):
            n += 1
            if n == 1:
                df_final = pd.read_csv(tsv, delimiter='\t')
                del df_final['id.exposure']
                del df_final['id.outcome']
            else:
                df = pd.read_csv(tsv, delimiter='\t')
                del df['id.exposure']
                del df['id.outcome']
                df_final = df_final.append(df, ignore_index=True)
    ## FinnGen
    for finngen in finngen_pheno:
        tsv = os.path.join(output_dir, organ, "MR_" + organ + "_2_" + finngen + "_OR.tsv")
        if os.path.exists(tsv):
            n += 1
            if n == 1:
                df_final = pd.read_csv(tsv, delimiter='\t')
                del df_final['id.exposure']
                del df_final['id.outcome']
            else:
                df = pd.read_csv(tsv, delimiter='\t')
                del df['id.exposure']
                del df['id.outcome']
                df_final = df_final.append(df, ignore_index=True)

n_disease = len(df_final.outcome.unique())
n_organ = len(df_final.exposure.unique())
print("Unique MRIBAG tested: %d" % n_organ)
print("Unique organ disease tested: %d" % n_disease)
df_final['P_bon_n_disease'] = 0.05/n_disease
df_final['P_bon_n_organ'] = 0.05/n_organ
df_final['P_bon_n_MRIBAG_de'] = 0.05/n_organ/n_disease
df_final_ivw = df_final.loc[df_final['method'].isin(['Inverse variance weighted'])]
df_final_de_sig = df_final.loc[(df_final['pval'] <= df_final['P_bon_n_disease'])]
df_final_organ_sig = df_final.loc[(df_final['pval'] <= df_final['P_bon_n_organ'])]
df_final_pheno_de_sig = df_final.loc[(df_final['pval'] <= df_final['P_bon_n_MRIBAG_de'])]
df_final.to_csv(os.path.join(output_dir_results, '2SampleMR_MRIClock2DE_all.tsv'), index=False, sep='\t', encoding='utf-8')
df_final_ivw.to_csv(os.path.join(output_dir_results, '2SampleMR_MRIClock2DE_ivw.tsv'), index=False, sep='\t', encoding='utf-8')
df_final_de_sig.to_csv(os.path.join(output_dir_results, '2SampleMR_MRIClock2DE_all_sig_by_DE.tsv'), index=False, sep='\t', encoding='utf-8')
df_final_organ_sig.to_csv(os.path.join(output_dir_results, '2SampleMR_MRIClock2DE_all_sig_by_MRIBAG.tsv'), index=False, sep='\t', encoding='utf-8')
df_final_pheno_de_sig.to_csv(os.path.join(output_dir_results, '2SampleMR_MRIClock2DE_all_sig_by_MRIBAG_DE.tsv'), index=False, sep='\t', encoding='utf-8')

print("STOP here...")
