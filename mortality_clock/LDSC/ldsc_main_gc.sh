#!/bin/bash
module load python/anaconda/3

stat_gz1=$1
output_file=$2
ref_ld_chr="/cbica/home/wenju/Project/ldsc/pre-computed_Eu_GWAS/eur_w_ld_chr/"  ### the / at the end is very important, corresponding to the issues here on github: https://github.com/bulik/ldsc/issues/250
w_ld_chr="/cbica/home/wenju/Project/ldsc/pre-computed_Eu_GWAS/eur_w_ld_chr/"

source activate ldsc
echo "Start training"
### softlink to make all smumstats.gz in the same current work_dir

### let's loop all the MAEs for all the organs from our UKBB analysis + FinnGenn + PGC

#### 11 ProtBAG
for bag in Reproductive_female Pulmonary Heart Brain Eye Hepatic Renal Reproductive_male Endocrine Immune Skin
do
  stat_gz2=/cbica/home/wenju/Reproducibile_paper/UKBB_Proteomics/gc/${bag}/${bag}.sumstats.gz
  python /cbica/home/wenju/Project/ldsc/ldsc.py --rg ${stat_gz1},${stat_gz2} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}_ProtBAG_${bag}
done

#### 5 MetBAG
for bag in Endocrine Digestive Hepatic Immune Metabolic
do
  stat_gz2=/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/gc/${bag}/${bag}.sumstats.gz
  python /cbica/home/wenju/Project/ldsc/ldsc.py --rg ${stat_gz1},${stat_gz2} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}_MetBAG_${bag}
done

#### with the 11 pan-disease MAE themselves
### brain
for mae in r1 r2 r3 r4 r5 r6
do
  stat_gz2=/cbica/home/wenju/Reproducibile_paper/Multiorgan_Subtype/output/gc/brain/${mae}/${mae}.sumstats.gz
  python /cbica/home/wenju/Project/ldsc/ldsc.py --rg ${stat_gz1},${stat_gz2} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}_${mae}_brain
done

### eye
for mae in r1 r2 r3
do
  stat_gz2=/cbica/home/wenju/Reproducibile_paper/Multiorgan_Subtype/output/gc/eye/${mae}/${mae}.sumstats.gz
  python /cbica/home/wenju/Project/ldsc/ldsc.py --rg ${stat_gz1},${stat_gz2} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}_${mae}_eye
done

### heart
for mae in r1 r2
do
  stat_gz2=/cbica/home/wenju/Reproducibile_paper/Multiorgan_Subtype/output/gc/heart/${mae}/${mae}.sumstats.gz
  python /cbica/home/wenju/Project/ldsc/ldsc.py --rg ${stat_gz1},${stat_gz2} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}_${mae}_heart
done

### 9DNEs
for mae in AD_SurrealGAN_1 AD_SurrealGAN_2 ASD_1 ASD_2 ASD_3 LLD_1 LLD_2 SCZ_1 SCZ_2
do
  stat_gz2=/cbica/home/wenju/Reproducibile_paper/DNE/output/GWAS/${mae}/gc/${mae}.sumstats.gz
  python /cbica/home/wenju/Project/ldsc/ldsc.py --rg ${stat_gz1},${stat_gz2} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}_${mae}
done

### 3 multi-modal brain BAG
for mae in muse dmri_fullmetric_fa fmri
do
  if [ ${mae} == "dmri_fullmetric_fa" ]; then
    stat_gz2="/cbica/home/wenju/Reproducibile_paper/BrainAge/output/GWAS/${mae}/output/gc/dmri.sumstats.gz"
  else
    stat_gz2="/cbica/home/wenju/Reproducibile_paper/BrainAge/output/GWAS/${mae}/output/gc/${mae}.sumstats.gz"
  fi
  python /cbica/home/wenju/Project/ldsc/ldsc.py --rg ${stat_gz1},${stat_gz2} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}_${mae}
done

### FinnGen
summary_tsv='/cbica/home/wenju/Dataset/FinnGen/GWAS_summary_stats/summary_stats_R9_manifest_5000_cases.tsv'
variable=`awk '{print $1}' ${summary_tsv}`
for de in $variable
do
  if [ ${de} != "phenocode" ]; then
    stat_gz2="/cbica/projects/MULTI/processed/FinnGen/LDSC/finngen_R9_${de}.sumstats.gz"
    python /cbica/home/wenju/Project/ldsc/ldsc.py --rg ${stat_gz1},${stat_gz2} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}_${de}
  fi
done

### PGC
summary_tsv='/cbica/home/wenju/Dataset/PGC/PGC_MUTATE.tsv'
variable=`awk -F'\t' '{print $2}' ${summary_tsv}`
for de in $variable
do
  if [ ${de} != "Phenotype" ]; then
    stat_gz2="/cbica/projects/MULTI/processed/PGC/GWAS_summary_stats/Result/gc/ldsc/${de}.sumstats.gz"
    python /cbica/home/wenju/Project/ldsc/ldsc.py --rg ${stat_gz1},${stat_gz2} --ref-ld-chr ${ref_ld_chr} --w-ld-chr ${w_ld_chr} --out ${output_file}_${de}
  fi
done

echo "Finish!"
conda deactivate