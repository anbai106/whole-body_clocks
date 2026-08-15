#!/bin/bash
#SBATCH --partition=all
#SBATCH --job-name=BWAS
#SBATCH --array=0-119
#SBATCH --mem-per-cpu=12G
#SBATCH --time=0-02:59:00
#SBATCH --output=/cbica/home/wenju/output/BWAS_%A_%a.out
#SBATCH --error=/cbica/home/wenju/output/BWAS_%A_%a.err


idp_array=( MUSE_Volume_23 MUSE_Volume_30 MUSE_Volume_31 MUSE_Volume_32 MUSE_Volume_36 MUSE_Volume_37 MUSE_Volume_38 MUSE_Volume_39 MUSE_Volume_47 MUSE_Volume_48 MUSE_Volume_55 MUSE_Volume_56 MUSE_Volume_57 MUSE_Volume_58 MUSE_Volume_59 MUSE_Volume_60 MUSE_Volume_71 MUSE_Volume_72 MUSE_Volume_73 MUSE_Volume_75 MUSE_Volume_76 MUSE_Volume_100 MUSE_Volume_101 MUSE_Volume_102 MUSE_Volume_103 MUSE_Volume_104 MUSE_Volume_105 MUSE_Volume_106 MUSE_Volume_107 MUSE_Volume_108 MUSE_Volume_109 MUSE_Volume_112 MUSE_Volume_113 MUSE_Volume_114 MUSE_Volume_115 MUSE_Volume_116 MUSE_Volume_117 MUSE_Volume_118 MUSE_Volume_119 MUSE_Volume_120 MUSE_Volume_121 MUSE_Volume_122 MUSE_Volume_123 MUSE_Volume_124 MUSE_Volume_125 MUSE_Volume_128 MUSE_Volume_129 MUSE_Volume_132 MUSE_Volume_133 MUSE_Volume_134 MUSE_Volume_135 MUSE_Volume_136 MUSE_Volume_137 MUSE_Volume_138 MUSE_Volume_139 MUSE_Volume_140 MUSE_Volume_141 MUSE_Volume_142 MUSE_Volume_143 MUSE_Volume_144 MUSE_Volume_145 MUSE_Volume_146 MUSE_Volume_147 MUSE_Volume_148 MUSE_Volume_149 MUSE_Volume_150 MUSE_Volume_151 MUSE_Volume_152 MUSE_Volume_153 MUSE_Volume_154 MUSE_Volume_155 MUSE_Volume_156 MUSE_Volume_157 MUSE_Volume_160 MUSE_Volume_161 MUSE_Volume_162 MUSE_Volume_163 MUSE_Volume_164 MUSE_Volume_165 MUSE_Volume_166 MUSE_Volume_167 MUSE_Volume_168 MUSE_Volume_169 MUSE_Volume_170 MUSE_Volume_171 MUSE_Volume_172 MUSE_Volume_173 MUSE_Volume_174 MUSE_Volume_175 MUSE_Volume_176 MUSE_Volume_177 MUSE_Volume_178 MUSE_Volume_179 MUSE_Volume_180 MUSE_Volume_181 MUSE_Volume_182 MUSE_Volume_183 MUSE_Volume_184 MUSE_Volume_185 MUSE_Volume_186 MUSE_Volume_187 MUSE_Volume_190 MUSE_Volume_191 MUSE_Volume_192 MUSE_Volume_193 MUSE_Volume_194 MUSE_Volume_195 MUSE_Volume_196 MUSE_Volume_197 MUSE_Volume_198 MUSE_Volume_199 MUSE_Volume_200 MUSE_Volume_201 MUSE_Volume_202 MUSE_Volume_203 MUSE_Volume_204 MUSE_Volume_205 MUSE_Volume_206 MUSE_Volume_207 )
idp=${idp_array[$SLURM_ARRAY_TASK_ID]}

sleep_dir=/cbica/home/wenju/Reproducibile_paper/SleepAging/NSS/data
idp_tsv=/cbica/home/wenju/Reproducibile_paper/BrainAge/data/imaging/T1_MUSE_GM.tsv
cov_tsv=/cbica/home/wenju/Reproducibile_paper/PRS_UKBB/prediction/data/UKBB_fullsample_covariate.csv
output_dir=/cbica/home/wenju/Reproducibile_paper/SleepAging/NSS/BWAS/Brain/T1
mkdir -p $output_dir

output_file="${output_dir}/NSS_logistic_results_${idp}.tsv"
if [ ! -f ${output_file} ]; then
  echo "Run GAMLSS for: ${idp}..."
  bash /cbica/home/wenju/Project/SleepAging/NSS/BWAS/Brain/T1/fit_logistic.sh ${sleep_dir} ${idp_tsv} ${output_dir} ${idp} ${cov_tsv}
else
  :
fi
