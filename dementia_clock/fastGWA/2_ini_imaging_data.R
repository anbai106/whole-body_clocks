library(stringr)
### remove the R environment each time you rerun the code
rm(list = ls())


lm_imaging_data <- function(data_src, cov_src, output_dir_final){
  imaging_data = as.data.frame(data.table::fread(data_src, header=T,stringsAsFactors=FALSE))
  covar_data = as.data.frame(data.table::fread(cov_src, header=T,stringsAsFactors=FALSE))

  imaging.norm = imaging.resid = imaging.resid.norm = imaging_data
  for(i in 3:3){
    imaging.norm[,i] = qnorm((rank(imaging_data[,i],na.last="keep")-0.5)/sum(!is.na(imaging_data[,i])));print(i) ### see post here: https://www.biostars.org/p/80597/ and this: https://www.biostars.org/p/312945/
  }

  fin = cbind(imaging.norm,covar_data)
  for(i in 3:3){
    imaging.resid[,i] = residuals(lm(data=fin, fin[,i] ~ sex_f31_0_0+age_when_attended_assessment_centre_f21003_0_0+genetic_principal_components_f22009_0_1+genetic_principal_components_f22009_0_2+genetic_principal_components_f22009_0_3+genetic_principal_components_f22009_0_4+genetic_principal_components_f22009_0_5+genetic_principal_components_f22009_0_6+genetic_principal_components_f22009_0_7+genetic_principal_components_f22009_0_8+genetic_principal_components_f22009_0_9+genetic_principal_components_f22009_0_10+genetic_principal_components_f22009_0_11+genetic_principal_components_f22009_0_12+genetic_principal_components_f22009_0_13+genetic_principal_components_f22009_0_14+genetic_principal_components_f22009_0_15+genetic_principal_components_f22009_0_16+genetic_principal_components_f22009_0_17+genetic_principal_components_f22009_0_18+genetic_principal_components_f22009_0_19+genetic_principal_components_f22009_0_20+genetic_principal_components_f22009_0_21+genetic_principal_components_f22009_0_22+genetic_principal_components_f22009_0_23+genetic_principal_components_f22009_0_24+genetic_principal_components_f22009_0_25+genetic_principal_components_f22009_0_26+genetic_principal_components_f22009_0_27+genetic_principal_components_f22009_0_28+genetic_principal_components_f22009_0_29+genetic_principal_components_f22009_0_30+genetic_principal_components_f22009_0_31+genetic_principal_components_f22009_0_32+genetic_principal_components_f22009_0_33+genetic_principal_components_f22009_0_34+genetic_principal_components_f22009_0_35+genetic_principal_components_f22009_0_36+genetic_principal_components_f22009_0_37+genetic_principal_components_f22009_0_38+genetic_principal_components_f22009_0_39+genetic_principal_components_f22009_0_40+weight_f21002_0_0+standing_height_f50_0_0+waist_circumference_f48_0_0+body_mass_index_bmi_f23104_0_0+diastolic_blood_pressure_automated_reading_f4079_0_0+systolic_blood_pressure_automated_reading_f4080_0_0 + I(age_when_attended_assessment_centre_f21003_0_0^2) +
                                       age_when_attended_assessment_centre_f21003_0_0*sex_f31_0_0 + I(age_when_attended_assessment_centre_f21003_0_0^2)*sex_f31_0_0, na.action=na.exclude))
    imaging.resid.norm[,i] = qnorm((rank(imaging.resid[,i],na.last="keep")-0.5)/sum(!is.na(imaging.resid[,i])))
    outlier = mean(imaging_data[,i],na.rm=T)+6*sd(imaging_data[,i],na.rm=T)
    tmp=which(imaging_data[,i]>outlier) ########### Remove outliers
    if(length(tmp)>0) imaging.resid.norm[tmp,i] = NA
    print(i)
  }

  write.table(imaging.resid.norm,file=paste(output_dir_final, "/EPOCH_pheno_normalized_residualized_with_related_indi.phen", sep=""),sep = " ",quote = F,col.names = F,row.names = F)
}

### 4 MetBAG
organ_list = c('Endocrine', 'Digestive', 'Hepatic', 'Metabolic')
output_dir = '/Users/hao/cubic-home/Reproducibile_paper/WholeBodyClock/'
for (organ in organ_list) {
  output_dir_final <- paste(output_dir, paste0(organ, "_metabolomics_dementia_clock"), 'fastGWA', 'data', sep = "/")
  data_src = paste(output_dir_final, 'EPOCH_pheno.txt', sep='/')
  cov_src = paste(output_dir_final, 'EPOCH_cov.txt', sep='/')
  lm_imaging_data(data_src, cov_src, output_dir_final)
}

### 7 ProtBAG
organ_list = c('Brain', 'Hepatic', 'Endocrine', 'Heart', 'Immune', 'Reproductive_female', 'Reproductive_male')
for (organ in organ_list) {
  output_dir_final <- paste(output_dir, paste0(organ, "_proteomics_dementia_clock"), 'fastGWA', 'data', sep = "/")
  data_src = paste(output_dir_final, 'EPOCH_pheno.txt', sep='/')
  cov_src = paste(output_dir_final, 'EPOCH_cov.txt', sep='/')
  lm_imaging_data(data_src, cov_src, output_dir_final)
}
