# Compress the fastGWA file.
# -j stores only the file name inside the zip, not the full path.
# -9 uses maximum compression.
zip_file=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/Brain_proteomics_mortality_clock/EPOCH_BAG_residual/organ_pheno_normalized_residualized.fastGWA.zip
fastgwa=/cbica/home/wenju/Reproducibile_paper/WholeBodyClock/mortality_clock/fastGWA/output/Brain_proteomics_mortality_clock/EPOCH_BAG_residual/organ_pheno_normalized_residualized.fastGWA
zip -9 -j "${zip_file}" "${fastgwa}"

# Test zip integrity.
zip -T "${zip_file}"

echo "Zip created successfully:"
ls -lh "${zip_file}"

echo "Done."
