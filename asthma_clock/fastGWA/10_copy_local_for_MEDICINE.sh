#!/bin/bash

cd /Users/hao
for organ in Endocrine Digestive Hepatic Immune Metabolic
do
    echo ${organ}
    ### GWAS SSD, Manhattan plot and QQ plot
#    rsync -avz wenju@cubic-login.uphs.upenn.edu:/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/${organ}/organ_pheno_normalized_residualized.fastGWA.zip /Users/hao/${organ}_fastGWA_EUR.zip
    rsync -avz wenju@cubic-login.uphs.upenn.edu:/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/${organ}/QQ_plot.png /Users/hao/Downloads/sex_specific_bag/MetBAG/${organ}_MetBAG_QQ_plot.png
    rsync -avz wenju@cubic-login.uphs.upenn.edu:/cbica/home/wenju/Reproducibile_paper/UKBB_metabolomics/fastGWA_MetBAG/output/${organ}/manhattan_qmplot.png /Users/hao/Downloads/sex_specific_bag/MetBAG/${organ}_MetBAG_manhattan_qmplot.png
done

#for organ in Endocrine Digestive Hepatic Immune Metabolic
#do
#    input_zip="${organ}_fastGWA_EUR.zip"
#    output_txt="${organ}_fastGWA_EUR_noAF1.txt"
#    output_zip="${organ}_fastGWA_EUR_noAF1.zip"
#
#    # Process and remove AF1 column
#    unzip -p "$input_zip" | awk -F'\t' '
#        NR==1 {
#            for (i=1; i<=NF; i++) if ($i=="AF1") col=i;
#        }
#        {
#            for (i=1; i<=NF; i++) if (i!=col) printf "%s%s", $i, (i<NF?OFS:"\n")
#        }' OFS='\t' > "$output_txt"
#
#    # Compress the cleaned file back into a zip archive
#    zip -j "$output_zip" "$output_txt"
#
#    # Remove the temporary text file
#    rm "$output_txt"
#
#    echo "Processed $input_zip -> $output_zip"
#done
