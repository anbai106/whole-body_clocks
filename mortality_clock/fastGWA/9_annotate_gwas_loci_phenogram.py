import numpy as np
import os
import pandas as pd
from scipy import stats
def annotate_genomic_loci(output_dir_metbag, output_dir_result):
    """
    https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ to annotate the cytogenetic band regions
    :param output_dir:
    :param output_dir_result:
    :param Zhao_loci_results:
    :return:
    """

    ### read the cytogenetic band
    df_cyto = pd.read_csv('/Users/hao/cubic-home/Dataset/GRch37_cytoband/cytoBand.txt.gz', sep='\t', header=None)

    ### brain
    df_final_metbag = pd.DataFrame(columns=['GenomicLocus', 'TopLeadSNP', 'Chromosome', 'Position', 'P-value', 'MappedGene', 'Phenotype', 'Cytogenetic_region_fuma', 'Cytogenetic_region_UCSC', 'Novel', 'Organ', 'colorgroup'])
    for pheno in ['Endocrine', 'Digestive', 'Hepatic', 'Immune', "Metabolic"]:
        df_loci = pd.read_csv(os.path.join(output_dir_metbag, pheno, 'GenomicRiskLoci.txt'), sep='\t')[['GenomicLocus', 'rsID', 'chr', 'pos', 'p']]
        df_loci = df_loci.loc[(df_loci['p'] < (5e-8)/5)]
        df_loci['pos'] = df_loci['pos'].astype(str)
        df_gene = pd.read_csv(os.path.join(output_dir_metbag, pheno, 'genes.txt'), sep='\t')[['symbol', 'GenomicLocus']]
        df_gene = df_gene.groupby('GenomicLocus', as_index=False).agg(list)
        df_gene['GenomicLocus'] = df_gene['GenomicLocus'].astype(str)
        df_gene['GenomicLocus'] = df_gene['GenomicLocus'].str.split(':').str[0].tolist()
        df_gene['GenomicLocus'] = df_gene['GenomicLocus'].astype(int)
        df_loci = df_loci.merge(df_gene, how='left', left_on='GenomicLocus', right_on='GenomicLocus')
        df_loci['Phenotype'] = pheno
        df_loci.rename({'rsID': 'TopLeadSNP'}, axis=1, inplace=True)
        df_loci.rename({'chr': 'Chromosome'}, axis=1, inplace=True)
        df_loci.rename({'pos': 'Position'}, axis=1, inplace=True)
        df_loci.rename({'p': 'P-value'}, axis=1, inplace=True)
        df_loci.rename({'symbol': 'MappedGene'}, axis=1, inplace=True)
        topleadsnp_list = list(df_loci['TopLeadSNP'])
        novel_list = []
        pos_list = list(df_loci['Position'])
        region_list = []
        region_ucsc_list = []
        chr_list = list(df_loci['Chromosome'])

        ### read GWAS Catalog file to define the cytogenetic region and the
        df_gwas_catalog = pd.read_csv(os.path.join(output_dir_metbag, pheno, 'gwascatalog.txt'), sep='\t')
        for i in range(len(topleadsnp_list)):
            if int(pos_list[i]) in df_gwas_catalog['bp'].to_list():
                novel = "N"
                novel_list.append(novel)
            else:
                novel= "Y"
                novel_list.append(novel)
            ### find rows whose "bp" value is closest to Position value
            df_gwas_catalog_chr = df_gwas_catalog.loc[df_gwas_catalog['chr'].isin([chr_list[i]])]
            df_closest = df_gwas_catalog_chr.iloc[(df_gwas_catalog_chr['bp'] - int(pos_list[i])).abs().argsort()[:1]]
            if df_closest.shape[0] != 1:
                region = 'https://www.ncbi.nlm.nih.gov/genome/gdv/ or https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/'
            else:
                region = df_closest['Region'].values[0]
            region_list.append(region)
            chr = chr_list[i]
            pos = pos_list[i]
            pos_row = df_cyto.loc[(df_cyto[1] <= int(pos)) & (df_cyto[2] >= int(pos))]
            pos_row = pos_row.loc[pos_row[0].isin(['chr' + str(chr)])]
            if pos_row.shape[0] == 1:
                region_ucsc = str(chr) + pos_row[3].values[0]
                region_ucsc_list.append(region_ucsc)
            else:
                raise Exception('Sth is wrong here...')
            if region_ucsc != region:
                print("For top lead SNP: %s" % topleadsnp_list[i])
                print("Cytogenetic region from UCSC is not the same from what I defined from FUMA via GWAS Catalog...")

        df_loci['Cytogenetic_region_fuma'] = region_list
        df_loci['Cytogenetic_region_UCSC'] = region_ucsc_list
        df_loci['Novel'] = novel_list
        df_loci['Organ'] = 'PhenoBAG'
        if pheno == 'Endocrine':
            df_loci['colorgroup'] = 'Endocrine'
        elif pheno == 'Digestive':
            df_loci['colorgroup'] = 'Digestive'
        elif pheno == 'Hepatic':
            df_loci['colorgroup'] = 'Hepatic'
        elif pheno == 'Immune':
            df_loci['colorgroup'] = 'Immune'
        elif pheno == 'Metabolic':
            df_loci['colorgroup'] = 'Metabolic'
        else:
            raise Exception('Stop here... sth is wrong...')

        df_final_metbag = df_final_metbag.append(df_loci, ignore_index=True)

    ### Prepare data for PhenoGram
    df_phenogram = df_final_metbag[['TopLeadSNP', 'Chromosome', 'Position', 'Phenotype', 'Cytogenetic_region_UCSC', 'Organ', 'colorgroup']]
    df_phenogram.rename({'TopLeadSNP': 'snp'}, axis=1, inplace=True)
    df_phenogram.rename({'Chromosome': 'chr'}, axis=1, inplace=True)
    df_phenogram.rename({'Position': 'pos'}, axis=1, inplace=True)
    df_phenogram.rename({'Phenotype': 'phenotype'}, axis=1, inplace=True)
    df_phenogram.rename({'Cytogenetic_region_UCSC': 'annotation'}, axis=1, inplace=True)
    df_phenogram.rename({'Organ': 'ethnicity'}, axis=1, inplace=True)

    df_final_metbag.to_csv(os.path.join(output_dir_result, 'Fuma_loci_annotation_fastGWA.tsv'), index=False, sep='\t', encoding='utf-8')
    df_phenogram.to_csv(os.path.join(output_dir_result, 'PhenoGram_loci_annotation_fastGWA.tsv'), index=False, sep='\t', encoding='utf-8')

output_dir_metbag = '/Users/hao/cubic-home/Reproducibile_paper/UKBB_metabolomics/fuma'
output_dir_result = '/Users/hao/cubic-home/Reproducibile_paper/UKBB_metabolomics/Result'
annotate_genomic_loci(output_dir_metbag, output_dir_result)