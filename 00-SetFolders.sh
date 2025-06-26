#!/bin/bash

#### Tools needed
# bcftools https://www.htslib.org/download/
# R
# GCTA-1.94.1 

#### Set folders
mkdir -p ${METAresPath}
mkdir -p ${METAresPath}/logs
mkdir -p ${METAresPath}/SampleScheme
mkdir -p ${METAresPath}/StandErrScheme
mkdir -p ${METAresPath}/MergedMeta

mkdir -p ${METAresPath}/FilteredMeta
mkdir -p ${METAresPath}/FilteredMeta/plots

mkdir -p ${COJOresPath}
mkdir -p ${COJOresPath}/logs
mkdir -p ${COJOresPath}/Cojo_${phenotype}

mkdir -p ${MstatresPath}
mkdir -p ${MstatresPath}/logs
mkdir -p ${MstatresPath}/plots
mkdir -p ${MstatresPath}/Mstat_${phenotype}

mkdir -p ${METAregrePath}
mkdir -p ${METAregrePath}/logs
mkdir -p ${METAregrePath}/plots
mkdir -p ${METAregrePath}/METAreg_${phenotype}

#### Download Reference files if needed
# 2024 version of dbSNP
# if [ -f ${HomePath}/Resources/GRCh37_latest_dbSNP_all.vcf.gz ]; then
#     echo "GRCh37_latest_dbSNP_all.vcf.gz already exists."
# else
#     echo "GRCh37_latest_dbSNP_all.vcf.gz does not exist. Downloading and processing..."
#     cd ${HomePath}/Resources || exit
#     wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_dbSNP_all.vcf.gz

#     # 38
#     zcat GRCh37_latest_dbSNP_all.vcf.gz | awk 'NR==1 || /^NC/ {print $1, $2, $3, $4, $5}' | gzip > GRCh37_latest_dbSNP_all.rsid.gz
#     zcat GRCh37_latest_dbSNP_all.rsid.gz | awk 'NR>1 { print $1 }' | uniq > original_chr_name.txt
#     awk '{ print $0, NR }' original_chr_name.txt > chr_name_convert.txt

#     bcftools annotate --rename-chrs chr_name_convert.txt GRCh37_latest_dbSNP_all.rsid.gz | gzip > GRCh37_latest_dbSNP_all.rsid.rename.gz
#     bcftools annotate --rename-chrs chr_name_convert.txt GRCh37_latest_dbSNP_all.vcf.gz | gzip > GRCh37_latest_dbSNP_all.rsid.rename.gz
    
#     # to do 
# fi

#### Resources files
if [ -f ${HomePath}/Resources/1000Geur/1000G_hg19_eur.bim ]; then
    echo "1000G eur reference files are exist."
else
    echo "1000G eur reference files does not exist. Please check the file path."
    exit 1
fi

if [ -f ${HomePath}/Resources/cohort_lambda.tsv ]; then
    echo "Cohort lambda file exists."
    cohortLambda=${HomePath}/Resources/cohort_lambda.tsv
else
    echo "Cohort lambda file does not exist. Please check the file path."
    exit 1
fi

if [ -f ${HomePath}/Resources/cohort_age.tsv ]; then
    echo "Cohort age file exists."
else
    echo "Cohort age file does not exist. Please check the file path."
    exit 1
fi