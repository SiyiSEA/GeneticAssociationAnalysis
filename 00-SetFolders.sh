#!/bin/bash

#### Tools
# bcftools https://www.htslib.org/download/
# R
# GCTA-1.94.1 

#### Set folders
mkdir -p ${HomePath}/METALresults
mkdir -p ${HomePath}/METALresults/logs
mkdir -p ${HomePath}/METALresults/SampleScheme
mkdir -p ${HomePath}/METALresults/StandErrScheme
mkdir -p ${HomePath}/METALresults/MergedMeta
mkdir -p ${HomePath}/METALresults/FilteredMeta
mkdir -p ${HomePath}/METALresults/FilteredMeta/plots

mkdir -p ${HomePath}/COJOresults
mkdir -p ${HomePath}/COJOresults/logs

mkdir -p ${HomePath}/Mstatresults/
mkdir -p ${HomePath}/Mstatresults/logs
mkdir -p ${HomePath}/Mstatresults/plots

mkdir -p ${HomePath}/METAregresults/
mkdir -p ${HomePath}/METAregresults/logs
mkdir -p ${HomePath}/METAregresults/plots

#### Download Reference files
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