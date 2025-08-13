#!/usr/bin/env bash

#### Tools needed
# bcftools https://www.htslib.org/download/
# R
# GCTA-1.94.1 
# ldsc
# gemma-0.98.5

source $1
#### Set folders
mkdir -p "${METAresPath}"
mkdir -p "${METAresPath}/logs"
mkdir -p "${METAresPath}/SampleScheme"
mkdir -p "${METAresPath}/StandErrScheme"
mkdir -p "${METAresPath}/MergedMeta"

mkdir -p "${METAresPath}/FilteredMeta"
mkdir -p "${METAresPath}/FilteredMeta/plots"

mkdir -p "${COJOresPath}"
mkdir -p "${COJOresPath}/logs"

mkdir -p "${MstatresPath}"
mkdir -p "${MstatresPath}/logs"
mkdir -p "${MstatresPath}/plots"

mkdir -p "${METAregrePath}"
mkdir -p "${METAregrePath}/logs"
mkdir -p "${METAregrePath}/plots"

mkdir -p "${LDSCresPath}"
mkdir -p "${LDSCresPath}/logs"
mkdir -p "${LDSCresPath}/plots"



for phenotype in DNAmAgeSD DNAmAgessSD PhenoAgeSD PhenoAgessSD DunedinPACESD DunedinPACEssSD gwas_smoking; do
    mkdir -p "${COJOresPath}/Cojo_${phenotype}"
    mkdir -p "${MstatresPath}/Mstat_${phenotype}"
    mkdir -p "${MstatresPath}/Mstat_${phenotype}/plots"
    mkdir -p "${METAregrePath}/METAreg_${phenotype}"
    mkdir -p "${LDSCresPath}/LDSC_${phenotype}"
done

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
if [ -f ${ScriptsPath}/Resources/1000Geur/1000G_hg19_eur.bim ]; then
    echo "1000G eur reference files are exist."
else
    echo "1000G eur reference files does not exist. Please check the file path."
    exit 1
fi

if [ -f ${HomePath}/Resources/cohort_lambda.tsv ]; then
    echo "Cohort lambda file exists."
else
    echo "Cohort lambda file does not exist. Please check the file path."
    exit 1
fi

#if [ -f ${HomePath}/Resources/cohort_age.tsv ]; then
#    echo "Cohort age file exists."
#else
#    echo "Cohort age file does not exist. Please check the file path."
#    exit 1
#fi

if [ -f ${HomePath}/Resources/rsid.txt ]; then
    echo "rsid.txt file exists."
else
    echo "rsid.txt file does not exist. Please check the file path."
    exit 1
fi

if [ -f ${HomePath}/Resources/w_hm3.snplist ]; then
    echo "w_hm3.snplist file exists."
else
    cp /lustre/projects/Research_Project-MRC190311/references/LDScore/resources/w_hm3.snplist ${HomePath}/Resources/
fi

if [ -f ${HomePath}/Resources/baselineLD_v2.2 ]; then
    echo "baselineLD_v2.2 file exists."
else
    cp -r /lustre/projects/Research_Project-MRC190311/references/LDScore/resources/1000_genomes_phase3/baselineLD_v2.2 ${HomePath}/Resources/
fi

if [ -f ${HomePath}/Resources/1000G_Phase3_weights_hm3_no_MHC ]; then
    echo "1000G_Phase3_weights_hm3_no_MHC file exists."
else
    cp -r /lustre/projects/Research_Project-MRC190311/references/LDScore/resources/1000_genomes_phase3/1000G_Phase3_weights_hm3_no_MHC ${HomePath}/Resources/
fi



if [ -f ${HomePath}/Resources/common_all_ref_hg37.BED.uni.final ]; then
    echo "filtered_w_hm3.BED file exists"
else
    # Download the hg37 vcf file from NCBI

    cd ${HomePath}/Resources/ || exit
    wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
    gunzip common_all_20180423.vcf.gz
    vcf2bed --max-mem=50G < common_all_20180423.vcf > common_all_hg37.BED
    awk '{print $1, $2, $3, $4}' common_all_hg37.BED > common_all_ref_hg37.BED.temp
    awk '!/Y/' common_all_ref_hg37.BED.temp > common_all_ref_hg37.BED.temp2
    sed -e '1i Chr Start End SNP' common_all_ref_hg37.BED.temp2 > common_all_ref_hg37.BED.temp3

    # only care about the chr:bp format
    echo "Extracting SNPs from w_hm3.snplist and common_all_ref_hg37.BED.temp3"
    awk 'NR==FNR {snp_map[$1]; next} 
     FNR==1 || $4 in snp_map' w_hm3.snplist common_all_ref_hg37.BED.temp3 > common_all_ref_hg37.BED.dup.temp4 # 1460900
    
    echo "Number of SNPs in the common_all_ref_hg37.BED.dup.temp4: "
    wc -l common_all_ref_hg37.BED.dup.temp4 # 1460900

    echo "Number of duplicated SNPs in common_all_ref_hg37.BED.dup.temp4: "
    awk '{print $4}' common_all_ref_hg37.BED.dup.temp4 | uniq -d | wc -l # 229024

    echo "Extracting unique SNPs from common_all_ref_hg37.BED.dup.temp4"
    head -n 1 common_all_ref_hg37.BED.dup.temp4 > common_all_ref_hg37.BED.uni.final
    tail -n +2 common_all_ref_hg37.BED.dup.temp4 | sort -k1,1n -k2,2n | uniq >> common_all_ref_hg37.BED.uni.final #1213746
    echo "Number of SNPs in the common_all_ref_hg37.BED.uni.final: "
    wc -l common_all_ref_hg37.BED.uni.final # 1213746
    echo "Successfully created common_all_ref_hg37.BED.uni.final file."
    mv common_all_ref_hg37.BED.temp3 common_all_ref_hg37.BED
    rm common_all_ref_hg37.BED.temp* common_all_ref_hg37.BED.dup.temp4
fi
