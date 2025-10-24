#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/06LSDCALL-%j-%a.out
#SBATCH --error=JobReports/06LSDCALL-%j-%a.err
#SBATCH --job-name=06LSDCALL
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --ntasks=16
#SBATCH --time=0-10:00:00
#SBATCH --array=1-7

#####################################################################################################################

#####################################################################################################################

module purge
module load R/4.2.1-foss-2022a
source $1 "${SLURM_ARRAY_TASK_ID}"
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ${LDSCresPath}/logs/logs06LDSCres_log${timetemp}_${phenotype}.log)


cd ${LDSCresPath}/LDSC_${phenotype}/ || exit

# echo "Manually Reformatting the LDSC input files for ${phenotype}..."
# Rscript ${RscriptsPath}/reformat_LDSCinput.R \
#                             ${LDSCresPath}/LDSC_${phenotype} \
#                             ${METAresPath}/FilteredMeta/${phenotype}.LDSCinput \
#                             ${HomePath}/Resources/common_all_ref_hg37.BED.uni.final \
#                             ${HomePath}/Resources/w_hm3.snplist \
#                             ${phenotype}

module purge
source "/gpfs/ts0/shared/software/Miniconda3/23.5.2-0/etc/profile.d/conda.sh"
conda activate RGreatSquareRoot
source $1 "${SLURM_ARRAY_TASK_ID}"
echo "R package Reformatting the LDSC input files for ${phenotype}..."
Rscript ${RscriptsPath}/reformat_MungeSumstats.R \
                            ${LDSCresPath}/LDSC_${phenotype} \
                            ${METAresPath}/FilteredMeta/${phenotype}.LocusZoom \
                            ${phenotype}
conda deactivate

# output is *_mungeinputbyR.RData and *_mungeinputbyR.tsv.gz

# if [ -f ${phenotype}.mungeinput ] && [ -f ${phenotype}_mungeinputbyR.tsv.gz ]; then
#     echo "LDSC input file for ${phenotype} is generated successfully."
# else
#     echo "LDSC input file for ${phenotype} is not generated. Please check the previous steps."
#     exit 1
# fi

echo "Doing the munge_sumstats and h2 estimation for ${phenotype}..."
module purge
source "/gpfs/ts0/shared/software/Miniconda3/23.5.2-0/etc/profile.d/conda.sh"
conda activate ldsc
source $1 "${SLURM_ARRAY_TASK_ID}"
LDSC="/lustre/home/sww208/Software/ldsc"

cd ${LDSCresPath}/LDSC_${phenotype} || exit

# ${LDSC}/munge_sumstats.py \
#         --sumstats ${phenotype}.mungeinput \
#         --snp SNP \
#         --N-col N \
#         --a1 Allele1 \
#         --a2 Allele2 \
#         --frq Freq1 \
#         --p P-value \
#         --chunksize 50000 \
#         --signed-sumstats Zscore,0 \
#         --merge-alleles ${HomePath}/Resources/w_hm3.snplist \
#         --out ${phenotype}_manually

# ${LDSC}/ldsc.py \
#         --h2 ${phenotype}_manually.sumstats.gz \
#         --ref-ld-chr ${HomePath}/Resources/baselineLD_v2.2/baselineLD. \
#         --w-ld-chr ${HomePath}/Resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
#         --out ${phenotype}_manually_h2


###
${LDSC}/munge_sumstats.py \
        --sumstats ${phenotype}_mungeinputbyR.tsv.gz \
        --snp SNP \
        --N-col N \
        --a1 A1 \
        --a2 A2 \
        --frq FRQ \
        --p P \
        --chunksize 50000 \
        --signed-sumstats Z,0 \
        --merge-alleles ${HomePath}/Resources/w_hm3.snplist \
        --out ${phenotype}_byR

${LDSC}/ldsc.py \
        --h2 ${phenotype}_byR.sumstats.gz \
        --ref-ld-chr ${HomePath}/Resources/baselineLD_v2.2/baselineLD. \
        --w-ld-chr ${HomePath}/Resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
        --out ${phenotype}_byR_h2

# information
echo " "
# echo "=====================The Estimated Hertibitliy on ${phenotype}===================="
# echo "Using the LDSC input file ${phenotype}.mungeinput"
# echo "How many SNPs used to estimated the hertibility:"
# grep "Writing summary statistics" ${phenotype}_manually.log
# grep "SNPs remain" ${phenotype}_manually_h2.log
# grep "Total Observed scale" ${phenotype}_manually_h2.log
# grep "Lambda GC" ${phenotype}_manually_h2.log
# grep "Mean Chi" ${phenotype}_manually_h2.log
# grep "Intercept" ${phenotype}_manually_h2.log
# grep "Ratio" ${phenotype}_manually_h2.log
# echo "Any Warnings or Errors in the log file:"
# grep -i "error" ${phenotype}_manually_h2.log
# grep -i "warning" ${phenotype}_manually_h2.log
# echo "================================================================================="

echo "=====================The Estimated Hertibitliy on ${phenotype}===================="
echo "Using the LDSC input file ${phenotype}_mungeinputbyR.tsv.gz"
echo "How many SNPs used to estimated the hertibility:"
grep "Writing summary statistics" ${phenotype}_byR.log
grep "SNPs remain" ${phenotype}_byR_h2.log
grep "Total Observed scale" ${phenotype}_byR_h2.log
grep "Lambda GC" ${phenotype}_byR_h2.log
grep "Mean Chi" ${phenotype}_byR_h2.log
grep "Intercept" ${phenotype}_byR_h2.log
grep "Ratio" ${phenotype}_byR_h2.log
echo "Any Warnings or Errors in the log file:"
grep -i "error" ${phenotype}_byR_h2.log
grep -i "warning" ${phenotype}_byR_h2.log
echo "================================================================================="

echo "Successfully completed the M-stats and Heterogeneity Analysis for ${phenotype}."