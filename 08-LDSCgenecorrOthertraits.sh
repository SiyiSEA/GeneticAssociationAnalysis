#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/08LSDCgenecorrOthertraits-%j.out
#SBATCH --error=JobReports/08LSDCgenecorrOthertraits-%j.err
#SBATCH --job-name=08LSDCgenecorr
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks=16
#SBATCH --time=0-10:00:00
#SBATCH --array=2-6

module purge
source "/gpfs/ts0/shared/software/Miniconda3/23.5.2-0/etc/profile.d/conda.sh"
conda activate ldsc
LDSC="/lustre/home/sww208/Software/ldsc"

source $1 ${SLURM_ARRAY_TASK_ID}
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ${GENECORresPath}/logs/logs08GENECORres_log${timetemp}_${phenotype}.log)
cd ${GENECORresPath}/GENECOR_${phenotype} || exit 1

TraitPath="/lustre/home/sww208/GoDMC/GADatasets/OtherGWASstats"

epiAge="${LDSCresPath}/LDSC_${phenotype}/${phenotype}_byR.sumstats.gz"

for ID in CRP_2022 INT_2017 NEA_2021 CEA_2021;
do
    trait="${TraitPath}/$ID/$ID.sumstats.gz"
    outFile="${phenotype}_vs_${ID}"

    ${LDSC}/ldsc.py \
        --rg ${epiAge},${trait} \
        --ref-ld-chr ${HomePath}/Resources/baselineLD_v2.2/baselineLD. \
        --w-ld-chr ${HomePath}/Resources/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
        --out ${outFile}
done

# information
echo " "
echo "=========The Estimated Genetic Correlation between ${phenotype} and below traits======="
for ID in CRP_2022 INT_2017 NEA_2021 CEA_2021;
do
    outFile="${phenotype}_vs_${ID}"
    echo "${phenotype} and ${ID} ----------------------------------------------"

    grep "Genetic Correlation:" ${outFile}.log
    grep "Z-score:" ${outFile}.log
    grep "P:" ${outFile}.log
done
echo "============================================================================"

