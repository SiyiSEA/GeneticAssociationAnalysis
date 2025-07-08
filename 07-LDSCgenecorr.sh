#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/05MetaRegre-%j-a%.out
#SBATCH --error=JobReports/05MetaRegre-%j-a%.err
#SBATCH --job-name=05MetaRegre
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --ntasks=16
#SBATCH --time=0-20:00:00
#SBATCH --array=1-7

#####################################################################################################################

#####################################################################################################################

module purge
module load R/4.2.1-foss-2022a
# source ./config "${SLURM_ARRAY_TASK_ID}"
source ./config 1
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ./METAregresults/logs/logs05METAregre_log${timetemp}_${phenotype}.log)



#### Dealing with the tsv for cohort description
if [ -f ${HomePath}/Resources/TableS1_cohortdescriptives.xlsx ]; then
    echo "TableS1_cohortdescriptives.xlsx exists, dealing with it..."
    Rscript ${RscriptsPath}/CohortDescription.R \
                ${HomePath}/Resources/ \
                ${HomePath}/Resources/TableS1_cohortdescriptives.xlsx \
                "${MstatresPath}"/Mstat_"${phenotype}"/${phenotype}_chorts.fastGWA_N
                
else
    echo "Error!!!TableS1_cohortdescriptives.xlsx does not exist, please upload it"
    exit 1
fi


##### Meta-regression
Rscript ${RscriptsPath}/Meta_regression.R \
            "${METAregrePath}"/METAreg_"${phenotype}" \
            "${MstatresPath}"/Mstat_"${phenotype}"/"${phenotype}"_MstatResults.RData \
            "${MstatresPath}"/Mstat_"${phenotype}"/${metaRData} \
            ${HomePath}/Resources/cohort_age.tsv \
            ${HomePath}/Resources/cohort_array.tsv \
            ${HomePath}/Resources/cohort_imputedpanel.tsv \
            ${phenotype}

echo "Successfully completed the M-stats and Heterogeneity Analysis for ${phenotype}."