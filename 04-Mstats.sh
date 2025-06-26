#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/04Mstats-%j.out
#SBATCH --error=JobReports/04Mstats-%j.err
#SBATCH --job-name=04Mstats
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --ntasks=16
#SBATCH --time=0-20:00:00
#SBATCH --array=7

#####################################################################################################################

#####################################################################################################################

module purge
module load R/4.2.1-foss-2022a
source ./config "${SLURM_ARRAY_TASK_ID}"
# source ./config 7
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ./Mstatresults/logs/04Mstat_log${timetemp}_${phenotype}.log)

cd "${MstatresPath}"/Mstat_"${phenotype}" || exit

##### M-stats
awk '{print $1, $5, $6}' ${METAresPath}/FilteredMeta/"${phenotype}"_input.cojo > ./effectsize.temp
grep -Ff "${COJOresPath}"/Cojo_"${phenotype}"/SignalCojoSNP.Mpre ./effectsize.temp > ./effectsize_MstatSignal.pre
rm ./effectsize.temp

jamcojo="${COJOresPath}"/Cojo_"${phenotype}"/${phenotype}.jma.cojo
head -n1 ${jamcojo} && tail -n +2 ${jamcojo} | sort -k8,8g > ${phenotype}.jma.cojo.sorted
head -n3 "${phenotype}".jma.cojo.sorted | awk '{print $2, $8}' > Top3SNPs.txt

echo "Running the M-stats with Signal SNPs for ${phenotype}..."
Rscript ${RscriptsPath}/M_stats.R \
            "${MstatresPath}"/Mstat_"${phenotype}" \
            SignalSNP.Mstat.pre \
            effectsize_MstatSignal.pre \
            Top3SNPs.txt \
            ${HomePath}/Resources/cohort_lambda.tsv \
            ${phenotype}


#### Forest Plot based on the METAL results
grep -Ff "${COJOresPath}"/Cojo_"${phenotype}"/SignalCojoSNP.forestpre ${METAresPath}/MergedMeta/${phenotype}_merged.tbl.filteredN > "${phenotype}".tbl.pre
echo "MarkerName Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq N Direction Effect StdErr P-value" > "${phenotype}".forestSignal.tbl.pre
awk '{print $1, $2, $3, $4, $5, $6, $7, $10, $12, $9, $18, $11}' "${phenotype}".tbl.pre >> "${phenotype}".forestSignal.tbl.pre

rm "${phenotype}".tbl.pre

Rscript ${RscriptsPath}/METAL_forest.R \
            "${MstatresPath}"/Mstat_"${phenotype}" \
            SignalSNP.forest.pre \
            "${phenotype}".forestSignal.tbl.pre \
            ${phenotype}


##### Heterogeneity Analysis
echo "Running the Heterogeneity Analysis with Signal SNPs for ${phenotype}..."
Rscript ${RscriptsPath}/plot_heter.R \
            "${MstatresPath}"/Mstat_"${phenotype}" \
            "$metaRData"

# information
echo "=====================The direction of topSNPs of ${phenotype}===================="
echo "MarkerName Effect P-value Direction"
for t in "${tartgetSNPs[@]}"; do
    grep "$t" "${phenotype}".forestSignal.tbl.pre | awk '{print $1, $10, $12, $9}'
done
echo "The corresponding study of the direction is:"
cohorts=$(grep 'StatsAgeSmoke' ${HomePath}/METALresults/SampleScheme/${phenotype}_SS1.tbl.info | awk -F'/' '{ print $2 }')
echo $cohorts
echo "================================================================================="

echo "Successfully completed the M-stats and Heterogeneity Analysis for ${phenotype}."