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
#SBATCH --array=2-7

#####################################################################################################################

#####################################################################################################################

module purge
module load R/4.2.1-foss-2022a
source $1 "${SLURM_ARRAY_TASK_ID}"
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ${MstatresPath}/logs/04Mstat_log${timetemp}_${phenotype}.log)

cd "${MstatresPath}"/Mstat_"${phenotype}" || exit

##### M-stats
awk '{print $1, $5, $6, $7}' ${METAresPath}/FilteredMeta/"${phenotype}"_input.cojo > ./effectsize.temp
grep -Ff "${COJOresPath}"/Cojo_"${phenotype}"/SignalCojoSNPlist.txt ./effectsize.temp > ./effectsize_MstatSignal.pre
rm ./effectsize.temp

echo "Running the M-stats with Signal SNPs for ${phenotype}..."
Rscript ${RscriptsPath}/M_stats.R \
            "${MstatresPath}"/Mstat_"${phenotype}" \
            SignalSNPs.pre \
            effectsize_MstatSignal.pre \
            ${HomePath}/Resources/FastGWA_lambda.txt \
            ${phenotype} #\
            #"${tartgetSNPs[@]}"



#### Forest Plot based on the METAL results
grep -Ff "${COJOresPath}"/Cojo_"${phenotype}"/SignalCojoSNPlist.txt ${METAresPath}/MergedMeta/${phenotype}_merged.tbl.filteredN > "${phenotype}".tbl.pre
echo "MarkerName Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq N Direction Effect StdErr P-value" > "${phenotype}".forestSignal.tbl.pre
awk '{print $1, $2, $3, $4, $5, $6, $7, $10, $12, $9, $18, $11}' "${phenotype}".tbl.pre >> "${phenotype}".forestSignal.tbl.pre

rm "${phenotype}".tbl.pre

Rscript ${RscriptsPath}/METAL_forest.R \
            "${MstatresPath}"/Mstat_"${phenotype}" \
            SignalSNPs.pre \
            "${phenotype}".forestSignal.tbl.pre \
            "${COJOresPath}"/Cojo_"${phenotype}"/SignalCojoSNPlist.txt.rsid \
            ${phenotype} \
            ${phenotype}_SignificantSNPs.RData

##### Heterogeneity Analysis
echo "Running the Heterogeneity Analysis with Signal SNPs for ${phenotype}..."
# Rscript ${RscriptsPath}/plot_heter.R \
#             "${MstatresPath}"/Mstat_"${phenotype}" \
#              ${phenotype}_SignificantSNPs.RData \
#              ${phenotype}

# information
echo " "
echo "=====================The direction of topSNPs of ${phenotype}===================="
echo "Please check the forest plots for each significant SNP in the ${MstatresPath}/Mstat_${phenotype}/plots/ directory."
echo "The corresponding study of the direction is:"
grep "Input File" ${METAresPath}/SampleScheme/${phenotype}_SS1.tbl.info | cut -d':' -f2- | sed 's/^ *//' | awk -F'/' '{print $(NF-3)}'
echo "================================================================================="

echo "Successfully completed the M-stats and Heterogeneity Analysis for ${phenotype}."