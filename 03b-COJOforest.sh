#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/03bCOJO-%j-%a.out
#SBATCH --error=JobReports/03bCOJO-%j-%a.err
#SBATCH --job-name=03bCOJO
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --ntasks=16
#SBATCH --time=0-20:00:00
#SBATCH --array=1-2

#####################################################################################################################

#####################################################################################################################

module purge
module load R/4.2.1-foss-2022a
source $1 "${SLURM_ARRAY_TASK_ID}"
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ${COJOresPath}/logs/03bCOJO_log${timetemp}_${phenotype}.log)

#### Organise the cojo output for forest plot
cd  "${COJOresPath}"/Cojo_"${phenotype}"/03bForest || exit


# extract the Top signals from every cohorts.fastGWA
cohortGWAs=$(grep "Input File" ${METAresPath}/SampleScheme/${phenotype}_SS1.tbl.info | cut -d':' -f2- | sed 's/^ *//')

for cohortGWA in $cohortGWAs
do
  echo "Grepping signal SNPs from ${cohortGWA}"
  cohortName=$(echo "${cohortGWA}" | awk -F'/' '{print $(NF-3)}')
  echo "Cohort name: ${cohortName}"
  awk -F $'\t' -v OFS=$'\t' -v name=${cohortName} 'NR==1  { print $0, "Study"; next } { print $0, name } ' ${HomePath}/${cohortGWA} > temp.fastGWA
  # append any lines matching exactly one of the SNPs
  grep -F -f "${COJOresPath}"/Cojo_"${phenotype}"/${phenotype}.top.signals temp.fastGWA | tail -n +2 > "${cohortName}"_"${phenotype}"_temp.fastGWA
done

head -n1 temp.fastGWA > "${phenotype}".TopSiganlsStats
cat ./*"${phenotype}"_temp.fastGWA >> "${phenotype}".TopSiganlsStats
rm ./*temp*.fastGWA


#### Forest Plot based on the METAL results
grep -F -f "${COJOresPath}"/Cojo_"${phenotype}"/"${phenotype}".signals ${METAresPath}/MergedMeta/${phenotype}_merged.tbl.filteredN > "${phenotype}".tbl.pre
echo "MarkerName Allele1 Allele2 Freq1 FreqSE MinFreq MaxFreq N Direction Effect StdErr P-value" > "${phenotype}".forestSignal.tbl.pre
awk '{print $1, $2, $3, $4, $5, $6, $7, $10, $12, $9, $18, $11}' "${phenotype}".tbl.pre >> "${phenotype}".forestSignal.tbl.pre

rm "${phenotype}".tbl.pre

Rscript ${RscriptsPath}/METAL_forest.R \
            ${COJOresPath}"/Cojo_"${phenotype}/03bForest/\
            "${phenotype}".TopSiganlsStats \
            "${phenotype}".forestSignal.tbl.pre \
            ${HomePath}/Resources/rsid.txt \
            ${phenotype} \
            "${tartgetSNPs[@]}" 



# information
echo " "
echo "==========================Information of ${phenotype}======================"
# echo "The total number of Sig SNPs detected by COJO of ${phenotype} is "
# wc -l "${phenotype}".jma.cojo
# echo "How many bad SNPs doesn't match to the 1000G eur reference file? "
# wc -l ${phenotype}.freq.badsnps
# echo "The output files from the script are:"
# echo "1.${phenotype}.jma.cojo: independent SNPs selected by COJO"
# echo "2.${phenotype}.cma.cojo: the rest of the SNPs"
# echo "3.${phenotype}.ldr.cojo: LD correlation matrix between all pairwise SNPs listed in ${phenotype}.jma.cojo"
# echo "4.SignalCojoSNPlist.txt: the list of independent SNPs selected by COJO and the target SNPs"
# echo "5.Mstat_"${phenotype}"/SignalSNPs.pre: the list of signal SNPs from all cohorts"
# echo "6.Mstat_"${phenotype}"/"${phenotype}"_chorts.fastGWA_N: the list of cohorts and their sample sizes"
echo "============================================================================"

echo "Successfully completed the COJO for ${phenotype}."