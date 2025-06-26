#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/03COJO-%j-%a.out
#SBATCH --error=JobReports/03COJO-%j-%a.err
#SBATCH --job-name=03COJO
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
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ./COJOresults/logs/03COJO_log${timetemp}_${phenotype}.log)


#### COJO 
# select independently associated SNPs
cd  "${COJOresPath}"/Cojo_"${phenotype}" || exit

gcta-1.94.1 --bfile ${HomePath}/Resources/1000Geur/1000G_hg19_eur \
            --cojo-file ${METAresPath}/FilteredMeta/"${phenotype}"_input.cojo \
            --cojo-slct \
            --cojo-p 0.00001 \
            --out "${phenotype}"

# output files descriptions
# "${phenotype}".ldr.cojo: LD correlation matrix between all pairwise SNPs listed in "${phenotype}".jma.cojo
# "${phenotype}".jma.cojo: independent SNPs selected by COJO
# "${phenotype}".cma.cojo: the rest of the SNPs

# plot LD from COJO results
Rscript ${RscriptsPath}/plot_ld.R \
                        "${COJOresPath}"/Cojo_"${phenotype}" \
                        "${phenotype}".ldr.cojo \
                        "${phenotype}"


# # Version-1: Pick the SNPs based on the BP(+-5kb) of independent SNPs identified by COJO from metaed data
# singal_Chr=$(tail -n +2 "${phenotype}".jma.cojo | cut -f1 | sort -nu)
# for chr in $singal_Chr; do
#   echo "split the cma.cojo by Chromosome: $chr"
#   echo "SNP BP" > Chr"${chr}"_cma.cojo
#   awk '{print $2,$3}' "${phenotype}".cma.cojo | grep "^$chr:" >> Chr"${chr}"_cma.cojo
#   echo "split the jma.cojo by Chromosome: $chr"
#   echo "SNP BP" > Chr"${chr}"_jma.cojo
#   awk '{print $2,$3}' "${phenotype}".jma.cojo | grep "^$chr:" >> Chr"${chr}"_jma.cojo

#   Rscript ${RscriptsPath}/pick_SNPs.R \
#                 ${COJOresPath}/Cojo_"${phenotype}" \
#                 Chr"${chr}"_cma.cojo \
#                 Chr"${chr}"_jma.cojo
# done

# cat Chr*CojoSNPlist.txt > SuperCojoSNP.txt
# rm Chr*CojoSNPlist.txt
# rm Chr*.cojo

# Version-2: Pick the independent SNPs from the meta-analysis results
awk '{print $2}' $phenotype.jma.cojo | tail -n +2  > SignalCojoSNP.Mpre
awk '{print $2}' $phenotype.jma.cojo | tail -n +2  > SignalCojoSNP.forestpre.temp

if [ "$phenotype" = "DNAmAgeSD" ] || [ "$phenotype" = "DNAmAgessSD" ]; then
    {
        echo "6:18130918_C_T"
        echo "6:18120400_C_T"
    } >> SignalCojoSNP.forestpre.temp
elif [ "$phenotype" = "PhenoAgeSD" ] || [ "$phenotype" = "PhenoAgessSD" ]; then
    {
        echo "6:18130918_C_T"
        echo "6:18109815_C_G"
    } >> SignalCojoSNP.forestpre.temp
     
elif [ "$phenotype" = "DunedinPACESD" ] || [ "$phenotype" = "DunedinPACEssSD" ]; then
    {
        echo "14:74141429_A_G"
        echo "14:74176679_C_T"
        echo "12:125342972_A_G"
    } >> SignalCojoSNP.forestpre.temp

elif [ "$phenotype" = "gwas_smoking" ]; then
    {
        echo "3:128321226_A_G"
        echo "3:22585598_C_T"
        echo "15:75052495_C_T"
        echo "19:58776215_A_T"
    } >> SignalCojoSNP.forestpre.temp

else
    echo "No specific SNPs for forest plot in ${phenotype}."
    exit
fi

sort -u SignalCojoSNP.forestpre.temp > SignalCojoSNP.forestpre


# From each cohort, extract the SNPs in the SuperCojoSNP.txt
# the cohort sets can be modified in the config file
echo "Cohort fastGWA_N" > "${phenotype}"_chorts.fastGWA_N

cohorts=$(grep 'StatsAgeSmoke' ${HomePath}/METALresults/SampleScheme/${phenotype}_SS1.tbl.info | awk -F'/' '{ print $2 }')
for cohort in $cohorts
do
  echo "Grepping signal SNPs from ${cohort}"
  cohortGWA="${HomePath}/StatsAgeSmoke/${cohort}/${fastGWA}"
  awk -F $'\t' -v OFS=$'\t' -v name=${cohort} 'NR==1  { print $0, "Study"; next } { print $0, name } ' ${cohortGWA} > temp.fastGWA
  
  # append any lines matching exactly one of the SNPs
  grep -Ff SignalCojoSNP.Mpre temp.fastGWA | tail -n +2 > "${cohort}"_"${phenotype}"_tempMpre.fastGWA
  grep -Ff SignalCojoSNP.forestpre temp.fastGWA | tail -n +2 > "${cohort}"_"${phenotype}"_tempforestpre.fastGWA

  samplesize=$(awk '{print $6}' temp.fastGWA | tail -n 1)
  echo "Counting the sample size: ${cohort} ${samplesize}"
  echo "${cohort} ${samplesize}" >> "${phenotype}"_chorts.fastGWA_N
done

head -n1 temp.fastGWA > "${MstatresPath}"/Mstat_"${phenotype}"/SignalSNP.Mstat.pre
head -n1 temp.fastGWA > "${MstatresPath}"/Mstat_"${phenotype}"/SignalSNP.forest.pre

cat *${phenotype}_tempMpre.fastGWA >> "${MstatresPath}"/Mstat_"${phenotype}"/SignalSNP.Mstat.pre
cat *${phenotype}_tempforestpre.fastGWA >> "${MstatresPath}"/Mstat_"${phenotype}"/SignalSNP.forest.pre
mv "${phenotype}"_chorts.fastGWA_N "${MstatresPath}"/Mstat_"${phenotype}"/"${phenotype}"_chorts.fastGWA_N
rm *${phenotype}_temp*.fastGWA


echo "==========================Information of ${phenotype}======================"
echo "The total number of Sig SNPs detected by COJO of ${phenotype} is "
wc -l "${phenotype}".jma.cojo
echo "How many bad SNPs doesn't match to the 1000G eur reference file? "
wc -l ${phenotype}.freq.badsnps
echo "============================================================================"

echo "Successfully completed the COJO for ${phenotype}."