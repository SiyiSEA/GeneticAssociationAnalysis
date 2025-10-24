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
#SBATCH --array=1-6

#####################################################################################################################

#####################################################################################################################

module purge
module load R/4.2.1-foss-2022a
source $1 "${SLURM_ARRAY_TASK_ID}"
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ${COJOresPath}/logs/03COJO_log${timetemp}_${phenotype}.log)

#### COJO 
# select independently associated SNPs
cd  "${COJOresPath}"/Cojo_"${phenotype}" || exit

# !!!!!! 1000Geur no ChrX!!!!
gcta-1.94.1 --bfile ${ScriptsPath}/Resources/1000Geur/1000G_hg19_eur \
            --cojo-file ${METAresPath}/FilteredMeta/"${phenotype}"_input.cojo \
            --cojo-slct \
            --cojo-p 5e-8 \
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

echo "Please check the LD plots for each ${phenotype} in the ${COJOresPath}/Cojo_${phenotype} directory."
echo "Making sure there are common signals for each pair of phenotypes."

# Pick the independent SNPs from the meta-analysis results
awk '{print $2}' $phenotype.jma.cojo | tail -n +2  > SignalCojoSNP.txt.temp

for snp in "${tartgetSNPs[@]}"; do
    echo "$snp" >> SignalCojoSNP.txt.temp
done

sort -u SignalCojoSNP.txt.temp > SignalCojoSNPlist.txt
rm SignalCojoSNP.txt.temp

# From each cohort, extract the SNPs in the SuperCojoSNP.txt
# the cohort sets can be modified in the config file
echo "Cohort fastGWA_N" > "${phenotype}"_chorts.fastGWA_N

cohortGWAs=$(grep "Input File" ${METAresPath}/SampleScheme/${phenotype}_SS1.tbl.info | cut -d':' -f2- | sed 's/^ *//')

for cohortGWA in $cohortGWAs
do
  echo "Grepping signal SNPs from ${cohortGWA}"
  cohortName=$(echo "${cohortGWA}" | awk -F'/' '{print $(NF-3)}')
  echo "Cohort name: ${cohortName}"
  awk -F $'\t' -v OFS=$'\t' -v name=${cohortName} 'NR==1  { print $0, "Study"; next } { print $0, name } ' ${cohortGWA} > temp.fastGWA
  
  # append any lines matching exactly one of the SNPs
  grep -Ff SignalCojoSNPlist.txt temp.fastGWA | tail -n +2 > "${cohortName}"_"${phenotype}"_temp.fastGWA

  samplesize=$(awk '{print $6}' temp.fastGWA | tail -n 1)
  echo "Counting the sample size: ${cohortName} ${samplesize}"
  echo "${cohortName} ${samplesize}" >> "${phenotype}"_chorts.fastGWA_N
done

head -n1 temp.fastGWA > "${MstatresPath}"/Mstat_"${phenotype}"/SignalSNPs.pre

cat *${phenotype}_temp.fastGWA >> "${MstatresPath}"/Mstat_"${phenotype}"/SignalSNPs.pre
mv "${phenotype}"_chorts.fastGWA_N "${MstatresPath}"/Mstat_"${phenotype}"/"${phenotype}"_chorts.fastGWA_N
rm *${phenotype}_temp*.fastGWA
rm temp.fastGWA

echo "Checking if the SignalSNPs.pre file is created successfully..."
num_lines=$(wc -l < "${MstatresPath}"/Mstat_"${phenotype}"/SignalSNPs.pre)
num_cohort=$(wc -l < "${MstatresPath}"/Mstat_"${phenotype}"/"${phenotype}"_chorts.fastGWA_N)
if [ $num_lines -gt $num_cohort ]; then
    echo "SignalSNPs.pre file created successfully with ${num_lines} lines."
else
    echo "Error: SignalSNPs.pre is failed !"
    exit 1
fi


# information
echo " "
echo "==========================Information of ${phenotype}======================"
echo "The total number of Sig SNPs detected by COJO of ${phenotype} is "
wc -l "${phenotype}".jma.cojo
echo "How many bad SNPs doesn't match to the 1000G eur reference file? "
wc -l ${phenotype}.freq.badsnps
echo "The output files from the script are:"
echo "1.${phenotype}.jma.cojo: independent SNPs selected by COJO"
echo "2.${phenotype}.cma.cojo: the rest of the SNPs"
echo "3.${phenotype}.ldr.cojo: LD correlation matrix between all pairwise SNPs listed in ${phenotype}.jma.cojo"
echo "4.SignalCojoSNPlist.txt: the list of independent SNPs selected by COJO and the target SNPs"
echo "5.Mstat_"${phenotype}"/SignalSNPs.pre: the list of signal SNPs from all cohorts"
echo "6.Mstat_"${phenotype}"/"${phenotype}"_chorts.fastGWA_N: the list of cohorts and their sample sizes"
echo "============================================================================"

echo "Successfully completed the COJO for ${phenotype}."