#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/01META-%j-%a.out
#SBATCH --error=JobReports/01META-%j-%a.err
#SBATCH --job-name=01MediMETA
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --ntasks=16
#SBATCH --time=0-20:00:00
#SBATCH --array=1-7

#####################################################################################################################
# This script is for mate-analysis on all the GWAS results collected from the GoDMC II.
# Although there are six phenotypes, only the ${MetalPath}/metal_DNAmAgeSD.txt need to be updated manually,
# The rest of the metal scripts will be auto freshed up. This may around 4 hours for all and PCA...
#####################################################################################################################

module purge
module load R/4.2.1-foss-2022a
source $1 "${SLURM_ARRAY_TASK_ID}"
cd ${HomePath} || exit
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ${METAresPath}/logs/01META-log${timetemp}_${phenotype}.log)

#### Prepare METAL scripts
for study in ${AllCohortList[@]}; do
    echo "PROCESS ${DataPath}/${study}/results/10/DNAmAgeSD.fastGWA" >> ${MetalPath}/DNAmAgeSD_SAMPLE_header.txt
    echo "PROCESS ${DataPath}/${study}/results/10/DNAmAgeSD.fastGWA" >> ${MetalPath}/DNAmAgeSD_STERR_header.txt
done

cat ${MetalPath}/DNAmAgeSD_SAMPLE_header.txt ${MetalPath}/DNAmAgeSD_SAMPLE_tail.txt > ${MetalPath}/DNAmAgeSD_SAMPLE.txt
cat ${MetalPath}/DNAmAgeSD_STERR_header.txt ${MetalPath}/DNAmAgeSD_STERR_tail.txt > ${MetalPath}/DNAmAgeSD_STERR.txt

#### check the DNAmAgeSD.txt file
if [ ! -f ${MetalPath}/DNAmAgeSD_SAMPLE.txt ]; then
    echo "Error: ${MetalPath}/DNAmAgeSD.txt does not exist. Please check the file path."
    exit 1
fi
if [ ! -f ${MetalPath}/DNAmAgeSD_STERR.txt ]; then
    echo "Error: ${MetalPath}/DNAmAgeSD_SAMPLE.txt does not exist. Please check the file path."
    exit 1
fi



#### prepare for the meta scripts for sample size scheme
echo "Preparing meta-analysis scripts..."
sed "s/DNAmAgeSD/DNAmAgessSD/g" ${MetalPath}/DNAmAgeSD_SAMPLE.txt > ${MetalPath}/DNAmAgessSD_SAMPLE.txt

sed "s/DNAmAge/PhenoAge/g" ${MetalPath}/DNAmAgeSD_SAMPLE.txt > ${MetalPath}/PhenoAgeSD_SAMPLE.txt
sed "s/PhenoAgeSD/PhenoAgessSD/g" ${MetalPath}/PhenoAgeSD_SAMPLE.txt > ${MetalPath}/PhenoAgessSD_SAMPLE.txt

sed "s/DNAmAge/DunedinPACE/g" ${MetalPath}/DNAmAgeSD_SAMPLE.txt > ${MetalPath}/DunedinPACESD_SAMPLE.txt
sed "s/DunedinPACESD/DunedinPACEssSD/g" ${MetalPath}/DunedinPACESD_SAMPLE.txt > ${MetalPath}/DunedinPACEssSD_SAMPLE.txt


#### prepare for the meta scripts for standard error scheme
echo "Preparing meta-analysis scripts..."
sed "s/DNAmAgeSD/DNAmAgessSD/g" ${MetalPath}/DNAmAgeSD_STERR.txt > ${MetalPath}/DNAmAgessSD_STERR.txt

sed "s/DNAmAge/PhenoAge/g" ${MetalPath}/DNAmAgeSD_STERR.txt > ${MetalPath}/PhenoAgeSD_STERR.txt
sed "s/PhenoAgeSD/PhenoAgessSD/g" ${MetalPath}/PhenoAgeSD_STERR.txt > ${MetalPath}/PhenoAgessSD_STERR.txt

sed "s/DNAmAge/DunedinPACE/g" ${MetalPath}/DNAmAgeSD_STERR.txt > ${MetalPath}/DunedinPACESD_STERR.txt
sed "s/DunedinPACESD/DunedinPACEssSD/g" ${MetalPath}/DunedinPACESD_STERR.txt > ${MetalPath}/DunedinPACEssSD_STERR.txt


#### prepare for the meta scripts for smoking
if [ ${phenotype} == "gwas_smoking" ]; then
    echo "Preparing meta-analysis scripts for smoking..."
    sed "s/10/11/g" ${MetalPath}/DNAmAgeSD_SAMPLE.txt > ${MetalPath}/smoking_SAMPLE.temp
    sed "s/DNAmAgeSD/gwas_smoking/g" ${MetalPath}/smoking_SAMPLE.temp > ${MetalPath}/gwas_smoking_SAMPLE.txt

    sed "s/10/11/g" ${MetalPath}/DNAmAgeSD_STERR.txt > ${MetalPath}/smoking_STERR.temp
    sed "s/DNAmAgeSD/gwas_smoking/g" ${MetalPath}/smoking_STERR.temp > ${MetalPath}/gwas_smoking_STERR.txt

    rm ${MetalPath}/smoking_SAMPLE.temp
    rm ${MetalPath}/smoking_STERR.temp
fi

#### Execute meta-analysis for .fastGWA based on sample size scheme and standard error scheme
metal ${MetalPath}/${phenotype}_SAMPLE.txt
metal ${MetalPath}/${phenotype}_STERR.txt

#### Check the cohorts in the meta-analysis results
echo " "
echo "==========================Information of ${phenotype}====================="
echo "The Top SNPs for ${phenotype} are:"
grep "Smallest p-value" ${METAresPath}/logs/01META-log${timetemp}_${phenotype}.log

nSS=$(grep "Input File" ${METAresPath}/SampleScheme/${phenotype}_SS1.tbl.info | wc -l)
nSE=$(grep "Input File" ${METAresPath}/StandErrScheme/${phenotype}_SE1.tbl.info | wc -l)

echo "Number of cohorts in Sample Size Scheme: ${nSS}"
echo "Number of cohorts in Standard Error Scheme: ${nSE}"

if [ ${nSS} -ne ${nSE} ]; then
    echo "Error: The number of cohorts in Sample Size Scheme and Standard Error Scheme do not match."
    echo "Please check the meta-analysis results for ${phenotype}."
    exit 1
else
    echo "Successfully completed the meta-analysis for ${phenotype}."
fi

echo "=========================================================================="
