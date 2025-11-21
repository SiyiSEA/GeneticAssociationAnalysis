#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/09LDSCplot-%j.out
#SBATCH --error=JobReports/09LDSCplot-%j.err
#SBATCH --job-name=09LDSCplot
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks=16
#SBATCH --time=0-10:00:00
#SBATCH --array=8

module purge
source $1
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ${GENECORresPath}/logs/logs06cLDSCplots_log${timetemp}.log)
# cd ${GENECORresPath}/GENECOR_${phenotype}/plots || exit 1


if [ ${SLURM_ARRAY_TASK_ID} == 8 ]; then
    echo "Generating SNP heritability plot for all phenotypes..."
    cd ${GENECORresPath}/plots
    echo "Phenotype SNPHeritability SE" > HerOverall.txt
    for phenotype in DNAmAgeSD DNAmAgessSD PhenoAgeSD PhenoAgessSD DunedinPACESD DunedinPACEssSD gwas_smoking
    do
        her=$(tail ${LDSCresPath}/logs/logs06LDSCres_*_${phenotype}.log | grep "h2" | awk '{print $5}')
        se=$(tail ${LDSCresPath}/logs/logs06LDSCres_*_${phenotype}.log | grep "h2" | awk '{print $6}' | sed 's/[()]//g')
        echo "$phenotype $her $se" >> HerOverall.txt
    done

    echo "SNP heritability results collected:"
    Rscript ${RscriptsPath}/plot_SNPher.R \
                        ${GENECORresPath}/plots
                        ${GENECORresPath}/plots/HerOverall.txt \
                        "All"

    echo "Generating genetic correlation plots among phenotypes..."
    awk '/=========The Estimated Genetic Correlation across all the phenotypes=======/{flag=1} flag' ${LDSCresPath}/logs/logs06aLDSCres.log > ${GENECORresPath}/plots/GenCorOverall.txt
    Rscript ${RscriptsPath}/plot_GeneCorrClocks.R \
                        ${GENECORresPath}/plots
                        ${GENECORresPath}/plots/GenCorOverall.txt

fi



#### Check the Number of SNPs in the filtered results
echo " "
if [ ${SLURM_ARRAY_TASK_ID} -lt 8 ]; then
    echo "==========================Information of LDSC plots======================"

    if [ -f ${GENECORresPath}/plots/OverAll_her.png ]; then
        echo "SNP heritability plot for all phenotypes generated successfully."
    else
        echo "Error: SNP heritability plot for all phenotypes generation failed."
        exit 1
    fi

    if [ -f ${GENECORresPath}/plots/OverAllClocks_GeneCorr.png ]; then
        echo "Genetic correlation plots among phenotypes generated successfully."
    else
        echo "Error: Genetic correlation plots among phenotypes generation failed."
        exit 1
    fi

    echo "============================================================================"
    
    echo "Successfully completed the LDSC plots for the data under the $DataPath."
fi