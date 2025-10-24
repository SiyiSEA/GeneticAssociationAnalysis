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

module purge

timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ${GENECORresPath}/logs/logs09LDSCplots_log${timetemp}_${phenotype}.log)
cd ${GENECORresPath}/GENECOR_${phenotype}/plots || exit 1

# make plot for heritability

# orgaizing results grep from log file
if [ ${SLURM_ARRAY_TASK_ID} == 8 ]; then
    cd ${GENECORresPath}/plots
    echo "Phenotype SNPHeritability SE" > HerOverall.txt
    for phenotype in DNAmAgeSD DNAmAgessSD PhenoAgeSD PhenoAgessSD DunedinPACESD DunedinPACEssSD gwas_smoking
    do
        her=$(tail ${LDSCresPath}/logs/logs06LDSCres_log2025-09-15_10-17_${phenotype}.log | grep "h2" | awk '{print $5}')
        se=$(tail ${LDSCresPath}/logs/logs06LDSCres_log2025-09-15_10-17_${phenotype}.log | grep "h2" | awk '{print $6}' | sed 's/[()]//g')
        echo "$phenotype $her $se" >> HerOverall.txt
    done

    Rscript ${RscriptsPath}/plot_heter.R \
                        ${GENECORresPath}/plots
                        ${GENECORresPath}/plots/HerOverall.txt \
                        "All"

fi