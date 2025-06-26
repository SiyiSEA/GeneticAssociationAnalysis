#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/02QCMETA-matched-%a.out
#SBATCH --error=JobReports/02QCMETA-matched-%a.err
#SBATCH --job-name=02QCMETAmatched
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --ntasks=16
#SBATCH --time=0-20:00:00
#SBATCH --array=8

####################################################################################################################
# This script is for quality control of meta-analysis results on all the GWAS results collected from the GoDMC II.
# It merges the sample size scheme and standard error scheme results,
# and generates the filtered results for further analysis.
# The script also generates Manhattan and QQ plots for the filtered results.
#####################################################################################################################

module purge
module load R/4.2.1-foss-2022a
source ./config_matched "${SLURM_ARRAY_TASK_ID}"
timetemp=$(date -u +%Y-%m-%d_%H-%M)
exec &> >(tee ${METAresPath}/logs/02QCMETA-log${timetemp}_${phenotype}.log)

#### Merge function
Merge_SSSE_META() {
    echo "Merging ${1}_SE1.tbl and ${1}_SS1.tbl files..."

    phenotypeSE=${1}_SE1.tbl
    phenotypeSS=${1}_SS1.tbl

    cd ${METAresPath} || exit
    # Extract the FreqSE, Standard Error from SE scheme 
    awk -F'\t' '{print $1,$5,$8,$9,$10}' ./StandErrScheme/${phenotypeSE} > ./MergedMeta/${phenotypeSE}
    sed -i 's/FreqSE/FreqSE_SE/g' ./MergedMeta/${phenotypeSE}
    sed -i 's/P-value/P-value_SE/g' ./MergedMeta/${phenotypeSE}

    # Sort the files
    (head -n 1 ./SampleScheme/${phenotypeSS} && tail -n +2 ./SampleScheme/${phenotypeSS} | sort -k1,1) > ./MergedMeta/"${phenotypeSS}".sorted
    (head -n 1 ./MergedMeta/${phenotypeSE} && tail -n +2 ./MergedMeta/${phenotypeSE} | sort -k1,1) > ./MergedMeta/"${phenotypeSE}".sorted

    # Build the merged header
    h1=$(head -n1 ./MergedMeta/"${phenotypeSS}".sorted)
    h2=$(head -n1 ./MergedMeta/"${phenotypeSE}".sorted | cut -d' ' -f2-)
    echo "$h1 $h2" > ./MergedMeta/${1}_merged.tbl

    # Build the merged body
    join <(tail -n +2 ./MergedMeta/"${phenotypeSS}".sorted) <(tail -n +2 ./MergedMeta/"${phenotypeSE}".sorted) >> ./MergedMeta/${1}_merged.tbl
    rm ./MergedMeta/"${phenotypeSS}".sorted
    rm ./MergedMeta/"${phenotypeSE}".sorted

    echo "Merged ${1}_merged.tbl file created successfully."
}

#### Filter function
Filter_forN() {
    echo "Filtering ${1}_merged.tbl file for N > mean(N)..."

    mergedTBL=${1}_merged.tbl

    cd ${METAresPath} || exit
    # Calculate the mean, min, and max of the N column (10th column)
    echo "The N of variants in ${mergedTBL} is:"
    awk 'NR > 1 {
    sum += $10
    count++
    if (min == "" || $10 < min) min = $10
    if ($10 > max) max = $10
    } END {
    mean = sum / count
    print "Min:", min
    print "Max:", max
    print "Mean:", mean
    }' ./MergedMeta/${mergedTBL}

    # Input the mean value based on the output
    # echo "Please input the mean value from the above output:"
    # read -p 'Mean is ' mean
    
    # Filter the merged table based on the mean value
    awk -F' ' -v m=3000 'NR==1 || $10 >= m' ./MergedMeta/${mergedTBL} > ./MergedMeta/${mergedTBL}.filteredN.temp
    awk -F' ' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' ./MergedMeta/${mergedTBL}.filteredN.temp > ./MergedMeta/${mergedTBL}.filteredN
    rm ./MergedMeta/${mergedTBL}.filteredN.temp

    echo "Filtered ${1}_merged.tbl.filteredN file created successfully."

}

if [ ${SLURM_ARRAY_TASK_ID} -lt 8 ]; then

    #### Merge and filter for smoking
    Merge_SSSE_META ${phenotype}
    Filter_forN ${phenotype}

    #### Filter the merged table based on the freSE, HetQ, HetI and HetP values
    echo "Filtering ${phenotype}_merged.tbl.filteredN file for freSE, HetQ, HetI and HetP values..."
    Rscript ${RscriptsPath}/filter_meta.R \
                            ${METAresPath}/FilteredMeta \
                            ${METAresPath}/MergedMeta/${phenotype}_merged.tbl.filteredN \
                            ${phenotype}

    # Output files in the ${METAresPath}/FilteredMeta:
    # 1. _sigSNPlist.txt: list of significant SNPs with P < 1e-5;
    # 2. _input.cojo: input file for COJO analysis; ["MarkerName", "Allele1", "Allele2","Freq1","Effect","StdErr","P-value","N","Weight","P-value_SE"]
    # 3. .LocusZoom: LocusZoom file for visualization; ["Chr","BP","MarkerName", "Allele1", "Allele2","Freq1","Effect","StdErr","P-value","N","Weight","P-value_SE"]
    # 4. plots of I^2 before and after filtering;

    echo "Don't forget to upload the *.LocusZoom file into LocusZoom server for visualization."
    echo "Successfully completed the QC for ${phenotype}."
    echo "please set the array task id == 8 and re-run the script"

fi

#### Manhattan plot and QQ plot
if [ ${SLURM_ARRAY_TASK_ID} == 8 ]; then
    for phenotype in "DNAmAge" "PhenoAge" "DunedinPACE" # "gwas_smoking"
    do
        echo "Generating plots for ${phenotype}..."
        Rscript ${RscriptsPath}/plot_meta.R \
                            ${METAresPath}/FilteredMeta \
                            ${METAresPath}/FilteredMeta/${phenotype}.LocusZoom \
                            ${phenotype}

    done
fi

#### Check the Number of SNPs in the filtered results
if [ ${SLURM_ARRAY_TASK_ID} -lt 8 ]; then
    echo "==========================Information of ${phenotype}======================"
    echo "The total number of the SNPs in the raw SS meta-analysis of ${phenotype} is "
    wc -l ${HomePath}/METALresults/SampleScheme/${phenotype}_SS1.tbl
    echo "How many bad SNPs doesn't match to the 1000G eur reference file? "
    wc -l ${HomePath}/METALresults/StandErrScheme/${phenotype}_SE1.tbl
    echo "The total number of the SNPs in the merged meta-analysis of ${phenotype} is "
    wc -l ${METAresPath}/MergedMeta/${phenotype}_merged.tbl
    echo "The total number of the SNPs in the filterN meta-analysis of ${phenotype} is "
    wc -l ${METAresPath}/MergedMeta/${phenotype}_merged.tbl.filteredN
    echo "The total number of the SNPs in to plot for ${phenotype} is "
    wc -l ${METAresPath}/FilteredMeta/${phenotype}.LocusZoom
    echo "============================================================================"
fi