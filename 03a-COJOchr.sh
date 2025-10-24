#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH -A Research_Project-MRC190311 # research project to submit under.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --output=JobReports/03aCOJO-%j-%a.out
#SBATCH --error=JobReports/03aCOJO-%j-%a.err
#SBATCH --job-name=03aCOJO
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
exec &> >(tee "${COJOresPath}"/logs/03aCOJO_log"${timetemp}"_"${phenotype}".log)

cd ${COJOresPath}/Cojo_"${phenotype}" || exit 1
cojoinput="${METAresPath}"/FilteredMeta/"${phenotype}"_input.cojo

for chr in {1..22}; do
    echo "Looking into the chromosome $chr ======================================================================="
    gcta-1.94.1 --bfile ${ScriptsPath}/Resources/1000Geur/1000G_hg19_eur \
                --chr "$chr" \
                --maf 0.01 \
                --cojo-file "$cojoinput" \
                --cojo-slct \
                --cojo-p 5e-8 \
                --out "${phenotype}"_"${chr}"

    if [ -f "${phenotype}"_"${chr}".jma.cojo ];
    then

        if [ "$nSignals" -lt 2 ];
        then
            echo "will skip the Step 1 (plotting) as there are only one signal"
            echo "will skip the Step 2 (condition analysis) as there are only one signal"
            echo "Step 3: Formatting the cojo output for LocusZoom"
            Rscript ${RscriptsPath}/format_cojo.R \
                    "${COJOresPath}"/Cojo_"${phenotype}" \
                    "${phenotype}"_"${chr}".jma.cojo \
                    "${phenotype}"_"${chr}"/03aLocusZoom
        else
            echo "Step 1: Plotting the LD for the selected Signals"
            # Rscript ${RscriptsPath}/plot_ld.R \
            #                 "${COJOresPath}"/Cojo_"${phenotype}" \
            #                 "${phenotype}_${chr}".ldr.cojo \
            #                 "${phenotype}_${chr}"
            if [[ "$chr" -eq 6 && ( "$phenotype" == "DNAmAgeSD" || "$phenotype" == "DNAmAgessSD" ) ]]; then
                echo "Prepare the files for conditional analysis for chr 6 for DNAmAgeSD and DNAmAgessSD"
                echo "extract the indenpendent SNPs in Chr6"
                Signals=("6:18104322_A_G" "6:91012867_C_T" "6:38629697_A_G" "6:2776087_A_C")
            else
                Signals=$(awk 'NR>1' "${phenotype}"_"${chr}".jma.cojo | sort -k8,8g | awk '{print $2}')
            fi

            nSignals=$(tail -n +2 "${phenotype}"_"${chr}".jma.cojo | wc -l)
            TopSignal=$(awk 'NR>1' "${phenotype}"_"${chr}".jma.cojo | sort -k8,8g | awk '{print $2}' | head -1)

            # sort the signal list
            awk 'NR>1' "${phenotype}"_"${chr}".jma.cojo | sort -k8,8g | awk '{print $2}' > "${phenotype}"_"${chr}".signals.list
            awk 'NR>1' "${phenotype}"_"${chr}".jma.cojo | sort -k8,8g | head -1 | awk '{print $2}' > "${phenotype}"_"${chr}".top.signals
            
            echo "The top Sig SNPs of $chr is $TopSignal  -------------------------------------------------"
            echo "There are $nSignals have been selected by the COJO -------------------------------------------------"


            echo "Step 2: Runing COJO for the top 5 signals while condition out the rest of the signal"
            for eachSignal in $Signals;
            do
                echo "For $eachSignal"
                grep -v "$eachSignal" "${phenotype}"_"${chr}".signals.list > out.txt
                # gcta-1.94.1 --bfile ${ScriptsPath}/Resources/1000Geur/1000G_hg19_eur\
                #     --cojo-file "$cojoinput" \
                #     --chr "$chr" \
                #     --maf 0.01 \
                #     --cojo-cond ./out.txt \
                #     --cojo-slct \
                #     --cojo-p 5e-8 \
                #     --out "${phenotype}"_"${eachSignal}"

                echo "Step 3: Formatting the cojo output on every singals for LocusZoom"
                echo "For "${phenotype}"_"${eachSignal}".cma.cojo"
                # Rscript ${RscriptsPath}/format_cojo.R \
                #     "${COJOresPath}"/Cojo_"${phenotype}" \
                #     "${phenotype}"_"${eachSignal}".cma.cojo \
                #     "${phenotype}"_"${eachSignal}"
            done
        fi

    else 
        echo "Skip as no signal."
    fi
    echo " "

    # rm out.txt
    # rm *.freq.badsnps
    # rm "${phenotype}"_"${chr}":*.jma.cojo
done
cat "${phenotype}"_*.signals.list > "${phenotype}".signals.temp
sort -n "${phenotype}".signals.temp | uniq > "${phenotype}".signals
cat "${phenotype}"_*.top.signals > "${phenotype}".top.signals.temp
sort -n "${phenotype}".top.signals.temp | uniq > "${phenotype}".top.signals
# rm *temp

echo " "
echo "Please upload the file into the LocusZoom by following the instruction:"
echo "1. Chromosome --> Chr"
echo "2. Position --> bp"
echo "3. Ref allele --> refA"
echo "4. Alt allele --> altA"
echo "5. p-value column --> pC"
echo "6. Std.Err --> bC_se"
echo "7. Effect allele --> Ref"
echo "8. Frequency --> freq"

# information
echo " "
echo "==========================Information of ${phenotype}======================"
echo "The total number of Sig SNPs detected by COJO of ${phenotype} is "
sort -n "${phenotype}".signals | uniq | wc -l
echo "How many chr has significantly SNPs?"
sort -n "${phenotype}".top.signals | uniq | wc -l
echo "A summary of number of Sig independent SNPs across each chr"
echo "chr    number"
cut -d: -f1 "${phenotype}".signals | sort -n | uniq -c | awk '{print $2 "\t" $1}'
echo "The outputs from this script needed to be checked"
echo "1. ${phenotype}_chr_LD.pdf: if the selected signals are related with each other."
echo "2. ${phenotype}.top.signals: a list of the top signals across every chr."
echo "3. ${phenotype}.signals: a list of the signals detected by Cojo."
echo "============================================================================"