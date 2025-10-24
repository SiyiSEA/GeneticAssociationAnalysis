###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(MungeSumstats))
options(bitmapType = "cairo")

# need to be installed:
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
filename <- arguments[2]
phenotype <- arguments[3]

setwd(setpath)
message("Reading in files")
#filename = "/lustre/home/sww208/GoDMC/GADatasets/CorrectGWAS10/METALresults/FilteredMeta/DNAmAgeSD.LocusZoom"
METAresult = read_sumstats(filename)

print("Reading in the sumstats file")
print(head(METAresult))

METAresult = METAresult[,c("Chr","BP","MarkerName", "Allele1", "Allele2","Freq1","Effect",
"StdErr","P-value_SE","N","Zscore")]

METAresult$Chr[grep('23', METAresult$Chr, ignore.case = TRUE)] <- "X"
METAresult$Chr[grep('24', METAresult$Chr, ignore.case = TRUE)] <- "Y"

colnames(METAresult) <- c("CHR", "BP", "MarkerName", "effect_allele", "non_effect_allele", "FRQ",
                        "BETA", "SE", "P-value", "N", "Z")


reformatted <- MungeSumstats::format_sumstats(path = METAresult,
                                            ref_genome="GRCH37",
                                            on_ref_genome=TRUE,
                                            INFO_filter = 0.6,
                                            FRQ_filter = 0.01,
                                            snp_ids_are_rs_ids = FALSE,
                                            bi_allelic_filter = FALSE,
                                            flip_frq_as_biallelic = TRUE,
                                            nThread = 2,
                                            save_path = paste0(setpath, "/",phenotype,"_mungeinputbyR.tsv.gz"),
                                            log_mungesumstats_msgs = TRUE, 
                                            force_new=TRUE,
                                            save_format='LDSC',
                                            log_folder = setpath)



save(reformatted, file = paste0(setpath, "/", phenotype, "_mungeinputbyR.RData"))

message("Saving the reformatted results for LDSC")
message("Please check the file: ", paste0(setpath, "/",phenotype,"_mungeinputbyR.tsv.gz"))