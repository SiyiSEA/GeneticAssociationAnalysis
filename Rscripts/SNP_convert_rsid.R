###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(colochelpR))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP155.GRCh37))
options(bitmapType = "cairo")


###############################functions###################################



#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
snplistfile <- arguments[2]
out <- arguments[3]

setwd(setpath)

SNPlist = read.table(snplistfile, header = F, stringsAsFactors = F)
colnames(SNPlist) = c("MarkerName")

SNPlist[c('CHR', 'BP_ALT_REF')] <- str_split_fixed(SNPlist$MarkerName, ':', 2)
SNPlist[c('BP', 'ALT_REF')] <- str_split_fixed(SNPlist$BP_ALT_REF, '_', 2)

snps = SNPlocs.Hsapiens.dbSNP155.GRCh37
SNPlist_rsid = colochelpR::convert_loc_to_rs(SNPlist,dbSNP=snps)

message("After converting, there are total ", nrow(SNPlist_rsid), " SNPs with rsid.")
message("There are ", length(duplicated(SNPlist_rsid$SNP) == TRUE), " duplicated SNPs after converting to rsid, removing...")
SNPlist_rsid = SNPlist_rsid[!duplicated(SNPlist_rsid$SNP),c("MarkerName","SNP") ]

if (dim(SNPlist_rsid)[1] == 0){
    message("Warning: NO RSID found after converting!")
}

fwrite(SNPlist_rsid, file = out, quote = F, sep = "\t", col.names = T, row.name = F )
message(paste0("Please check the file name ", out))