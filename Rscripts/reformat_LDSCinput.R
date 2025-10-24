###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
options(bitmapType = "cairo")


#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
filename <- arguments[2]
BEDfile <- arguments[3]
whm3snplist <- arguments[4]
phenotype <- arguments[5]

setwd(setpath)
message("Reading in files")
METAresult = fread(filename, header = T, data.table=F, stringsAsFactors=F)
METAresult[grep('X', METAresult[,c("Chr")], ignore.case = TRUE),c("Chr")] = as.numeric(23)
METAresult$Chr = as.numeric(METAresult$Chr)
BED = fread(BEDfile, header = T, data.table=F, stringsAsFactors=F)
snplist = fread(whm3snplist, header = T, data.table=F, stringsAsFactors=F)

METAresult$"CHR:POS" = paste0(METAresult$Chr, ":", METAresult$BP)
BED$"CHR:POS" = paste0(BED$Chr, ":", BED$End)

METAsubset = subset(METAresult, METAresult$"CHR:POS" %in% BED$"CHR:POS")
message(paste("There are", length(METAsubset$"CHR:POS"),"SNPs","out of",length(snplist$SNP),"have been successfully converted into rsid."))

METAmerged = merge(METAsubset, BED[,c("SNP", "CHR:POS")], by.x = "CHR:POS", by.y = "CHR:POS")

METAout = METAmerged[,c("Chr","SNP", "BP","CHR:POS","Freq1","Allele1","Allele2","Zscore","N", "P-value","Direction")]
METAout = subset(METAout, METAout$Freq1 > 0.01)

# A1 as the reference (genome build) allele and A2 as the effect allele need to be reversed to match LDSC and the zscore needs to be recalculated
# where the Z score is calcuated by metal which is hard to recalcualate
message(paste0("Please check the formatted output:",setpath,"/", phenotype,".mungeinput"))
write.table(METAout, file = paste0(phenotype,".mungeinput"), row.names = F, quote = F, sep = "\t")

