###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(meta))
# suppressPackageStartupMessages(library(gap))
suppressPackageStartupMessages(library(dplyr))
options(bitmapType = "cairo")

#################################functions####################################
source ("/lustre/home/sww208/GoDMC/GeneticAssociationAnalysis/Rscripts/function_METAL_forestplot.R")

#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
MetaAll <- arguments[2]
MetaTbl <- arguments[3]
rsid <- arguments[4]
phenotype <- arguments[5]
targetSNPs <- arguments[6:length(arguments)]


setwd(setpath)
message("Reading in the files")
MetaAll = fread(MetaAll, header = T, data.table=F, stringsAsFactors=F)
colnames(MetaAll)[which(colnames(MetaAll) == "SNP")] <- "MarkerName"
colnames(MetaAll)[which(colnames(MetaAll) == "A1")] <- "EFFECT_ALLELE"
colnames(MetaAll)[which(colnames(MetaAll) == "A2")] <- "REFERENCE_ALLELE"
colnames(MetaAll)[which(colnames(MetaAll) == "P")] <- "PVAL"
colnames(MetaAll)[which(colnames(MetaAll) == "AF1")] <- "CODE_ALL_FQ"
MetaAll$prot <- phenotype

MetaTbl = fread(MetaTbl, header = T, data.table=F, stringsAsFactors=F, sep = " ")
MetaTbl$Chromosome <- str_split_fixed(MetaTbl$MarkerName, ":", 2)[,1]
MetaTbl$Position <- str_split_fixed(str_split_fixed(MetaTbl$MarkerName, ":", 2)[,2], "_",2)[,1]
MetaTbl$SNPout <- paste0("Chr",MetaTbl$Chromosome, "_", MetaTbl$Position)
MetaTbl$prot <- phenotype

message("Subset the Rsid file to include only the target SNPs")
message("Target SNPs: ", targetSNPs)
SNPrsid = fread(rsid, header = T, data.table=F, stringsAsFactors=F)
SNPrsid = subset(SNPrsid, MarkerName %in% targetSNPs)

# subset the MetaTbl and MetaAll to inlcude only the top SNPs
Tblinput = subset(MetaTbl, MarkerName %in% SNPrsid$MarkerName)
Allinput = subset(MetaAll, MarkerName %in% SNPrsid$MarkerName)
print(Tblinput)
# print(Allinput)
# print(SNPrsid)

METAL_forestplot_update(Tblinput, Allinput, SNPrsid,
                 digits.TE=2,digits.se=2,
                 col.diamond="green",col.inside="black", col.square="black", 
                 split = TRUE,
                 package="meta",method="REML",outpath="./plots/")