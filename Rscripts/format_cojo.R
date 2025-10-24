###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))
options(bitmapType = "cairo")


###############################functions###################################



#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
cmafile <- arguments[2]
out <- arguments[3]

setwd(setpath)

cma = fread(cmafile, header = T)

cmaAA <- as.data.frame(str_split_fixed(cma$SNP, '_', 3))
colnames(cmaAA) = c("Chr:BP", "A1", "A2")
cmaAA$SNP = cma$SNP
# Join and assign altA
result <- cma %>%
  left_join(cmaAA, by = "SNP") %>%
  mutate(altA = ifelse(refA == A1, A2, A1))  # if refA==A1 use A2 else A1

result = result[,-c("Chr:BP","A1","A2")]
fwrite(result, file = paste0(out,"/LocusZoom.Cojo"), quote = F, sep = "\t", col.names = T, row.name = F )
message(paste0("Please check the file name", paste0(phenotype,"/LocusZoom.Cojo")))