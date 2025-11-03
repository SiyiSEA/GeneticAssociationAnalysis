###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(meta))
suppressPackageStartupMessages(library(dmetar))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggrepel))
options(bitmapType = "cairo")

#################################functions####################################
source ("/lustre/home/sww208/GoDMC/GeneticAssociationAnalysis/Rscripts/function_influence.analysis.R")

#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
targetSNPs <- arguments[2]
phenotype <- arguments[3]

setwd(setpath)
message("Reading in the files")
load(file = targetSNPs)
message("The significant SNPs to be analyzed for heterogeneity is: ", SigSNPvector[1])
chr = str_split_fixed(SigSNPvector[1],":",2)[1]
bp = str_split_fixed(str_split_fixed(SigSNPvector[1],"_",2)[1],":",2)[2]

metaobject= paste0("./plots/",phenotype,"-Chr",chr,"_",bp,".RData")
load(file=metaobject)

Outname = tools::file_path_sans_ext(metaobject)
sink(paste0(Outname, "_hemo_heter.record"))

message("------------Listing the outliers based on the upper and lower CI limits---------------")
print(find.outliers(mg))

# Outliers and Influence Cases
metaoutput_inf = InfluenceAnalysis_update(mg, random = TRUE, return.separate.plots = TRUE)
message("Plotting baujat plot for SNP ", Outname, " ----------------------------")
pdf(file = paste0(Outname,"_baujat.pdf"),width=5, height=5)
plot(metaoutput_inf, "baujat")
dev.off()
message("Plotting influence plot for ", Outname, " ----------------------------")
pdf(file = paste0(Outname,"_influence.pdf"), width = 15, height = 10)
plot(metaoutput_inf, "influence")
dev.off()
message("Plotting forestplot_es.pdf plot for ", Outname, " ----------------------------")
png(file = paste0(Outname,"_forestplot_es.png"), width = 600, height = 600)
plot(metaoutput_inf, "es")
dev.off()
message("Plotting forestplot_i2.pdf plot for ", Outname, " ----------------------------")
png(file = paste0(Outname,"_forestplot_i2.png"), width = 600, height = 600)
plot(metaoutput_inf, "i2")
dev.off()

# GOSH test for heterogeneity
# library(metafor)
# m.rma <- rma(yi = metaoutput$TE,
#             sei = metaoutput$seTE,
#             method = metaoutput$method.tau,
#             test = "knha")
# res.gosh <- gosh(m.rma)
# pdf(file = paste0(SNPout,"gosh.pdf"), width = 10, height = 13)
# plot(res.gosh, alpha = 0.01)
# dev.off()

sink() 


