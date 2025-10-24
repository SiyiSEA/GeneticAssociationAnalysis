###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hexbin))
options(bitmapType = "cairo")


#################################functions####################################

Visluazaion <- function(meta, plotname){
    # Z-score vs P-value
    ZPplot = ggplot(meta, aes(x = Zscore, y = -log10(`P-value`))) +
                geom_bin2d(bins = 100) +
                scale_fill_viridis_c(trans = "log10", option = "D") +
                geom_smooth(method = "loess", color = "white", se = FALSE, linewidth = 0.8) +
                labs(x    = "Z-score",
                    y    = expression(-log[10]*"(P-value)"),
                    fill = "Count") +
                theme_minimal(base_size = 11)

    # Effect vs StdErr, colored by P-value_SE
    ESPplot = ggplot(meta, aes(x = Effect, y = StdErr, color = `P-value_SE`)) +
                geom_point(size = 0.5) +
                scale_color_gradient(low = "steelblue", high = "firebrick") +
                labs(x     = "Effect Size",
                    y     = "Standard Error",
                    color = "SE of P-value") +
                theme_minimal(base_size = 11) 

    # Heterogeneity metrics: HetISq, HetChiSq & HetPVal
    HHHplot = ggplot(meta, aes(x = HetISq, y = HetChiSq, color = -log10(HetPVal))) +
                geom_point(size = 0.5) +
                scale_color_viridis_c() +
                labs(x = expression(I^2), 
                    y = "Chi-square", 
                    color = expression(-log[10]*"(HetPVal)")) +
                theme_minimal(base_size = 10)

    # FreqSE vs P-value vs P-value_SE
    FPPplot = ggplot(meta, aes(x = -log10(`P-value_SE`), y = -log10(`P-value`), color = `FreqSE` )) +
                geom_point(size = 0.5) +
                scale_y_continuous(expand = c(0,0)) +
                scale_color_viridis_c(option = "magma") +
                labs(
                    x = paste(expression(-log[10]*"(P-value_SE)"),'Standard Error Scheme'),
                    y = paste(expression(-log[10]*"(P-value)"),'Sample Size Scheme'),
                    color = "FreqSE") +
                theme_minimal(base_size = 11)

  
    plotlist = list(ESPplot, HHHplot, FPPplot)
    namelist = c("ESP", "HHH", "FPP")
    for (Nplot in 1:length(plotlist)) {
        outname = paste0("./plots/",plotname, "_", namelist[Nplot], ".tiff")
        message("Saving plot: ", outname)
        ggsave(plot = plotlist[[Nplot]], filename = outname, width = 8, 
                height = 6, dpi = 400)
    }
    

}



#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
filename <- arguments[2]
phenotype <- arguments[3]

setwd(setpath)
METAresult = fread(filename, header = T, data.table=F, stringsAsFactors=F)

# SE > 0.2 will be regrads as large difference, but in this case, I will remove the SNPs with SE > 0.3
message("How many SNPs with FreqSE in the sample size scheme > 0.3?")
dim(METAresult[METAresult$FreqSE > 0.3,])[1]
FreqSE_SS = dim(METAresult[METAresult$FreqSE > 0.3,])[1]
message("How many SNPs with FreqSE_SE in the standard error scheme > 0.3?")
dim(METAresult[METAresult$FreqSE_SE > 0.3,])[1]
FreqSE_SE = dim(METAresult[METAresult$FreqSE_SE > 0.3,])[1]

# if (FreqSE_SS > FreqSE_SE){
#     message("FreqSE_SS > FreqSE_SE, so I will use FreqSE to filter out the SNPs")
#     METAresult_filter <- METAresult[METAresult$FreqSE < 0.3,]
# }else{
#     message("FreqSE_SE > FreqSE_SS, so I will use FreqSE_SE to filter out the SNPs")
#     METAresult_filter <- METAresult[METAresult$FreqSE_SE < 0.3,]
# }

# HetIsq can detact the heterogeneity better than QC, so I will use HetIsq to filter out the SNPs
message("How many SNPs with HetIsq > 80 and with HetPval < 0.05?")
dim(METAresult[which(METAresult$HetISq > 80 & METAresult$HetPVal < 0.05),])[1]
message("Filter out SNPs with HetIsq > 80 and with HetPval < 0.05")
METAresult <- METAresult[METAresult$HetISq < 80 & METAresult$HetPVal > 0.05,]

# Visluazaion
# Visluazaion(METAresult, phenotype)

# message("Saving the filtered results")
# fwrite(METAresult, file = paste0(phenotype,".filtered.tbl"), quote = F, sep = "\t", col.names = T, row.name = F)

# Make a list of SNPs with p-value < 1e-5
message("Saving the SNPs with p-value < 1e-5")
SNPlist <- METAresult[METAresult$`P-value` < 5e-8,]
write.table(SNPlist, file = paste0(phenotype , "_sigSNPlist.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

# Formatted the filter results for upload into COJO
message("Saving the SNPs for COJO")
METACOJO <- METAresult[,c("MarkerName", "Allele1", "Allele2","Freq1","Effect","StdErr","P-value_SE","N","Weight","Zscore","P-value")]
fwrite(METACOJO, file = paste0(phenotype,"_input.cojo"), quote = F, sep = "\t", col.names = T, row.name = F )

# Formatted the filter results for upload into locuszoom
METACOJO[c('Chr', 'BP_ALT_REF')] <- str_split_fixed(METACOJO$MarkerName, ':', 2)
METACOJO[c('BP', 'ALT_REF')] <- str_split_fixed(METACOJO$BP_ALT_REF, '_', 2)
METACOJO[grep('X', METACOJO[,c("Chr")], ignore.case = TRUE),c("Chr")] = as.numeric(23)
message("Saving the SNPs for LocusZoom")
METALoczoom = METACOJO[,c("Chr","BP","MarkerName", "Allele1", "Allele2","Freq1","Effect","StdErr","P-value_SE","N","Weight","Zscore","P-value")]
METALoczoom$Chr = as.numeric(METALoczoom$Chr)
METALoczoom$BP = as.numeric(METALoczoom$BP)
METALoczoom$`P-value` = as.numeric(METALoczoom$`P-value`)
METALoczoom = dplyr::arrange(METALoczoom, Chr, BP)
fwrite(METALoczoom, file = paste0(phenotype,".LocusZoom"), quote = F, sep = "\t", col.names = T, row.name = F )
METACOJO[grep('23', METACOJO[,c("Chr")], ignore.case = TRUE),c("Chr")] = as.character("X")
fwrite(METALoczoom, file = paste0(phenotype,".LocusZoom.X"), quote = F, sep = "\t", col.names = T, row.name = F )

# Formatted the filter results for upload into LDSC
message("Saving the SNPs for LDSC")
LDSCout = METAresult[,c("MarkerName", "Allele1", "Allele2", "Freq1", "FreqSE","MinFreq","MaxFreq","N", "Weight", "Zscore", "P-value", "Direction")]
LDSCout[c('Chr', 'BP_ALT_REF')] <- str_split_fixed(LDSCout$MarkerName, ':', 2)
LDSCout[c('BP', 'ALT_REF')] <- str_split_fixed(LDSCout$BP_ALT_REF, '_', 2)
LDSCout = LDSCout[,c("Chr","BP","MarkerName", "Allele1", "Allele2","Freq1","FreqSE","MinFreq","MaxFreq","N","Weight","Zscore","P-value", "Direction")]
LDSCout$BP = as.numeric(LDSCout$BP)
LDSCout$`P-value` = as.numeric(LDSCout$`P-value`)
LDSCout = dplyr::arrange(LDSCout, Chr, BP)
fwrite(LDSCout, file = paste0(phenotype,".LDSCinput"), quote = F, sep = "\t", col.names = T, row.name = F )