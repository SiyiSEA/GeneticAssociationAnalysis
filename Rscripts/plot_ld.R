###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))

options(bitmapType = "cairo")


#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
filename <- arguments[2]
phenotype <- arguments[3]

setwd(setpath)
LDresult = fread(filename, header = T, data.table=F, stringsAsFactors=F)

row.names(LDresult) = LDresult$SNP
LDresult = subset(LDresult, select = -c(SNP))

# check if the LDresult is symmetric
while ( ncol(LDresult) != nrow(LDresult) ) {
        message("The LD result is not symmetric, dealding with it...")
        LDresult = subset(LDresult, select = rownames(LDresult))
}

message("The LD result is symmetric, plotting it...")
pdf(file=paste0(phenotype, "_LD.pdf"), 
     width = 7, height = 6)
heatmap(as.matrix(LDresult), scale="column", col=heat.colors(5), margins = c(8, 13),
main = paste0(phenotype,": LD Correlation of Signal SNPs"))
legend(x="right", legend=c(0, 0.25,0.5, 0.75,1),fill=heat.colors(5), cex = 0.8)
dev.off()
