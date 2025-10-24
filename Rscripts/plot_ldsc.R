###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(CMplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
options(bitmapType = "cairo")


###############################functions###################################
TwoTraitsManPlot = function(gwas_sub, name, collist){
        CMplot(gwas_sub, col=collist, plot.type="m",multraits=TRUE,
        threshold=c(1e-5,5e-8),threshold.lty=c(1,2), 
        threshold.lwd=c(1,1), threshold.col=c("black","black"), amplify=TRUE,bin.size=1e6,
        chr.den.col=NULL, signal.col=NULL,signal.cex=1, 
        file="jpg",file.name=name,dpi=300,file.output=TRUE,verbose=TRUE,
        points.alpha=100,legend.ncol=1, legend.pos="left", width=18,height=6,chr.labels.angle=45)

        CMplot(gwas_sub, plot.type="q",col=collist,multraits=TRUE,
        threshold=1e-6,ylab.pos=2,signal.pch=c(19,6),signal.cex=1.2,signal.col="#ff0000",
        conf.int=TRUE,box=TRUE,axis.cex=1,file="jpg",file.name=name,dpi=300,
        file.output=TRUE, verbose=TRUE,ylim=c(0,8),width=5,height=5)

        message("Successfully generated the Manhattan plots for two traits ", name)
        message("Successfully generated the QQ plots for two traits ", name)

}

OneTraitManPlot = function(gwas_sub, name, onecol){
        CMplot(gwas_sub, col=c(onecol,"grey60"), plot.type="m", LOG10=TRUE, ylim=NULL, 
        threshold=c(1e-5,5e-8),threshold.lty=c(1,2),
        threshold.lwd=c(1,1), threshold.col=c("black","black"), amplify=TRUE,bin.size=1e6,
        signal.col=c("#ff0000","#a1ff0a"),signal.cex=c(1.5,1.5),signal.pch=c(19,19),
        file="jpg",file.name=name,dpi=300,file.output=TRUE,verbose=TRUE,
        width=18,height=6,chr.labels.angle=45)

        CMplot(gwas_sub, plot.type="q",col=onecol,
        box=TRUE,threshold.col="red",threshold.lty=2,
        conf.int=TRUE,conf.int.col=NULL, 
        file="jpg",file.name=name,dpi=300,
        file.output=TRUE,verbose=TRUE,width=5,height=5)

        message("Successfully generated the Manhattan plots for one trait ", name)
        message("Successfully generated the QQ plots for one trait ", name)

}


#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
filename <- arguments[2]
phenotype <- arguments[3]

setwd(setpath)

if (phenotype != "gwas_smoking"){
        filename = paste0(phenotype,"SD.LocusZoom")
        METAresult = fread(filename, header = T, data.table=F, stringsAsFactors=F)
        METAresult = subset(METAresult, select = c("MarkerName","Chr","BP","P-value"))
        colnames(METAresult) = c("MarkerName","Chr","BP",paste0(phenotype,"SD"))
        row.names(METAresult) = METAresult$MarkerName

        filenamess = paste0(phenotype,"ssSD.LocusZoom")
        METAresultSS = fread(filenamess, header = T, data.table=F, stringsAsFactors=F)
        METAresultSS = subset(METAresultSS, select = c("MarkerName","Chr","BP","P-value"))
        colnames(METAresultSS) = c("MarkerName","Chr","BP",paste0(phenotype,"ssSD"))
        row.names(METAresultSS) = METAresultSS$MarkerName
        
        # Merge the two data frames by MarkerName
        Clocks = merge(METAresult, METAresultSS[,c("MarkerName",paste0(phenotype,"ssSD"))], by.x = "MarkerName", by.y="MarkerName", all = FALSE)
        # colnames(Clocks) = c("MarkerName", "Chr", "BP", paste0(phenotype,"SD"), paste0(phenotype,"ssSD"))
}else{
        METAresult = fread(filename, header = T, data.table=F, stringsAsFactors=F)
        METAresult = METAresult[,c("MarkerName","Chr","BP","P-value")]
        row.names(METAresult) = METAresult$MarkerName
}



if (phenotype != "gwas_smoking"){

        if (phenotype == "DNAmAge"){
        collist = c("#f18f01","#048ba8")
        } else if (phenotype == "PhenoAge"){
        collist = c("#57737a","#68edc6")
        } else if (phenotype == "DunedinPACE"){
        collist = c("#e8e288","#08a045")
        } else { 
        message("ERROR:not matching phenotype")
        }

        OneTraitManPlot(METAresult, paste0(phenotype,"SD"), collist[1])
        OneTraitManPlot(METAresultSS, paste0(phenotype,"ssSD"), collist[2])
        TwoTraitsManPlot(Clocks, phenotype, collist)

} else {
        OneTraitManPlot(METAresult, "gwas_smoking", "grey30")
}




library(readxl)
EpiAge_Traits <- read_excel("EpiAge_Traits.xlsx")
View(EpiAge_Traits)

library(reshape2)
melted_EpiAge_Traits <- melt(EpiAge_Traits)
colnames(melted_EpiAge_Traits)=c("EpiAgeAcc", "Traits", "Rg")

library(ggplot2)
ggplot(data = melted_EpiAge_Traits, aes(x=EpiAgeAcc, y=Traits, fill=Rg)) + 
  geom_tile()+
  geom_text(aes(label = round(Rg, 2))) +
  scale_fill_gradient2(
    high = 'dodgerblue4',
    mid = 'white',
    low = 'firebrick2',
    limits = c(-0.5, 0.5),
    midpoint = 0
  )+theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+labs(
    title = "Genetic Correlation among Tratis",
  )

######################################################

EpiAge_Traits_P <- read_excel("EpiAge_Traits_Pvalue.xlsx")


library(reshape2)
melted_EpiAge_Traits_P <- melt(EpiAge_Traits_P)
colnames(melted_EpiAge_Traits_P)=c("EpiAgeAcc", "Traits", "Pvalue")

library(ggplot2)
ggplot(data = melted_EpiAge_Traits_P, aes(x=EpiAgeAcc, y=Traits, fill=Pvalue)) + 
  geom_text(aes(label = Pvalue)) +geom_tile()+
  geom_text(aes(label = Pvalue)) +
  scale_fill_gradient(
    high = "white", low = "red"
  )+theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+labs(
    title = "P value for Genetic Correlation among Tratis",
  )

# 6 epi age####################################################################

library(readxl)
EpiAge_Traits_o <- read_excel("EpiAge_Own.xlsx")
EpiAge_Traits_o <- as.data.frame(EpiAge_Traits_o)
EpiAge_Traits_o$DNAmAgeSD = as.numeric(EpiAge_Traits_o$DNAmAgeSD)
EpiAge_Traits_o$DNAmAgessSD = as.numeric(EpiAge_Traits_o$DNAmAgessSD)
rownames(EpiAge_Traits_o) = EpiAge_Traits_o$...1
EpiAge_Traits_o=EpiAge_Traits_o[,2:7]
library(ggcorrplot)


EpiAge_p <- read_excel("EpiAge_Own_Pvalue.xlsx")
EpiAge_p <- as.data.frame(EpiAge_p)
rownames(EpiAge_p) = EpiAge_p$colname
EpiAge_p=EpiAge_p[,2:7]

ggcorrplot(EpiAge_Traits_o,
           ggtheme = ggplot2::theme_gray,
           outline.color = "white",
           type = "lower",
           colors = c("#6D9EC1", "white", "#E46726"),
           lab=T,
           title = "Genetic correlations for the epigenetic age acc")

ggcorrplot(EpiAge_p,
           ggtheme = ggplot2::theme_gray,
           outline.color = "white",
           type = "lower",
           colors = c("white","red"),
           title = "P value for the epigenetic age acc")

######################################################

EpiAge_Traits_P <- read_excel("EpiAge_Traits_Pvalue.xlsx")


library(reshape2)
melted_EpiAge_Traits_P <- melt(EpiAge_Traits_P)
colnames(melted_EpiAge_Traits_P)=c("EpiAgeAcc", "Traits", "Pvalue")

library(ggplot2)
ggplot(data = melted_EpiAge_Traits_P, aes(x=EpiAgeAcc, y=Traits, fill=Pvalue)) + 
  geom_text(aes(label = Pvalue)) +geom_tile()+
  geom_text(aes(label = Pvalue)) +
  scale_fill_gradient(
    high = "white", low = "red"
  )+theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+labs(
    title = "P value for Genetic Correlation among Tratis",
  )

############################################