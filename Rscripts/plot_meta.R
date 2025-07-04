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

