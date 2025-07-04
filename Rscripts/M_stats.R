###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(getmstatistic))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
options(bitmapType = "cairo")

#################################functions####################################
forestfunnalplot = function (forestframe, markername, markerP, type, standererror, phenotype){
  chr = str_split_fixed(markername,":",2)[1]
  bp = str_split_fixed(str_split_fixed(markername,"_",2)[1],":",2)[2]
  SNPname = paste0("SNP Chr", chr, " ", bp)
  SNPout = paste0("Chr", chr, "_", bp)
  message("Making forest plot and funnal plot for ", SNPname)
  plotpath = paste0("./plots/",phenotype, "_", SNPout, "_", type)
  
  # Compute inverse-variance weighted fixed effects model
  # using fixed effects model to get variable box size
  metafor_results_fe <- metafor::rma.uni(yi = forestframe[, type], sei = forestframe[, standererror],
                                          weighted = T, slab = paste(forestframe$Study, forestframe$N, sep = ","), 
                                          method = "FE")
  moduleinfo = paste0(
      "RE Model (Q = ", round(metafor_results_fe$QE, 2),
      ", df = ", metafor_results_fe$k - metafor_results_fe$p, ", p < ",
      round(metafor_results_fe$QEp,3), "; ",
      "I^2 = ", round(metafor_results_fe$I2,2), "%, ",
      "tau^2 = ", round(metafor_results_fe$tau2,2), ")")

  pdf(file=paste(plotpath,'forest.pdf', sep="_"))
  par(mar=c(4,2,2,2))
  xtitle=paste0(type, "statistic based on ",SNPname," with P value ", markerP)
  forest(metafor_results_fe, xlim=c(-1.8, 1.6), at=c(-1, -0.5, 0, 0.5, 1), cex=0.66, 
          xlab = xtitle, order = "obs", ilab.xpos = c(-1.1), ilab.pos = c(2), addfit=F)
  mtext(paste(type,"forest plot for",phenotype), side = 3, line = 1, cex = 0.9, font = 2)
  mtext(moduleinfo, side = 3, cex = 0.7, font = 2)
  dev.off()
  
  pdf(file=paste(plotpath,'funnel.pdf', sep="_"))
  funnel(metafor_results_fe, xlab = paste0(type,"statistic based on ",SNPname),
        slab = metafor_results_fe$slab,label="all",legend = T)
        # label = 3, the most extreme, and the second and third most extreme points are labeled
  dev.off()

  message(paste(plotpath,'forest.pdf', sep="_")," has been generated successfully.")
  message(paste(plotpath,'funnel.pdf', sep="_")," has been generated successfully.")
}

#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
Mstat <- arguments[2]
effectsize <- arguments[3]
lambda <- arguments[4]
phenotype <- arguments[5]
SigSNPs <- arguments[6:length(arguments)]

setpath="/lustre/home/sww208/GoDMC/GADatasets/CorrectGWAS10/Mstatresults/Mstat_DNAmAgeSD"
Mstat="SignalSNPs.pre"
effectsize="effectsize_MstatSignal.pre"
lambda="/lustre/home/sww208/GoDMC/GADatasets/CorrectGWAS10/Resources/cohort_lambda.tsv"
phenotype="DNAmAgeSD"
SigSNPs="6:18141639_C_T 6:18104276_A_C 1:236518226_C_T"

setwd(setpath)
message("Reading in the files for Mstat")
Mresults = fread(Mstat, header = T, data.table=F, stringsAsFactors=F)
Effectsize = fread(effectsize, header = F, data.table=F, stringsAsFactors=F)
colnames(Effectsize) = c("SNP", "Effect", "StdErr","P")
Lambda = fread(lambda, header = T, data.table=F, stringsAsFactors=F)

if ( phenotype == "gwas_smoking") {
  lambdacol = "Lambda_Smoke"
} else {
  lambdacol = paste0("Lambda_", phenotype)
}


message("How many studies are in the Mstat file?")
message(length(unique(Mresults$Study)))
message("How many studies overlapped in Mstat file and lambda file?")
length(intersect(Mresults$Study, Lambda$Study))
message("How many studies are in Mstat file but not in the lambda file?")
setdiff(Mresults$Study,Lambda$Study)
message("How many studies are in lambda file but not in the Mstat file?")
setdiff(Lambda$Study,Mresults$Study)

venn.diagram(
  x = list(Mresults$Study, Lambda$Study),
  category.names = c("Mstat Studies" , "All Studies"),
  fill = c("#f6cfbe", "#b9dcf2"), cat.default.pos = "text",
  margin = 0.1,
  filename = paste0("./plots/",phenotype,"_CohortsVenn.png"),
  output=T)

# prepare the table for Mstat
MstatTable = merge(Mresults, Effectsize, by.x = "SNP", by.y = "SNP")
MstatTable = merge(MstatTable, Lambda, by.x = "Study", by.y = "Study")
MstatTable$gcse = MstatTable$SE * sqrt(MstatTable[,lambdacol])

getmstatistic_results = getmstatistic(as.numeric(MstatTable$BETA), 
                                      as.numeric(MstatTable$gcse), 
                                      MstatTable$SNP, 
                                      MstatTable$Study,
                                      save_dir = paste0(setpath,"/plots/"),)

print("Retrieve dataset of stronger than average studies (significant at 5% level)")
print(getmstatistic_results$influential_studies_0_05)

print("Retrieve dataset of weaker than average studies (significant at 5% level)")
print(getmstatistic_results$weaker_studies_0_05)

print("Retrieve number of studies and variants")
print(getmstatistic_results$number_studies)
print(getmstatistic_results$number_variants)

print("Retrieve expected mean, sd and critical M value at 5% significance level")
print(getmstatistic_results$M_expected_mean)
print(getmstatistic_results$M_expected_sd)
print(getmstatistic_results$M_crit_alpha_0_05)

print(paste("Saving M statistic results for ", phenotype))
save(getmstatistic_results, file = paste0(phenotype,"_MstatResults.RData"))
# load(paste0(phenotype,"_MstatResults.RData"))

#### plot-1
dframe = getmstatistic_results$M_dataset
dframe$study_names_in_n<-paste0(dframe$study_names_in," (n=",dframe$beta_n,"/",getmstatistic_results$number_variants,")")
# does beta_n is the sample size of the study? - No, beta_n is the number of variants in the study

drame_effectsize = merge(dframe[,c("study_names_in","variant_names_in","beta_in","study_names_in_n")], Effectsize, 
                                  by.x = "variant_names_in", by.y = "SNP")
drame_effectsize$labelA = "After Meta-analysis"
drame_effectsize$labelB = "Before Meta-analysis"
boxplot_effectsizeA = drame_effectsize[,c("study_names_in_n","Effect", "labelA")]
colnames(boxplot_effectsizeA) = c("study","value", "label")
boxplot_effectsizeB = drame_effectsize[,c("study_names_in_n","beta_in", "labelB")]
colnames(boxplot_effectsizeB) = c("study","value", "label")
boxplot_effectsizeC=rbind(boxplot_effectsizeA, boxplot_effectsizeB)
boxplot_effectsizeC$median = ave(boxplot_effectsizeC$value, boxplot_effectsizeC$study, boxplot_effectsizeC$label, FUN = median)

boxplot_effectsizeC$median <- ave(
  boxplot_effectsizeC$value,
  boxplot_effectsizeC$study,
  boxplot_effectsizeC$label,
  FUN = function(x) median(x, na.rm = TRUE)
)

p1 = ggplot(boxplot_effectsizeC, aes(x=study, y=value, fill=label)) + 
    geom_boxplot()+ theme_bw() + scale_fill_manual(values=c("#f8c7cc", "#81a684")) +
    labs(x = "Study", y = "Effect Size", title = "Boxplot of effect size across all studies") +
    theme(axis.text.x=element_text(angle=45,hjust=1))
ggsave(plot=p1, file=paste0("./plots/",phenotype,"_effectsize.pdf"), width=14, height=7, dpi = 100, limitsize=F)
message(paste0("./plots/",phenotype,"_effectsize.pdf has been generated successfully."))

### plot-2 M statistic vs avergae effect size
p2=ggplot(dframe, aes(x = oddsratio, y = M, color = M)) +
  geom_point(size = 3) +
  scale_color_gradientn(name = "M value", colors = rainbow(7)) +
  labs(x     = "Effect size",
       y     = "M statistic",
      title = paste("Effect size vs. M across studies for ", phenotype)) +
  theme_bw() +
  theme( plot.title  = element_text(size = 14, face = "bold"),
       axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_text( label=dframe$study_names_in,
              hjust = 0, nudge_x=0.002, size=2,
              check_overlap=T)
ggsave(p2,width=10, height=6, dpi = 100,file=paste0("./plots/",phenotype,"_cohortMstats.pdf"))
message(paste0("./plots/",phenotype,"_cohortMstats.pdf has been generated successfully."))

### plot-3+4 M forest plot + beta/effect size forest plot
getm <- dplyr::arrange(getmstatistic_results$M_dataset, M)

for (OneSigSNP in as.list(strsplit(SigSNPs,split=" "))[[1]]){
  print(paste0("Making M forest plot for ", OneSigSNP))
  markername = as.character(OneSigSNP)
  markerP = Effectsize[Effectsize$SNP == markername,"P"]
  MstatSNP = MstatTable[MstatTable$SNP == markername,]
  getmSNP = subset(getm, getm$variant_names_in == markername)
  Mforest = merge(MstatSNP[, c("Study", "N")],getmSNP,by.x="Study", by.y="study_names_in")
  Mforest = Mforest[!duplicated(Mforest), ]
  Mforest = dplyr::arrange(Mforest, M)
  forestfunnalplot(Mforest, markername, markerP, "M", "M_se", phenotype)
  BETAforest = MstatSNP[!duplicated(MstatSNP), ]
  BETAforest = dplyr::arrange(BETAforest, BETA)
  forestfunnalplot(BETAforest, markername, markerP, "BETA", "SE", phenotype)

}



