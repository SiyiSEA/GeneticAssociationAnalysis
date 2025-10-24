###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(meta))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
options(bitmapType = "cairo")

#################################functions####################################
check_overlap <- function(cohort1, cohort2, tobeCheck) {
    if (length(setdiff(cohort1,  cohort2)) > 0 | length(setdiff(cohort1, cohort2)) > 0) {
        message(paste0("Checking the overlapped cohorts for ", tobeCheck, "---------"))
        message("Studies not in the meta data, please check the following studies:")
        print(setdiff(cohort1, cohort2))
        message("Studies not in the checked data, please check the following studies:")
        print(setdiff(cohort2, cohort1))
    }
}

#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
Mstat <- arguments[2]
metaobject <- arguments[3]
age <- arguments[4]
array <- arguments[5]
imputation <- arguments[6]
phenotype <- arguments[7]

Mstat="/lustre/home/sww208/GoDMC/GeneticAssociationAnalysis/Mstatresults/Mstat_DNAmAgeSD/DNAmAgeSD_MstatResults.RData"
metaobject = "/lustre/home/sww208/GoDMC/GeneticAssociationAnalysis/Mstatresults/Mstat_DNAmAgeSD/DNAmAgeSD-Chr6_18130918.RData"
age = "/lustre/home/sww208/GoDMC/GeneticAssociationAnalysis/Resources/cohort_age.tsv"
array = "/lustre/home/sww208/GoDMC/GeneticAssociationAnalysis/Resources/cohort_array.tsv"
impuation = "/lustre/home/sww208/GoDMC/GeneticAssociationAnalysis/Resources/cohort_imputedpanel.tsv"
phenotype = "DNAmAgeSD"

setwd(setpath)
message("Reading in the files")
load(Mstat)
load(metaobject)
age = fread(age, header = T, data.table=F, stringsAsFactors=F)
array = fread(array, header = T, data.table=F, stringsAsFactors=F)
impuation = fread(impuation, header = T, data.table=F, stringsAsFactors=F)
SNPname = str_split_fixed(tools::file_path_sans_ext(basename(metaobject)), "-", 2)[,2]

#### [meta-analysis of perticular significantly SNP]
print("==========================================================================")
print(paste("*******************Based on the Meta result of one ", SNPname, "Significant SNP********************"))
print("==========================================================================")
#### Meta-regression for stratificated Age
agemgdata = mg$data
mgdata$Study = str_split_fixed(mgdata$.studlab , " ", 2)[,1]
mgdf = data.frame(mg = c(summary(mg)$tau2, summary(mg)$I2 * 100, summary(mg)$Q, summary(mg)$pval.Q))
rownames(mgdf) = c("tau2", "I2", "Q", "pval.Q")
# Split cohorts by the gap of MinAge, MaxAge, MeanAge, and MedianAge
age = age %>% 
    mutate(GroupbyMin = ifelse(MinAge < 22, "young", "adults")) %>%
    mutate(GroupbyMax = ifelse(MaxAge < 30, "young", "adults")) %>%
    mutate(GroupbyMean = ifelse(MeanAge < 32, "young", "adults")) %>%
    mutate(GroupbyMedian = ifelse(MedianAge < 39, "young", "adults"))

age_sort = subset(age, age$Studyname %in% mgdata$Study)
mgdata_merged = merge(mgdata, age_sort, by.x = "Study", by.y = "Studyname", all.x = TRUE)


pdf(file=paste0("AgeBubble_",tools::file_path_sans_ext(basename(metaobject)) ,".pdf"), width=8, height=12)
par(mfrow=c(3,2))
for (covar in c("MeanAge", "MedianAge", "GroupbyMin", "GroupbyMax", "GroupbyMean", "GroupbyMedian")) {
    if (is.null(mgdata_merged[[covar]])) {
        stop(paste0("The covariate ", covar, " is not found in the merged data. Please check the input files."))
    }else {
       covarEnv = mgdata_merged[, covar]
       m.gen.reg <- metareg(mg, ~covarEnv)
       mgdf[,c(covar)] = c(summary(m.gen.reg)$tau2, summary(m.gen.reg)$I2, summary(m.gen.reg)$QE, summary(m.gen.reg)$pval[1])
       print(paste0("Meta-regression results for ", covar, " for All cohorts=========================="))
       print(summary(m.gen.reg))
       bubble(m.gen.reg, studlab = TRUE, xlab = covar, ylab = "Standardized Mean Difference (SMD)", col.line = "red")
    }
}
dev.off()

#### Meta-regression for Array
check_overlap(array$Studyname, mgdata$Study, "Array")
array_sort = subset(array, array$Studyname %in% mgdata$Study)
mgdata_merged = merge(mgdata, array_sort, by.x = "Study", by.y = "Studyname", all.x = TRUE)
coverArray = mgdata_merged$Array
m.gen.reg <- metareg(mg, ~coverArray)
message("Meta-regression results for Array for All cohorts==========================")
print(summary(m.gen.reg))
mgdf[,c("Array")] = c(summary(m.gen.reg)$tau2, summary(m.gen.reg)$I2, summary(m.gen.reg)$QE, summary(m.gen.reg)$pval[1])
jpeg(file=paste0("ArrayBubble_",tools::file_path_sans_ext(basename(metaobject)) ,".pdf"))
bubble(m.gen.reg, studlab = TRUE, xlab = "Array Type", ylab = "Standardized Mean Difference (SMD)",
col.line = "red")
dev.off()


#### Meta-regression for Imputation Panel
check_overlap(impuation$Studyname, mgdata$Study, "ImputationPanel")
impuation_sort = subset(impuation, impuation$Studyname %in% mgdata$Study)
mgdata_merged = merge(mgdata, impuation_sort, by.x = "Study", by.y = "Studyname", all.x = TRUE)
coverImputation = mgdata_merged$ImputationPanel
m.gen.reg <- metareg(mg, ~coverImputation)
message("Meta-regression results for Imputation panels for All cohorts==========================")
print(summary(m.gen.reg))
mgdf[,c("Imputation")] = c(summary(m.gen.reg)$tau2, summary(m.gen.reg)$I2, summary(m.gen.reg)$QE, summary(m.gen.reg)$pval[1])
jpeg(file=paste0("ImputationBubble_",tools::file_path_sans_ext(basename(metaobject)) ,".pdf"))
bubble(m.gen.reg, studlab = TRUE, xlab = "Imputation Panels", ylab = "Standardized Mean Difference (SMD)",
col.line = "red")
dev.off()



# todo
# 还咩有完成的是MAF和INFOscore，这两是要从brisotl里面重新抓
# case and control的也没有完成
# number of CpG sites
# males

#### [meta-analysis of M-statistic across all significant SNPs]
print("==========================================================================")
print(paste("*******************Based on the M-stat of one ", SNPname, "Significant SNP********************"))
print("==========================================================================")
#### Meta-regression for stratificated Age
Mtable = getmstatistic_results$M_dataset
SigSNP = gsub("_", ":", SNPname)
SigSNP = gsub("Chr", "^", SigSNP)

# 我发现一个问题：
# SigSNP 是我后加的，很有可能不存在SignalCojoSNP-03里面
# 不能用同样的SNP来跑Mstatistic的meta-regression
# 解决方案：
# 1. 通过LD来找到和SigSNP相关的SNP在Mtable有的；
# 2. 直接从Mtable里面找一个pval最小的来跑；
# 3. 不管怎么样，要重新走一个script

Mtable %>% filter(variant_names_in %like% "%a%")

Mtable = Mtable[grepl(SigSNP, Mtable$variant_names_in, perl = TRUE),] 
grepl("^6:18130918", Mtable$variant_names_in, perl = TRUE)
grepl("^Si", c("SEA", "Siyi", "Wang"), perl = TRUE)
grepl("^6:18130918", unique(Mtable$variant_names_in), perl = TRUE)

age_sort = subset(age, age$Studyname %in% Mtable$study_names_in)
Mtable_merged = merge(Mtable, age_sort, by.x = "study_names_in", by.y = "Studyname", all.x = TRUE)

mrma.reg = rma.uni(yi=M, vi=0, sei=M_se, mods = ~ as.numeric(MedianAge), data = Mtable_merged)
message("Meta-regression results for MedianAge for All cohorts==========================")
print(rma.obj)

rma.obj = rma.uni(yi=M, vi=0, sei=M_se, mods = ~ as.numeric(MeanAge), data = MAgetable)
message("Meta-regression results for MeanAge for All cohorts:==========================")
print(rma.obj)

# subset cohorts with age
below20cohort = which(as.numeric(MAgetable$MedianAge) <= 20)
over20cohort = which(as.numeric(MAgetable$MedianAge) > 20)

rma.obj = rma.uni(yi=M, vi=0, sei=M_se, mods = ~ as.numeric(MedianAge), data = MAgetable[below20cohort,])
message("Meta-regression results for MedianAge for cohorts with Median Age <= 20==========================")
print(rma.obj)
rma.obj = rma.uni(yi=M, vi=0, sei=M_se, mods = ~ as.numeric(MeanAge), data = MAgetable[below20cohort,])
message("Meta-regression results for MeanAge for cohorts with Median Age <= 20==========================")
print(rma.obj)

rma.obj = rma.uni(yi=M, vi=0, sei=M_se, mods = ~ as.numeric(MedianAge), data = MAgetable[over20cohort,])
message("Meta-regression results for MedianAge for cohorts with Median Age > 20==========================")
print(rma.obj)
rma.obj = rma.uni(yi=M, vi=0, sei=M_se, mods = ~ as.numeric(MeanAge), data = MAgetable[over20cohort,])
message("Meta-regression results for MeanAge for cohorts with Median Age > 20==========================")
print(rma.obj)

# plot the relation between age and M statistic
p1 = ggplot(MAgetable,aes(x=MedianAge,y=M))+
    geom_point(aes(colour=factor(study_names_in), size =  MeanAge)) +
    xlab("MedianAge") +
    ylab("M Statistic") +
    theme(axis.text.x = element_text(face = "bold"))+ 
    labs(title=paste0("All cohorts - M statistic vs Median Age for ", TopSNPs$SNP[2]))
ggsave(p1,file=paste0(SNPout,"M_MedianAge_All.png"),width=12,height=8)


p2=ggplot(MAgetable[over20cohort,],aes(x=MedianAge,y=M))+
    geom_point(aes(colour=factor(study_names_in), size =  MeanAge)) +
    xlab("MedianAge") +
    ylab("M Statistic") +
    theme(axis.text.x = element_text(face = "bold"))+ 
    labs(title=paste0("Cohorts with Median Age over 20 - M statistic vs Median Age for ", TopSNPs$SNP[2]))
ggsave(p2,file=paste0(SNPout,"M_MedianAge_over20.png"),width=12,height=8)

