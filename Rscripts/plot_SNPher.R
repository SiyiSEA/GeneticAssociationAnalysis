###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
#display.brewer.all()
options(bitmapType = "cairo")


###############################functions###################################



#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
filename <- arguments[2]
phenotype <- arguments[3]

setwd(setpath)

if (phenotype == "All") {
        message("Ploting the SNP heritbility for meta results across all 7 phenotypes")
        her_table = read.table(filename, header=T, sep = " ")
        her_table$Phenotype <- factor(her_table$Phenotype, levels = her_table$Phenotype)
        p = ggplot(her_table, aes(x = Phenotype, y = SNPHeritability)) +
                geom_col(fill = c("#f18f01","#048ba8","#57737a","#68edc6","#e8e288","#08a045","grey30")) +  # bar plot
                geom_errorbar(aes(ymin = SNPHeritability - SE, ymax = SNPHeritability + SE), 
                                width = 0.2, color = "black") +  # error bars
                theme_light() +
                labs(
                title = "SNP Heritability with SE across 7 phenotypes",
                x = "Phenotypes",
                y = "SNP Heritability") + theme(
                axis.text.x = element_text(angle = 45, hjust = 1),  # tilt labels
                axis.title.x = element_text(margin = margin(t = 10)) # add space
                )+ scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
        ggsave("OverAll_her.png",  plot = p, width = 6, height = 6, dpi = 300)
}
