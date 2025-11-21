###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))
#display.brewer.all()
options(bitmapType = "cairo")


###############################functions###################################



#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
filename <- arguments[2]

setwd(setpath)

lines <- readLines(filename)

# Extract blocks
starts <- grep("^.* and .*----", lines)

extract_block <- function(i) {
  trait_line <- lines[starts[i]]
  parts <- strsplit(trait_line, " and |----")[[1]]
  trait1 <- parts[1]
  trait2 <- parts[2]

  gc_line <- lines[starts[i] + 1]
  z_line  <- lines[starts[i] + 2]
  p_line  <- lines[starts[i] + 3]

  # Extract GC and SE
  gc_num <- sub(".*Genetic Correlation: ([0-9.eE+-]+) \\(.*", "\\1", gc_line)
  gc_se  <- sub(".*\\(([0-9.eE+-]+)\\).*", "\\1", gc_line)

  # Z
  z <- sub("Z-score: ", "", z_line)

  # P-value: remove spaces and handle e-02, e-2, etc
  p <- sub("P: ", "", p_line)
  p <- gsub(" ", "", p)     # remove any stray spaces
  p <- gsub("e-0", "e-", p) # convert e-02 → e-2 (R accepts both)
  
  tibble(
    Trait1 = trait1,
    Trait2 = trait2,
    GeneticCorrelation = as.numeric(gc_num),
    SE = as.numeric(gc_se),
    Z = as.numeric(z),
    P = as.numeric(p)
  )
}

df <- bind_rows(lapply(seq_along(starts), extract_block))

message("Genetic Correlation Results among clocls :")
print(df)


message("Making the heatmap plots among clocks ...")
traits <- unique(c(df$Trait1, df$Trait2))
df <- df %>%
  mutate(Trait1 = factor(Trait1, levels = traits),
         Trait2 = factor(Trait2, levels = traits))

df <- df %>%
  mutate(label = paste0(round(GeneticCorrelation, 3), 
                        " (", round(SE, 3), ")\nP=", signif(P, 3)))

# ---- Plot heatmap ----
p = ggplot(df, aes(x = Trait1, y = Trait2, fill = GeneticCorrelation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0.5, limits = c(0,1), name = "GeneticCorr") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave("OverAllClocks_GeneCorr.png",  plot = p, width = 8, height = 8, dpi = 300)
message("Heatmap plot saved as OverAllClocks_GeneCorr.png")