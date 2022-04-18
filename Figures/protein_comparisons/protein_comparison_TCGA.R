library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(ggprism)

### Set variables
options(scipen = 100)
path <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Protein_comparison_TCGA"
proteome <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Proteome/output/protein_abundance.txt"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Proteome/output/metadata.txt"

### Read in data
proteome <- read.delim(proteome)
samples <- read.delim(samples)

### Subset samples
samples <- samples[samples$Type %in% c("Oligodendroglioma", "Glioblastoma", "Astrocytoma", "Normal"), ]
samples$Type <- factor(samples$Type, levels = c("Normal", "Oligodendroglioma", "Astrocytoma", "Glioblastoma"),
                       labels = c("Normal", "Type I", "Type II", "Type III"))

### Format proteome data
proteome <- proteome[complete.cases(proteome), ]
proteome <- proteome[!(duplicated(proteome$Gene)), ]
row.names(proteome) <- proteome$Gene
proteome <- proteome[, c(10:ncol(proteome))]
proteome <- as.data.frame(t(proteome))

### Center and normalize data
scaled_data <- as.data.frame(scale(proteome))

### Merge with metadata
scaled_data <- merge(scaled_data, samples, by.x = "row.names", by.y = c("set"))

### Calculate p-values for each subtype
group1 <- c("Normal", "Normal", "Normal", "Type I", "Type I", "Type II")
group2 <- c("Type I", "Type II", "Type III", "Type II" ,"Type III" ,"Type III")

### SYK
a <- t.test(scaled_data$SYK[scaled_data$Type=="Normal"], scaled_data$SYK[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$SYK[scaled_data$Type=="Normal"], scaled_data$SYK[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$SYK[scaled_data$Type=="Normal"], scaled_data$SYK[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$SYK[scaled_data$Type=="Type I"], scaled_data$SYK[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$SYK[scaled_data$Type=="Type I"], scaled_data$SYK[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$SYK[scaled_data$Type=="Type II"], scaled_data$SYK[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
SYK <- as.data.frame(cbind(group1, group2, t_test))
SYK$p <- ifelse(SYK$t_test < 0.05 & SYK$t_test > 0.01, "*",
                 ifelse(SYK$t_test < 0.01 & SYK$t_test > 0.001, "**",
                        ifelse(SYK$t_test < 0.001, "***", "")))
SYK <- SYK[SYK$p %in% c("*", "**", "***"), colnames(SYK) %in% c("group1", "group2", "p")]
SYK$y.position <- c(2.5)

### CDH1 not fully detected

### ANXA1
a <- t.test(scaled_data$ANXA1[scaled_data$Type=="Normal"], scaled_data$ANXA1[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$ANXA1[scaled_data$Type=="Normal"], scaled_data$ANXA1[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$ANXA1[scaled_data$Type=="Normal"], scaled_data$ANXA1[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$ANXA1[scaled_data$Type=="Type I"], scaled_data$ANXA1[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$ANXA1[scaled_data$Type=="Type I"], scaled_data$ANXA1[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$ANXA1[scaled_data$Type=="Type II"], scaled_data$ANXA1[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
ANXA1 <- as.data.frame(cbind(group1, group2, t_test))
ANXA1$p <- ifelse(ANXA1$t_test < 0.05 & ANXA1$t_test > 0.01, "*",
                  ifelse(ANXA1$t_test < 0.01 & ANXA1$t_test > 0.001, "**",
                         ifelse(ANXA1$t_test < 0.001, "***", "")))
ANXA1 <- ANXA1[ANXA1$p %in% c("*", "**", "***"), colnames(ANXA1) %in% c("group1", "group2", "p")]
ANXA1$y.position <- c(1, 1.5, 3, 2, 3.5)

### HER2/3 not detected

### Set Theme
my_theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 20), 
                  axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  legend.key = element_rect(fill = "white"),
                  legend.text = element_text(size = 15),
                  legend.position = "none",
                  axis.text = element_text(size = 15),
                  axis.title = element_text(size = 15),
                  axis.text.x = element_text(angle = 45, hjust = 1))

### SYK
SYK_plot <- ggplot(scaled_data, aes(Type, SYK)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(SYK, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("SYK") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-3, 3), expand = c(0,0))
SYK_plot

### ANXA1
ANXA1_plot <- ggplot(scaled_data, aes(Type, ANXA1)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(ANXA1, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("ANXA1") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2, 4), expand = c(0,0))
ANXA1_plot

### Save images
pdf(file.path(path, "cohort_comparison.pdf"), height = 4.5, width = 5)
grid.arrange(SYK_plot, ANXA1_plot, nrow = 1)
dev.off()
