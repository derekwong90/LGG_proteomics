library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(ggprism)

### Set variables
options(scipen = 100)
path <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Protein_comparison_2"
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

### Vimentin
a <- t.test(scaled_data$VIM[scaled_data$Type=="Normal"], scaled_data$VIM[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$VIM[scaled_data$Type=="Normal"], scaled_data$VIM[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$VIM[scaled_data$Type=="Normal"], scaled_data$VIM[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$VIM[scaled_data$Type=="Type I"], scaled_data$VIM[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$VIM[scaled_data$Type=="Type I"], scaled_data$VIM[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$VIM[scaled_data$Type=="Type II"], scaled_data$VIM[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
VIM <- as.data.frame(cbind(group1, group2, t_test))
VIM$p <- ifelse(VIM$t_test < 0.05 & VIM$t_test > 0.01, "*",
                 ifelse(VIM$t_test < 0.01 & VIM$t_test > 0.001, "**",
                        ifelse(VIM$t_test < 0.001, "***", "")))
VIM <- VIM[VIM$p %in% c("*", "**", "***"), colnames(VIM) %in% c("group1", "group2", "p")]
VIM$y.position <- c(2, 2.5, 3, 3.5, 4)

### OLIG2
a <- t.test(scaled_data$OLIG2[scaled_data$Type=="Normal"], scaled_data$OLIG2[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$OLIG2[scaled_data$Type=="Normal"], scaled_data$OLIG2[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$OLIG2[scaled_data$Type=="Normal"], scaled_data$OLIG2[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$OLIG2[scaled_data$Type=="Type I"], scaled_data$OLIG2[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$OLIG2[scaled_data$Type=="Type I"], scaled_data$OLIG2[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$OLIG2[scaled_data$Type=="Type II"], scaled_data$OLIG2[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
OLIG2 <- as.data.frame(cbind(group1, group2, t_test))
OLIG2$p <- ifelse(OLIG2$t_test < 0.05 & OLIG2$t_test > 0.01, "*",
                ifelse(OLIG2$t_test < 0.01 & OLIG2$t_test > 0.001, "**",
                       ifelse(OLIG2$t_test < 0.001, "***", "")))
OLIG2 <- OLIG2[OLIG2$p %in% c("*", "**", "***"), colnames(OLIG2) %in% c("group1", "group2", "p")]
OLIG2$y.position <- c(1, 1.35, 1.75)

### Nestin
a <- t.test(scaled_data$NES[scaled_data$Type=="Normal"], scaled_data$NES[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$NES[scaled_data$Type=="Normal"], scaled_data$NES[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$NES[scaled_data$Type=="Normal"], scaled_data$NES[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$NES[scaled_data$Type=="Type I"], scaled_data$NES[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$NES[scaled_data$Type=="Type I"], scaled_data$NES[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$NES[scaled_data$Type=="Type II"], scaled_data$NES[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
NES <- as.data.frame(cbind(group1, group2, t_test))
NES$p <- ifelse(NES$t_test < 0.05 & NES$t_test > 0.01, "*",
                  ifelse(NES$t_test < 0.01 & NES$t_test > 0.001, "**",
                         ifelse(NES$t_test < 0.001, "***", "")))
NES <- NES[NES$p %in% c("*", "**", "***"), colnames(NES) %in% c("group1", "group2", "p")]
NES$y.position <- c(1, 2.5, 3.5, 3, 4)

### BCAT1
a <- t.test(scaled_data$BCAT1[scaled_data$Type=="Normal"], scaled_data$BCAT1[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$BCAT1[scaled_data$Type=="Normal"], scaled_data$BCAT1[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$BCAT1[scaled_data$Type=="Normal"], scaled_data$BCAT1[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$BCAT1[scaled_data$Type=="Type I"], scaled_data$BCAT1[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$BCAT1[scaled_data$Type=="Type I"], scaled_data$BCAT1[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$BCAT1[scaled_data$Type=="Type II"], scaled_data$BCAT1[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
BCAT1 <- as.data.frame(cbind(group1, group2, t_test))
BCAT1$p <- ifelse(BCAT1$t_test < 0.05 & BCAT1$t_test > 0.01, "*",
                ifelse(BCAT1$t_test < 0.01 & BCAT1$t_test > 0.001, "**",
                       ifelse(BCAT1$t_test < 0.001, "***", "")))
BCAT1 <- BCAT1[BCAT1$p %in% c("*", "**", "***"), colnames(BCAT1) %in% c("group1", "group2", "p")]
BCAT1$y.position <- c(2, 4.5, 5.5)

### S100A1
a <- t.test(scaled_data$S100A1[scaled_data$Type=="Normal"], scaled_data$S100A1[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$S100A1[scaled_data$Type=="Normal"], scaled_data$S100A1[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$S100A1[scaled_data$Type=="Normal"], scaled_data$S100A1[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$S100A1[scaled_data$Type=="Type I"], scaled_data$S100A1[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$S100A1[scaled_data$Type=="Type I"], scaled_data$S100A1[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$S100A1[scaled_data$Type=="Type II"], scaled_data$S100A1[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
S100A1 <- as.data.frame(cbind(group1, group2, t_test))
S100A1$p <- ifelse(S100A1$t_test < 0.05 & S100A1$t_test > 0.01, "*",
                  ifelse(S100A1$t_test < 0.01 & S100A1$t_test > 0.001, "**",
                         ifelse(S100A1$t_test < 0.001, "***", "")))
S100A1 <- S100A1[S100A1$p %in% c("*", "**", "***"), colnames(S100A1) %in% c("group1", "group2", "p")]
S100A1$y.position <- c(2.25, 2.75, 3.25)

### KCNN3
a <- t.test(scaled_data$KCNN3[scaled_data$Type=="Normal"], scaled_data$KCNN3[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$KCNN3[scaled_data$Type=="Normal"], scaled_data$KCNN3[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$KCNN3[scaled_data$Type=="Normal"], scaled_data$KCNN3[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$KCNN3[scaled_data$Type=="Type I"], scaled_data$KCNN3[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$KCNN3[scaled_data$Type=="Type I"], scaled_data$KCNN3[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$KCNN3[scaled_data$Type=="Type II"], scaled_data$KCNN3[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
KCNN3 <- as.data.frame(cbind(group1, group2, t_test))
KCNN3$p <- ifelse(KCNN3$t_test < 0.05 & KCNN3$t_test > 0.01, "*",
                   ifelse(KCNN3$t_test < 0.01 & KCNN3$t_test > 0.001, "**",
                          ifelse(KCNN3$t_test < 0.001, "***", "")))
KCNN3 <- KCNN3[KCNN3$p %in% c("*", "**", "***"), colnames(KCNN3) %in% c("group1", "group2", "p")]
KCNN3$y.position <- c(3.5, 4, 4.5)

### MUCL1
a <- t.test(scaled_data$MUCL1[scaled_data$Type=="Normal"], scaled_data$MUCL1[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$MUCL1[scaled_data$Type=="Normal"], scaled_data$MUCL1[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$MUCL1[scaled_data$Type=="Normal"], scaled_data$MUCL1[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$MUCL1[scaled_data$Type=="Type I"], scaled_data$MUCL1[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$MUCL1[scaled_data$Type=="Type I"], scaled_data$MUCL1[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$MUCL1[scaled_data$Type=="Type II"], scaled_data$MUCL1[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
MUCL1 <- as.data.frame(cbind(group1, group2, t_test))
MUCL1$p <- ifelse(MUCL1$t_test < 0.05 & MUCL1$t_test > 0.01, "*",
                  ifelse(MUCL1$t_test < 0.01 & MUCL1$t_test > 0.001, "**",
                         ifelse(MUCL1$t_test < 0.001, "***", "")))
MUCL1 <- MUCL1[MUCL1$p %in% c("*", "**", "***"), colnames(MUCL1) %in% c("group1", "group2", "p")]
MUCL1$y.position <- c(4.75, 5.25, 5.75, 6.25)

### PRNP
a <- t.test(scaled_data$PRNP[scaled_data$Type=="Normal"], scaled_data$PRNP[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$PRNP[scaled_data$Type=="Normal"], scaled_data$PRNP[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$PRNP[scaled_data$Type=="Normal"], scaled_data$PRNP[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$PRNP[scaled_data$Type=="Type I"], scaled_data$PRNP[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$PRNP[scaled_data$Type=="Type I"], scaled_data$PRNP[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$PRNP[scaled_data$Type=="Type II"], scaled_data$PRNP[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
PRNP <- as.data.frame(cbind(group1, group2, t_test))
PRNP$p <- ifelse(PRNP$t_test < 0.05 & PRNP$t_test > 0.01, "*",
                  ifelse(PRNP$t_test < 0.01 & PRNP$t_test > 0.001, "**",
                         ifelse(PRNP$t_test < 0.001, "***", "")))
PRNP <- PRNP[PRNP$p %in% c("*", "**", "***"), colnames(PRNP) %in% c("group1", "group2", "p")]
PRNP$y.position <- c(2.25, 2.75, 3.25, 3.75, 4.25)

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

### VIM
VIM_plot <- ggplot(scaled_data, aes(Type, VIM)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(VIM, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("VIM") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2, 5), expand = c(0,0))
VIM_plot

### OLIG2
OLIG2_plot <- ggplot(scaled_data, aes(Type, OLIG2)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(OLIG2, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("OLIG2") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-1, 2), expand = c(0,0))
OLIG2_plot

### Nestin
NES_plot <- ggplot(scaled_data, aes(Type, NES)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(NES, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("NES") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2.5, 5), expand = c(0,0))
NES_plot

### BCAT1
BCAT1_plot <- ggplot(scaled_data, aes(Type, BCAT1)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(BCAT1, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("BCAT1") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2, 6), expand = c(0,0))
BCAT1_plot

### S100A1
S100A1_plot <- ggplot(scaled_data, aes(Type, S100A1)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(S100A1, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("S100A1") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2.5, 4), expand = c(0,0))
S100A1_plot

### KCNN3
KCNN3_plot <- ggplot(scaled_data, aes(Type, KCNN3)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(KCNN3, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("KCNN3") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-1.25, 5), expand = c(0,0))
KCNN3_plot

### MUCL1
MUCL1_plot <- ggplot(scaled_data, aes(Type, MUCL1)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(MUCL1, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("MUCL1") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2, 7), expand = c(0,0))
MUCL1_plot

### PRNP
PRNP_plot <- ggplot(scaled_data, aes(Type, PRNP)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(PRNP, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("PRNP") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2, 5), expand = c(0,0))
PRNP_plot

### Save images
pdf(file.path(path, "cohort_comparison.pdf"), height = 18, width = 7)
grid.arrange(VIM_plot, NES_plot, OLIG2_plot, S100A1_plot,
             BCAT1_plot, MUCL1_plot, KCNN3_plot, PRNP_plot, nrow = 4)
dev.off()
