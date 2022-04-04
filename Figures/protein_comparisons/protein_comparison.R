library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(ggprism)

### Set variables
options(scipen = 100)
path <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Protein_comparison"
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

### Seperate out necessary 
ODG <- scaled_data[scaled_data$Type == "Type I", ]

### Calculate p-values for each subtype
group1 <- c("Normal", "Normal", "Normal", "Type I", "Type I", "Type II")
group2 <- c("Type I", "Type II", "Type III", "Type II" ,"Type III" ,"Type III")

### IDH1
a <- t.test(scaled_data$IDH1[scaled_data$Type=="Normal"], scaled_data$IDH1[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$IDH1[scaled_data$Type=="Normal"], scaled_data$IDH1[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$IDH1[scaled_data$Type=="Normal"], scaled_data$IDH1[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$IDH1[scaled_data$Type=="Type I"], scaled_data$IDH1[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$IDH1[scaled_data$Type=="Type I"], scaled_data$IDH1[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$IDH1[scaled_data$Type=="Type II"], scaled_data$IDH1[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
IDH1 <- as.data.frame(cbind(group1, group2, t_test))
IDH1$p <- ifelse(IDH1$t_test < 0.05 & IDH1$t_test > 0.01, "*",
                           ifelse(IDH1$t_test < 0.01 & IDH1$t_test > 0.001, "**",
                                  ifelse(IDH1$t_test < 0.001, "***", "")))
IDH1 <- IDH1[IDH1$p %in% c("*", "**", "***"), colnames(IDH1) %in% c("group1", "group2", "p")]
IDH1$y.position <- c(2, 2.5, 3)

### IDH2
a <- t.test(scaled_data$IDH2[scaled_data$Type=="Normal"], scaled_data$IDH2[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$IDH2[scaled_data$Type=="Normal"], scaled_data$IDH2[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$IDH2[scaled_data$Type=="Normal"], scaled_data$IDH2[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$IDH2[scaled_data$Type=="Type I"], scaled_data$IDH2[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$IDH2[scaled_data$Type=="Type I"], scaled_data$IDH2[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$IDH2[scaled_data$Type=="Type II"], scaled_data$IDH2[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
IDH2 <- as.data.frame(cbind(group1, group2, t_test))
IDH2$p <- ifelse(IDH2$t_test < 0.05 & IDH2$t_test > 0.01, "*",
                 ifelse(IDH2$t_test < 0.01 & IDH2$t_test > 0.001, "**",
                        ifelse(IDH2$t_test < 0.001, "***", "")))
IDH2 <- IDH2[IDH2$p %in% c("*", "**", "***"), colnames(IDH2) %in% c("group1", "group2", "p")]
IDH2$y.position <- c(2.5, 3, 3.5, 4, 4.5)

### TP53
a <- t.test(scaled_data$TP53[scaled_data$Type=="Normal"], scaled_data$TP53[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$TP53[scaled_data$Type=="Normal"], scaled_data$TP53[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$TP53[scaled_data$Type=="Normal"], scaled_data$TP53[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$TP53[scaled_data$Type=="Type I"], scaled_data$TP53[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$TP53[scaled_data$Type=="Type I"], scaled_data$TP53[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$TP53[scaled_data$Type=="Type II"], scaled_data$TP53[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
TP53 <- as.data.frame(cbind(group1, group2, t_test))
TP53$p <- ifelse(TP53$t_test < 0.05 & TP53$t_test > 0.01, "*",
                 ifelse(TP53$t_test < 0.01 & TP53$t_test > 0.001, "**",
                        ifelse(TP53$t_test < 0.001, "***", "")))
TP53 <- TP53[TP53$p %in% c("*", "**", "***"), colnames(TP53) %in% c("group1", "group2", "p")]
TP53$y.position <- c(1.1, 1.4)

### ATRX
a <- t.test(scaled_data$ATRX[scaled_data$Type=="Normal"], scaled_data$ATRX[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$ATRX[scaled_data$Type=="Normal"], scaled_data$ATRX[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$ATRX[scaled_data$Type=="Normal"], scaled_data$ATRX[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$ATRX[scaled_data$Type=="Type I"], scaled_data$ATRX[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$ATRX[scaled_data$Type=="Type I"], scaled_data$ATRX[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$ATRX[scaled_data$Type=="Type II"], scaled_data$ATRX[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
ATRX <- as.data.frame(cbind(group1, group2, t_test))
ATRX$p <- ifelse(ATRX$t_test < 0.05 & ATRX$t_test > 0.01, "*",
                 ifelse(ATRX$t_test < 0.01 & ATRX$t_test > 0.001, "**",
                        ifelse(ATRX$t_test < 0.001, "***", "")))
ATRX <- ATRX[ATRX$p %in% c("*", "**", "***"), colnames(ATRX) %in% c("group1", "group2", "p")]
ATRX$y.position <- c()

### TERT
a <- t.test(scaled_data$TERT[scaled_data$Type=="Normal"], scaled_data$TERT[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$TERT[scaled_data$Type=="Normal"], scaled_data$TERT[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$TERT[scaled_data$Type=="Normal"], scaled_data$TERT[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$TERT[scaled_data$Type=="Type I"], scaled_data$TERT[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$TERT[scaled_data$Type=="Type I"], scaled_data$TERT[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$TERT[scaled_data$Type=="Type II"], scaled_data$TERT[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
TERT <- as.data.frame(cbind(group1, group2, t_test))
TERT$p <- ifelse(TERT$t_test < 0.05 & TERT$t_test > 0.01, "*",
                 ifelse(TERT$t_test < 0.01 & TERT$t_test > 0.001, "**",
                        ifelse(TERT$t_test < 0.001, "***", "")))
TERT <- TERT[TERT$p %in% c("*", "**", "***"), colnames(TERT) %in% c("group1", "group2", "p")]
TERT$y.position <- c()

### CIC
a <- t.test(scaled_data$CIC[scaled_data$Type=="Normal"], scaled_data$CIC[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$CIC[scaled_data$Type=="Normal"], scaled_data$CIC[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$CIC[scaled_data$Type=="Normal"], scaled_data$CIC[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$CIC[scaled_data$Type=="Type I"], scaled_data$CIC[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$CIC[scaled_data$Type=="Type I"], scaled_data$CIC[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$CIC[scaled_data$Type=="Type II"], scaled_data$CIC[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
CIC <- as.data.frame(cbind(group1, group2, t_test))
CIC$p <- ifelse(CIC$t_test < 0.05 & CIC$t_test > 0.01, "*",
                 ifelse(CIC$t_test < 0.01 & CIC$t_test > 0.001, "**",
                        ifelse(CIC$t_test < 0.001, "***", "")))
CIC <- CIC[CIC$p %in% c("*", "**", "***"), colnames(CIC) %in% c("group1", "group2", "p")]
CIC$y.position <- c()

### EGFR
a <- t.test(scaled_data$EGFR[scaled_data$Type=="Normal"], scaled_data$EGFR[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$EGFR[scaled_data$Type=="Normal"], scaled_data$EGFR[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$EGFR[scaled_data$Type=="Normal"], scaled_data$EGFR[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$EGFR[scaled_data$Type=="Type I"], scaled_data$EGFR[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$EGFR[scaled_data$Type=="Type I"], scaled_data$EGFR[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$EGFR[scaled_data$Type=="Type II"], scaled_data$EGFR[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
EGFR <- as.data.frame(cbind(group1, group2, t_test))
EGFR$p <- ifelse(EGFR$t_test < 0.05 & EGFR$t_test > 0.01, "*",
                 ifelse(EGFR$t_test < 0.01 & EGFR$t_test > 0.001, "**",
                        ifelse(EGFR$t_test < 0.001, "***", "")))
EGFR <- EGFR[EGFR$p %in% c("*", "**", "***"), colnames(EGFR) %in% c("group1", "group2", "p")]
EGFR$y.position <- c(3, 3.7, 4.2)

### PTEN
a <- t.test(scaled_data$PTEN[scaled_data$Type=="Normal"], scaled_data$PTEN[scaled_data$Type=="Type I"])$p.value
b <- t.test(scaled_data$PTEN[scaled_data$Type=="Normal"], scaled_data$PTEN[scaled_data$Type=="Type II"])$p.value
c <- t.test(scaled_data$PTEN[scaled_data$Type=="Normal"], scaled_data$PTEN[scaled_data$Type=="Type III"])$p.value
d <- t.test(scaled_data$PTEN[scaled_data$Type=="Type I"], scaled_data$PTEN[scaled_data$Type=="Type II"])$p.value
e <- t.test(scaled_data$PTEN[scaled_data$Type=="Type I"], scaled_data$PTEN[scaled_data$Type=="Type III"])$p.value
f <- t.test(scaled_data$PTEN[scaled_data$Type=="Type II"], scaled_data$PTEN[scaled_data$Type=="Type III"])$p.value
t_test <- c(a, b, c, d, e, f)
PTEN <- as.data.frame(cbind(group1, group2, t_test))
PTEN$p <- ifelse(PTEN$t_test < 0.05 & PTEN$t_test > 0.01, "*",
                 ifelse(PTEN$t_test < 0.01 & PTEN$t_test > 0.001, "**",
                        ifelse(PTEN$t_test < 0.001, "***", "")))
PTEN <- PTEN[PTEN$p %in% c("*", "**", "***"), colnames(PTEN) %in% c("group1", "group2", "p")]
PTEN$y.position <- c(1.25)

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

### Calculate p-values CIC within ODGs
ODG$cic_mutation <- factor(ODG$cic_mutation, levels = c("WT", "Missense", "LOF"),
                           labels = c("Wildtype", "Missense", "LOF"))
group1 <- c("Wildtype", "Wildtype", "Missense")
group2 <- c("Missense", "LOF", "LOF")

### CIC
a <- t.test(ODG$CIC[ODG$cic_mutation=="Wildtype"], ODG$CIC[ODG$cic_mutation=="Missense"])$p.value
b <- t.test(ODG$CIC[ODG$cic_mutation=="Wildtype"], ODG$CIC[ODG$cic_mutation=="LOF"])$p.value
c <- t.test(ODG$CIC[ODG$cic_mutation=="Missense"], ODG$CIC[ODG$cic_mutation=="LOF"])$p.value
t_test <- c(a, b, c)
ODG_CIC <- as.data.frame(cbind(group1, group2, t_test))
ODG_CIC$p <- ifelse(ODG_CIC$t_test < 0.05 & ODG_CIC$t_test > 0.01, "*",
                 ifelse(ODG_CIC$t_test < 0.01 & ODG_CIC$t_test > 0.001, "**",
                        ifelse(ODG_CIC$t_test < 0.001, "***", "")))
ODG_CIC <- ODG_CIC[ODG_CIC$p %in% c("*", "**", "***"), colnames(ODG_CIC) %in% c("group1", "group2", "p")]
ODG_CIC$y.position <- c()

### ETV4
a <- t.test(ODG$ETV4[ODG$cic_mutation=="Wildtype"], ODG$ETV4[ODG$cic_mutation=="Missense"])$p.value
b <- t.test(ODG$ETV4[ODG$cic_mutation=="Wildtype"], ODG$ETV4[ODG$cic_mutation=="LOF"])$p.value
c <- t.test(ODG$ETV4[ODG$cic_mutation=="Missense"], ODG$ETV4[ODG$cic_mutation=="LOF"])$p.value
t_test <- c(a, b, c)
ODG_ETV4 <- as.data.frame(cbind(group1, group2, t_test))
ODG_ETV4$p <- ifelse(ODG_ETV4$t_test < 0.05 & ODG_ETV4$t_test > 0.01, "*",
                    ifelse(ODG_ETV4$t_test < 0.01 & ODG_ETV4$t_test > 0.001, "**",
                           ifelse(ODG_ETV4$t_test < 0.001, "***", "")))
ODG_ETV4 <- ODG_ETV4[ODG_ETV4$p %in% c("*", "**", "***"), colnames(ODG_ETV4) %in% c("group1", "group2", "p")]
ODG_ETV4$y.position <- c()

### ETV5
a <- t.test(ODG$ETV5[ODG$cic_mutation=="Wildtype"], ODG$ETV5[ODG$cic_mutation=="Missense"])$p.value
b <- t.test(ODG$ETV5[ODG$cic_mutation=="Wildtype"], ODG$ETV5[ODG$cic_mutation=="LOF"])$p.value
c <- t.test(ODG$ETV5[ODG$cic_mutation=="Missense"], ODG$ETV5[ODG$cic_mutation=="LOF"])$p.value
t_test <- c(a, b, c)
ODG_ETV5 <- as.data.frame(cbind(group1, group2, t_test))
ODG_ETV5$p <- ifelse(ODG_ETV5$t_test < 0.05 & ODG_ETV5$t_test > 0.01, "*",
                     ifelse(ODG_ETV5$t_test < 0.01 & ODG_ETV5$t_test > 0.001, "**",
                            ifelse(ODG_ETV5$t_test < 0.001, "***", "")))
ODG_ETV5 <- ODG_ETV5[ODG_ETV5$p %in% c("*", "**", "***"), colnames(ODG_ETV5) %in% c("group1", "group2", "p")]
ODG_ETV5$y.position <- c()

### DUSP6
a <- t.test(ODG$DUSP6[ODG$cic_mutation=="Wildtype"], ODG$DUSP6[ODG$cic_mutation=="Missense"])$p.value
b <- t.test(ODG$DUSP6[ODG$cic_mutation=="Wildtype"], ODG$DUSP6[ODG$cic_mutation=="LOF"])$p.value
c <- t.test(ODG$DUSP6[ODG$cic_mutation=="Missense"], ODG$DUSP6[ODG$cic_mutation=="LOF"])$p.value
t_test <- c(a, b, c)
ODG_DUSP6 <- as.data.frame(cbind(group1, group2, t_test))
ODG_DUSP6$p <- ifelse(ODG_DUSP6$t_test < 0.05 & ODG_DUSP6$t_test > 0.01, "*",
                      ifelse(ODG_DUSP6$t_test < 0.01 & ODG_DUSP6$t_test > 0.001, "**",
                             ifelse(ODG_DUSP6$t_test < 0.001, "***", "")))
ODG_DUSP6 <- ODG_DUSP6[ODG_DUSP6$p %in% c("*", "**", "***"), colnames(ODG_DUSP6) %in% c("group1", "group2", "p")]
ODG_DUSP6$y.position <- c()

### SPRY4
a <- t.test(ODG$SPRY4[ODG$cic_mutation=="Wildtype"], ODG$SPRY4[ODG$cic_mutation=="Missense"])$p.value
b <- t.test(ODG$SPRY4[ODG$cic_mutation=="Wildtype"], ODG$SPRY4[ODG$cic_mutation=="LOF"])$p.value
c <- t.test(ODG$SPRY4[ODG$cic_mutation=="Missense"], ODG$SPRY4[ODG$cic_mutation=="LOF"])$p.value
t_test <- c(a, b, c)
ODG_SPRY4 <- as.data.frame(cbind(group1, group2, t_test))
ODG_SPRY4$p <- ifelse(ODG_SPRY4$t_test < 0.05 & ODG_SPRY4$t_test > 0.01, "*",
                      ifelse(ODG_SPRY4$t_test < 0.01 & ODG_SPRY4$t_test > 0.001, "**",
                             ifelse(ODG_SPRY4$t_test < 0.001, "***", "")))
ODG_SPRY4 <- ODG_SPRY4[ODG_SPRY4$p %in% c("*", "**", "***"), colnames(ODG_SPRY4) %in% c("group1", "group2", "p")]
ODG_SPRY4$y.position <- c()

### Graph Abundance of proteins
### IDH1
IDH1_plot <- ggplot(scaled_data, aes(Type, IDH1)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(IDH1, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("IDH1") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2, 3.5), expand = c(0,0))
IDH1_plot

### IDH2
IDH2_plot <- ggplot(scaled_data, aes(Type, IDH2)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(IDH2, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("IDH2") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2.5, 5), expand = c(0,0))
IDH2_plot

### ATRX
ATRX_plot <- ggplot(scaled_data, aes(Type, ATRX)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  #add_pvalue(ATRX, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("ATRX") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-1.5, 3), expand = c(0,0))
ATRX_plot

### CIC
CIC_plot <- ggplot(scaled_data, aes(Type, CIC)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  #add_pvalue(CIC, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("CIC") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-1, 2), expand = c(0,0))
CIC_plot

### EGFR
EGFR_plot <- ggplot(scaled_data, aes(Type, EGFR)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(EGFR, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("EGFR") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-2.5, 5), expand = c(0,0))
EGFR_plot

### TERT
TERT_plot <- ggplot(scaled_data, aes(Type, TERT)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  #add_pvalue(TERT, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("TERT") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-1.25, 2), expand = c(0,0))
TERT_plot

### TP53
TP53_plot <- ggplot(scaled_data, aes(Type, TP53)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(TP53, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("TP53") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-1, 1.75), expand = c(0,0))
TP53_plot

### PTEN
PTEN_plot <- ggplot(scaled_data, aes(Type, PTEN)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(PTEN, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("PTEN") + 
  scale_fill_manual(labels = c("Normal", "Type I", "Type II", "Type III"), 
                    values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-1.1, 1.5), expand = c(0,0))
PTEN_plot

### CIC_ODG
CIC_ODG_plot <- ggplot(ODG, aes(cic_mutation, CIC)) +
  geom_boxplot(aes(fill = cic_mutation), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  #add_pvalue(OLIG2, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("CIC") + 
  scale_fill_manual(labels = c("Wildtype", "Missense", "LOF"), 
                    values = c("grey", "#5eb956", "#c18743")) +
  my_theme +
  scale_y_continuous(limits=c(-1, 0.5), expand = c(0,0))
CIC_ODG_plot

### ETV4_ODG
ETV4_ODG_plot <- ggplot(ODG, aes(cic_mutation, ETV4)) +
  geom_boxplot(aes(fill = cic_mutation), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  #add_pvalue(OLIG2, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("ETV4") + 
  scale_fill_manual(labels = c("Wildtype", "Missense", "LOF"), 
                    values = c("grey", "#5eb956", "#c18743")) +
  my_theme +
  scale_y_continuous(limits=c(-1, 1.5), expand = c(0,0))
ETV4_ODG_plot

### ETV5_ODG
ETV5_ODG_plot <- ggplot(ODG, aes(cic_mutation, ETV5)) +
  geom_boxplot(aes(fill = cic_mutation), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  #add_pvalue(OLIG2, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("ETV5") + 
  scale_fill_manual(labels = c("Wildtype", "Missense", "LOF"), 
                    values = c("grey", "#5eb956", "#c18743")) +
  my_theme +
  scale_y_continuous(limits=c(-1, 0.75), expand = c(0,0))
ETV5_ODG_plot

### DUSP6_ODG
DUSP6_ODG_plot <- ggplot(ODG, aes(cic_mutation, DUSP6)) +
  geom_boxplot(aes(fill = cic_mutation), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  #add_pvalue(OLIG2, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("DUSP6") + 
  scale_fill_manual(labels = c("Wildtype", "Missense", "LOF"), 
                    values = c("grey", "#5eb956", "#c18743")) +
  my_theme +
  scale_y_continuous(limits=c(-0.75, 0.5), expand = c(0,0))
DUSP6_ODG_plot

### SPRY4_ODG
SPRY4_ODG_plot <- ggplot(ODG, aes(cic_mutation, SPRY4)) +
  geom_boxplot(aes(fill = cic_mutation), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  #add_pvalue(OLIG2, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("SPRY4") + 
  scale_fill_manual(labels = c("Wildtype", "Missense", "LOF"), 
                    values = c("grey", "#5eb956", "#c18743")) +
  my_theme +
  scale_y_continuous(limits=c(-0.75, 0.75), expand = c(0,0))
SPRY4_ODG_plot


### Save images
pdf(file.path(path, "cohort_comparison.pdf"), height = 8, width = 12)
grid.arrange(IDH1_plot, IDH2_plot, TERT_plot, CIC_plot,
             ATRX_plot, TP53_plot, EGFR_plot, PTEN_plot, nrow = 2)
dev.off()

pdf(file.path(path, "ODG_comparison.pdf"), height = 4, width = 10)
grid.arrange(CIC_ODG_plot, ETV4_ODG_plot, ETV5_ODG_plot, DUSP6_ODG_plot, SPRY4_ODG_plot, nrow = 1)
dev.off()
