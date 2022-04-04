library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales)

### Set Paths
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Transcriptome"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/TCGA_DE"

### Read in files
I_II <- read.delim(file.path(path, "IvsII.txt"))
I_III <- read.delim(file.path(path, "IvsIII.txt"))
II_III <- read.delim(file.path(path, "IIvsIII.txt"))

WT_LOF <- read.delim(file.path(path, "WTvsLOF.txt"))
WT_Mis <- read.delim(file.path(path, "WTvsMissense.txt"))
LOF_Mis <- read.delim(file.path(path, "LOFvsMissense.txt"))

### Remove lines of NA
I_II <- I_II[complete.cases(I_II[1:9]), 1:9]
I_III<- I_III[complete.cases(I_III[1:9]), 1:9]
II_III <- II_III[complete.cases(II_III[1:9]), 1:9]

WT_LOF <- WT_LOF[complete.cases(WT_LOF[1:9]), 1:9]
WT_Mis <- WT_Mis[complete.cases(WT_Mis)[1:9], 1:9]
LOF_Mis <- LOF_Mis[complete.cases(LOF_Mis[1:9]), 1:9]

### Set functions
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

color_fun <- function(x) {
  x$threshold <- ifelse(x$log2FoldChange >= 1 & x$padj <= 0.01 , "A",
                       ifelse(x$log2FoldChange <= -1 & x$padj <= 0.01, "B", "C"))
  x
}

color_fun2 <- function(x) {
  x$threshold <- ifelse(x$log2FoldChange >= 0 & x$padj <= 0.05 , "A",
                        ifelse(x$log2FoldChange <= 0 & x$padj <= 0.05, "B", "C"))
  x
}

### Annotate colors
I_II <- color_fun(I_II)
I_III <- color_fun(I_III)
II_III <- color_fun(II_III)
WT_LOF <- color_fun2(WT_LOF)
WT_Mis <- color_fun2(WT_Mis)
LOF_Mis <- color_fun2(LOF_Mis)

### Set theme
my_theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
                  axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  legend.position="none",
                  axis.text = element_text(size = 15),
                  axis.title = element_text(size = 15),
                  axis.text.x = element_text(angle = 45, hjust = 1))

### Plot volcano plots
a <- ggplot(I_II, aes(x=log2FoldChange, y=padj)) +
  geom_point(aes(colour = threshold), size=1) +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("") +
  xlab("Log2FoldChange") +
  ylab("p-adj") +
  my_theme +
  scale_y_continuous(trans=reverselog_trans(10))
a

b <- ggplot(I_III, aes(x=log2FoldChange, y=padj)) +
  geom_point(aes(colour = threshold), size=1) +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("") +
  xlab("Log2FoldChange") +
  ylab("p-adj") +
  my_theme +
  scale_y_continuous(trans=reverselog_trans(10))
b

c <- ggplot(II_III, aes(x=log2FoldChange, y=padj)) +
  geom_point(aes(colour = threshold), size=1) +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("") +
  xlab("Log2FoldChange") +
  ylab("p-adj") +
  my_theme +
  scale_y_continuous(trans=reverselog_trans(10))
c

d <- ggplot(WT_LOF, aes(x=log2FoldChange, y=padj)) +
  geom_point(aes(colour = threshold), size=1) +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("") +
  xlab("Log2FoldChange") +
  ylab("p-adj") +
  my_theme +
  scale_y_continuous(trans=reverselog_trans(10))
d

e <- ggplot(WT_Mis, aes(x=log2FoldChange, y=padj)) +
  geom_point(aes(colour = threshold), size=1) +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("") +
  xlab("Log2FoldChange") +
  ylab("p-adj") +
  my_theme +
  scale_y_continuous(trans=reverselog_trans(10))
e

f <- ggplot(LOF_Mis, aes(x=log2FoldChange, y=padj)) +
  geom_point(aes(colour = threshold), size=1) +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("") +
  xlab("Log2FoldChange") +
  ylab("p-adj") +
  my_theme +
  scale_y_continuous(trans=reverselog_trans(10))
f

### Save Volcano Plots
ggsave(file.path(outdir, "TypeI_vs_TypeII.pdf"), a, device = "pdf", width = 3, height = 5, units = "in")
ggsave(file.path(outdir, "TypeI_vs_TypeIII.pdf"), b, device = "pdf", width = 3, height = 5, units = "in")
ggsave(file.path(outdir, "TypeII_vs_TypeIII.pdf"), c, device = "pdf", width = 3, height = 5, units = "in")
ggsave(file.path(outdir, "WT_vs_LOF.pdf"), e, device = "pdf", width = 3, height = 5, units = "in")
ggsave(file.path(outdir, "WT_vs_Missense.pdf"), e, device = "pdf", width = 3, height = 5, units = "in")
ggsave(file.path(outdir, "LOF_vs_Missense.pdf"), f, device = "pdf", width = 3, height = 5, units = "in")

### Find overlap between CIC analyses
WT_LOF <- WT_LOF[WT_LOF$padj < 0.05, ]
WT_Mis <- WT_Mis[WT_Mis$padj < 0.05, ]
Overlap <- merge(WT_LOF, WT_Mis, by = "Gene")
Overlap$threshold <- ifelse(Overlap$log2FoldChange.x > 0 & Overlap$log2FoldChange.y > 0 , "A",
                            ifelse(Overlap$log2FoldChange.x < 0 & Overlap$log2FoldChange.y < 0, "B", "C"))
Overlap$label <- ifelse(Overlap$Gene %in% c("ETV1", "ETV4", "ETV5", "DUSP6", "DUSP4", "SPRY4", "SHC3", "SPRED1", "SPRED2"), Overlap$Gene, NA)

### Plot Overlap
g <- ggplot(Overlap, aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  stat_smooth(method = "lm", col = "grey", se = FALSE) +
  geom_point(aes(colour = threshold), size=1) +
  geom_text_repel(aes(label = label),
                   nudge_y       = -1.5,
                   size          = 3,
                   box.padding   = 0.5,
                   point.padding = 0.5,
                   force         = 100,
                   segment.size  = 0.2,
                   segment.color = "grey50") +
  stat_cor(label.x = -2.75, label.y = 1) +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("DE Genes") +
  xlab("CIC LOF") +
  ylab("CIC Missense") +
  my_theme +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0)
g
ggsave(file.path(outdir, "Overlap.pdf"), g, device = "pdf", width = 5, height = 5, units = "in")
