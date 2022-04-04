library(ggplot2)
library(dplyr)
library(ggrepel)
library(scales)

### Set Paths
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Proteome/output"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Protein_DE"

### Read in files
ODG_Astro <- read.delim(file.path(path, "DE_AstrocytomaVSOligodendroglioma.txt"))
ODG_GBM <- read.delim(file.path(path, "DE_GlioblastomaVSOligodendroglioma.txt"))
Astro_GBM <- read.delim(file.path(path, "DE_GlioblastomaVSAstrocytoma.txt"))

### Set functions
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

color_fun <- function(x) {
  x$threshold <- ifelse(x$log2FC >= 0.53 & x$adj.pval.deqms <= 0.05 , "A",
                       ifelse(x$log2FC <= -0.53 & x$adj.pval.deqms <= 0.05, "B", "C"))
  x
}

### Reverse the fold changes (analyses were opposite of the transcriptome ones)
ODG_Astro$log2FC <- ODG_Astro$log2FC*-1
ODG_GBM$log2FC <- ODG_GBM$log2FC*-1
Astro_GBM$log2FC <- Astro_GBM$log2FC*-1

### Annotate colors
ODG_Astro <- color_fun(ODG_Astro)
ODG_GBM <- color_fun(ODG_GBM)
Astro_GBM <- color_fun(Astro_GBM)

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
A <- ggplot(ODG_Astro, aes(x=log2FC, y=adj.pval.deqms)) +
  geom_point(aes(colour = threshold), size=1) +
  annotate("text", x = -0.5, y = 0.000001, label = "n = 24\nup = 13\ndown = 11") +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("") +
  xlab("Log2FoldChange") +
  ylab("p-adj") +
  my_theme +
  scale_y_continuous(trans=reverselog_trans(10))
A

B <- ggplot(ODG_GBM, aes(x=log2FC, y=adj.pval.deqms)) +
  geom_point(aes(colour = threshold), size=1) +
  annotate("text", x = -0.75, y = 0.00000003, label = "n = 197\nup = 93\ndown = 104") +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("") +
  xlab("Log2FoldChange") +
  ylab("p-adj") +
  my_theme +
  scale_y_continuous(trans=reverselog_trans(10))
B

C <- ggplot(Astro_GBM, aes(x=log2FC, y=adj.pval.deqms)) +
  geom_point(aes(colour = threshold), size=1) +
  annotate("text", x = -0, y = 0.0005, label = "n = 46\nup = 20\ndown = 26") +
  scale_colour_manual(values = c("A"= "#FE7568", "B"="#548FCC",  "C"= "grey")) +
  ggtitle("") +
  xlab("Log2FoldChange") +
  ylab("p-adj") +
  my_theme +
  scale_y_continuous(trans=reverselog_trans(10))
C

ggsave(file.path(outdir, "TypeI_vs_TypeII.pdf"), A, device = "pdf", width = 3, height = 5, units = "in")
ggsave(file.path(outdir, "TypeI_vs_TypeIII.pdf"), B, device = "pdf", width = 3, height = 5, units = "in")
ggsave(file.path(outdir, "TypeII_vs_TypeIII.pdf"), C, device = "pdf", width = 3, height = 5, units = "in")

