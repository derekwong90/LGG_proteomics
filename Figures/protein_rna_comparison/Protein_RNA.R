library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(matrixStats)
library(ggpubr)

### Set Paths
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Proteome/output"
path2 <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Transcriptome"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Protein_RNA_comparison"

#Get Chromosome information and chromPlot files
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = 'useast')
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
chromosomes <- getBM(attributes=c("external_gene_name", 'uniprot_gn_id', 'start_position', 'end_position', 'chromosome_name', 'band'), mart = ensembl)
chromosomes <- chromosomes[!(duplicated(chromosomes$external_gene_name)), ]

### Read in files
RNA <- read.delim(file.path(path2, "LGG.uncv2.mRNAseq_RSEM_all.txt"))
master <- read.delim(file.path(path2, "LGG-GBM Cohort.txt"))

data <- read.delim(file.path(path, "protein_abundance.txt"))
samples <- read.delim(file.path(path, "metadata.txt"))

### Format data
data <- data[complete.cases(data), ]
samples <- samples[samples$Type %in% c("Oligodendroglioma"), ]

#Merge chromosome information and make plot dataframe
chroms <- c("1", "19")
data <- merge(data, chromosomes, by.x = "Gene", by.y = "external_gene_name")
data <- data[data$chromosome_name  %in% chroms,]

### Order dataframe
data <- data[!(duplicated(data$Gene)), ]
data <- data[order(data$chromosome_name,
                   data$start_position), ]

### Order data and match samples
rownames(data) <- data$Gene
data <- data[ , colnames(data) %in% samples$set]
data <- data[ , samples$set]

### Keep only 1p 19q
codel1 <- chromosomes[chromosomes$chromosome_name == "1", ]
#codel1 <- codel1[grepl("p", codel1$band), ]
codel19 <- chromosomes[chromosomes$chromosome_name == "19", ]
#codel19 <- codel19[grepl("q", codel19$band), ]
codel <- rbind(codel1, codel19)
rm(codel1, codel19)
data <- data[rownames(data) %in% codel$external_gene_name, ]

### Format RNA expression like protein
master <- master[master$Type == "I", ]
RNA$Gene.ID <- strsplit(as.character(RNA$HYBRIDIZATION.R), split="[|]")
RNA$Gene <- sapply(RNA$Gene.ID, "[[", 1)
RNA <- RNA[!(duplicated(RNA$Gene)), ]
rownames(RNA) <- RNA$Gene
RNA <- RNA[row.names(RNA) %in% rownames(data), colnames(RNA) %in% master$TCGA.ID]
data <- data[row.names(data) %in% row.names(RNA), ]
data <- data[rownames(RNA), ]

### Calculate medians for RNA and protein
protein <- rowMeans(data)
expression <- rowMeans(RNA)
comparison <- as.data.frame(cbind(protein, expression))

### Scale data
comparison <- comparison[!(row.names(comparison) %in% c("APOE", "BCAN", "EEF2", "PEA15", "GLUL")), ] #Expression outliers
comparison$protein <- (comparison$protein - min(comparison$protein)) / (max(comparison$protein) - min(comparison$protein))
comparison$expression <- (comparison$expression - min(comparison$expression)) / (max(comparison$expression) - min(comparison$expression))
comparison <- merge(comparison, codel, by.x = "row.names", by.y = "external_gene_name")

### Calculate regression
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(R)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
regression <- comparison[, 2:3]
colnames(regression) <- c("x", "y")
R2 <- lm_eqn(regression)

### Format for each arm
comparison$arm <- substr(comparison$band, 1, 1)
comparison$chromosome <- paste0(comparison$chromosome_name, comparison$arm)
comparison$chromosome <- factor(comparison$chromosome, levels = c("1p", "1q", "19p", "19q"))

### Plot data
score_plot <- ggplot(comparison, aes(x = protein, y = expression)) +
  geom_point(aes(fill = "grey"), pch = 16, alpha = 0.5) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) + 
  stat_cor(aes(label=..rr.label..), label.x=0, label.y=0.95, size = 5) +
  facet_wrap(.~chromosome, scales = "fixed", nrow = 2) +
  xlab("Protein Expression (Relative Abundance)") + 
  ylab("mRNA Expression (RSEM)") +
  labs(fill = "Chromosome") +
  #scale_fill_manual(values = c("red", "blue"), labels = c("1p", "19q")) +
  ggtitle("Protein vs mRNA Expression") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 15, face = "bold"),
        legend.position = "none",
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 15))

score_plot
ggsave(file.path(outdir, "Comparison.pdf"), score_plot, device = "pdf", width = 5, height = 6, units = "in")


