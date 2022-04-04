library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(matrixStats)
library(ggpubr)

### Set Paths
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Proteome/output"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Protein_CNV"

#Get Chromosome information and chromPlot files
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = 'useast')
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
chromosomes <- getBM(attributes=c("external_gene_name", 'uniprot_gn_id', 'start_position', 'end_position', 'chromosome_name', 'band'), mart = ensembl)
chromosomes <- chromosomes[!(duplicated(chromosomes$external_gene_name)), ]

### Read in files
data <- read.delim(file.path(path, "protein_abundance.txt"))
samples <- read.delim(file.path(path, "metadata.txt"))

### Format data
data <- data[complete.cases(data), ]
samples <- samples[samples$Type %in% c("Oligodendroglioma", "Astrocytoma", "Glioblastoma", "Normal"), ]
samples <- samples[order(factor(samples$molecular, levels = c("", "I", "II", "III"))), ]
samples$molecular <- ifelse(samples$molecular == "", "Normal",
                            ifelse(samples$molecular == "I", "Type I",
                                   ifelse(samples$molecular == "II", "Type II", "Type III")))

#Merge chromosome information and make plot dataframe
#chroms <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
chroms <- c("1", "7", "10", "19")
data <- merge(data, chromosomes, by.x = "Gene", by.y = "external_gene_name")
data <- data[data$chromosome_name  %in% chroms,]

### Order dataframe
data <- data[!(duplicated(data$Gene)), ]
data <- data[order(data$chromosome_name,
                   data$start_position), ]

### Order data and match samples
rownames(data) <- data$Gene
chr_split <- data$chromosome_name
chr_split <- factor(chr_split, levels = chroms)
data <- data[ , colnames(data) %in% samples$set]
data <- data[ , samples$set]

### Normalize data
normal_median <- rowMeans(data[, 1:6])
normal_sd <- rowSds(as.matrix(data[, 1:6]))
data <- (data - normal_median)/normal_sd

### Calculate Z-scores
codel1 <- chromosomes[chromosomes$chromosome_name == "1", ]
codel1 <- codel1[grepl("p", codel1$band), ]
codel19 <- chromosomes[chromosomes$chromosome_name == "19", ]
codel19 <- codel19[grepl("q", codel19$band), ]
codel <- rbind(codel1, codel19)
rm(codel1, codel19)

codel_data <- data[rownames(data) %in% codel$external_gene_name, ]
codel_healthy <- codel_data[, 1:6]
codel_sums <- colSums(codel_healthy)
codel_mean <- mean(codel_sums)
codel_sd <- sd(codel_sums)

codel_sums <- colSums(codel_data)
codel_scores <- (codel_sums - codel_mean)/codel_sd
samples$scores <- codel_scores

### Set orders and splits
samples <- samples[order(samples$scores, decreasing = TRUE), ]
data <- data[, samples$set]
type_split <- samples$molecular
type_split <- factor(type_split, levels = c("Normal", "Type I", "Type II", "Type III"))

### Make annotations
top_annotation <- HeatmapAnnotation("1p19q\nZ-score" = anno_points(samples$scores))

### Legend Param
heatmap_legend_param <-list(title = "SD from Healthy",
                            border = FALSE,
                            at = c(-20, -10, 0, 10, 20))

### Make the Heatmap
pdf(file.path(outdir, "cnv_heatmap.pdf"), height = 6, width = 6)
Heatmap <- Heatmap(data,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   row_order = row.names(data),
                   column_order = colnames(data),
                   row_split = chr_split,
                   column_split = type_split,
                   top_annotation = top_annotation,
                   heatmap_legend_param = heatmap_legend_param)

draw(Heatmap)
dev.off()

### Compare scores
my_comparisons <- list( c("Normal", "Type I"), 
                        c("Normal", "Type II"),
                        c("Normal", "Type III"),
                        c("Type I", "Type II"),
                        c("Type I", "Type III"))

score_plot <- ggplot(samples, aes(x = molecular, y = scores, fill = molecular)) +
  geom_boxplot() +
  geom_jitter(pch = 20, size = 2, alpha = 0.5, width = 0.15) +
  xlab("Molecular Subtype") + 
  ylab("Z-score") +
  scale_fill_manual(values = c("grey", "#5eb956", "#c18743", "#6f7dcb")) +
  ggtitle("1p-19q Z-Score") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  scale_y_continuous(limits=c(-12, 15), expand = c(0,0)) + 
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     tip.length = 0,
                     step.increase = 0.15,
                     position = 0.3)
score_plot
ggsave(file.path(outdir, "1p19q_scores.pdf"), score_plot, device = "pdf", width = 4, height = 6, units = "in")


