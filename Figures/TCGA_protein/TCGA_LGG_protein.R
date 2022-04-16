library(DESeq2)
library(scales)
library(dplyr)
library(ggplot2)
library(ggprism)
library(gridExtra)

### Set variables
path <- "/Volumes/GoogleDrive/My Drive/Post-Doc/Yip_Proteomics/TCGA_protein"

data_protein <- file.path(path, "TCGA_LGG_protein.txt")
data_samples <- file.path(path, "LGG-GBM Cohort.txt")

### Read in tables
data <- read.table(data_protein, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples <- read.table(data_samples, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

### Format data frames
rownames(data) <- data$Sample.REF
data <- data[2:nrow(data), 2:ncol(data)]
colnames(data) <- gsub("\\.", "-", colnames(data))
colnames(data) <- substr(colnames(data), 1, 15) 
data <- data[complete.cases(data), ]

data <- data[, colnames(data) %in% samples$TCGA.ID]
samples <- samples[samples$TCGA.ID %in% colnames(data), ]
row.names(samples) <- samples$TCGA.ID

data <- data[, samples$TCGA.ID]
samples <- samples[colnames(data), ]
data2 <- data

### Scale data values
data <- mutate_all(data, function(x) as.numeric(as.character(x)))
data <- apply(data, 1, rescale)
data <- t(as.data.frame(data))

### Transform into integer count matrix
data <- data*100
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}
data <- round_df(data, 0)

### Create dataframes for each subtype
samples_odg <- samples[samples$Type == "I", ]
samples_ast <- samples[samples$Type == "II", ]
samples_gbm <- samples[samples$Type == "III", ]

data_odg <- data[, colnames(data) %in% samples_odg$TCGA.ID]
data_ast <- data[, colnames(data) %in% samples_ast$TCGA.ID]
data_gbm <- data[, colnames(data) %in% samples_gbm$TCGA.ID]

### Combine for pairwise comparisons and factor
data_IvII <- cbind(data_odg, data_ast)
data_IvIII <- cbind(data_odg, data_gbm)
data_IIvIII <- cbind(data_ast, data_gbm)

samples_IvII <- rbind(samples_odg, samples_ast)
samples_IvII$Type <- factor(samples_IvII$Type, levels = c("I", "II"))

samples_IvIII <- rbind(samples_odg, samples_gbm)
samples_IvIII$Type <- factor(samples_IvIII$Type, levels = c("I", "III"))

samples_IIvIII <- rbind(samples_ast, samples_gbm)
samples_IIvIII$Type <- factor(samples_IIvIII$Type, levels = c("II", "III"))

### Run DeSeq2
# IvII
DESeqDataSet <- DESeqDataSetFromMatrix(countData = data_IvII,
                                       colData = samples_IvII,
                                       design = ~ Type)
DE <- DESeq(DESeqDataSet)
DE_results <- results(DE)
DE_results <- DE_results[order(DE_results$padj),]
DE_results <- as.data.frame(DE_results)
write.table(DE_results, file.path(path, "IvII_pDE.txt"), sep = "\t")

# IvIII
DESeqDataSet <- DESeqDataSetFromMatrix(countData = data_IvIII,
                                       colData = samples_IvIII,
                                       design = ~ Type)
DE <- DESeq(DESeqDataSet)
DE_results <- results(DE)
DE_results <- DE_results[order(DE_results$padj),]
DE_results <- as.data.frame(DE_results)
write.table(DE_results, file.path(path, "IvIII_pDE.txt"), sep = "\t")

# IIvIII
DESeqDataSet <- DESeqDataSetFromMatrix(countData = data_IvII,
                                       colData = samples_IvII,
                                       design = ~ Type)
DE <- DESeq(DESeqDataSet)
DE_results <- results(DE)
DE_results <- DE_results[order(DE_results$padj),]
DE_results <- as.data.frame(DE_results)
write.table(DE_results, file.path(path, "IIvIII_pDE.txt"), sep = "\t")

### Format dataframe for protein comparisons
data2 <- t(data2)
data2 <- merge(data2, samples, by = "row.names")

### Calculate p-values for each subtype
group1 <- c("I", "I", "II")
group2 <- c("II" ,"III" ,"III")

### EGFR
a <- t.test(data2$EGFR[data2$Type=="I"], data2$EGFR[data2$Type=="II"])$p.value
b <- t.test(data2$EGFR[data2$Type=="I"], data2$EGFR[data2$Type=="III"])$p.value
c <- t.test(data2$EGFR[data2$Type=="II"], data2$EGFR[data2$Type=="III"])$p.value
t_test <- c(a, b, c)
EGFR <- as.data.frame(cbind(group1, group2, t_test))
EGFR$p <- c("***", "***", "***")
EGFR$y.position <- c(115, 135, 125)

### p53
a <- t.test(data2$p53[data2$Type=="I"], data2$p53[data2$Type=="II"])$p.value
b <- t.test(data2$p53[data2$Type=="I"], data2$p53[data2$Type=="III"])$p.value
c <- t.test(data2$p53[data2$Type=="II"], data2$p53[data2$Type=="III"])$p.value
t_test <- c(a, b, c)
p53 <- as.data.frame(cbind(group1, group2, t_test))
p53$p <- c("***", "***", "***")
p53$y.position <- c(105, 125, 115)

### PTEN
a <- t.test(data2$PTEN[data2$Type=="I"], data2$PTEN[data2$Type=="II"])$p.value
b <- t.test(data2$PTEN[data2$Type=="I"], data2$PTEN[data2$Type=="III"])$p.value
c <- t.test(data2$PTEN[data2$Type=="II"], data2$PTEN[data2$Type=="III"])$p.value
t_test <- c(a, b, c)
PTEN <- as.data.frame(cbind(group1, group2, t_test))
PTEN$p <- c("*", "***", "*")
PTEN$y.position <- c(105, 125, 115)

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

### Graph Abundance of proteins
### EGFR
EGFR_plot <- ggplot(data2, aes(Type, EGFR)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(EGFR, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("EGFR") + 
  scale_fill_manual(labels = c("I", "II", "III"), 
                    values = c("#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-10, 140), expand = c(0,0))
EGFR_plot

### p53
p53_plot <- ggplot(data2, aes(Type, p53)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(p53, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("TP53") + 
  scale_fill_manual(labels = c("I", "II", "III"), 
                    values = c("#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-10, 130), expand = c(0,0))
p53_plot

### PTEN
PTEN_plot <- ggplot(data2, aes(Type, PTEN)) +
  geom_boxplot(aes(fill = Type), outlier.shape = NA) +
  geom_jitter(pch = 20, size = 3, alpha = 0.5, width = 0.15) +
  add_pvalue(PTEN, tip.length = 0, label.size = 6) +
  xlab("") + 
  ylab("Scaled Protein Abundance") +
  ggtitle("PTEN") + 
  scale_fill_manual(labels = c("I", "II", "III"), 
                    values = c("#5eb956", "#c18743", "#6f7dcb")) +
  my_theme +
  scale_y_continuous(limits=c(-10, 130), expand = c(0,0))
PTEN_plot

### Save figure
pdf(file.path(path, "TCGA_protein.pdf"), height = 4, width = 7)
grid.arrange(EGFR_plot, p53_plot, PTEN_plot, nrow = 1)
dev.off()
