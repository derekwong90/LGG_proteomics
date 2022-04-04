setwd("C:/Users/user/Google Drive/Grad School/LGG-Proteomics/TCGA Analyses")

library(biomaRt)
library(chromPlot)
library(org.Hs.eg.db)
library(dplyr)

dat <- read.delim("data.txt")

dat <- dat[,1:9]

#Get Chromosome information and chromPlot files
ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
chromosomes <- getBM(attributes=c('entrezgene_id', 'start_position', 'end_position', 'chromosome_name'), mart = ensembl)
data(hg_gap)

#Merge chromosome information and make plot dataframe
chroms <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
dat_chrom <- merge(dat, chromosomes, by.x = "Entrez.ID", by.y = "entrezgene_id")
dat_chrom <- dat_chrom[dat_chrom$chromosome_name  %in% chroms,]
dat_plot <- data.frame(Gene = as.character(dat_chrom$Gene),
                       Chrom = as.character(dat_chrom$chromosome_name),
                       Start = as.integer(dat_chrom$start_position),
                       End = as.integer(dat_chrom$end_position),
                       padj = as.numeric(dat_chrom$padj),
                       log2foldchange = as.numeric(dat_chrom$log2FoldChange))


#Filter out genes for plot and assign placements
dat_chrom_plot <- dat_plot[dat_chrom$padj < 0.01, ]
dat_chrom_plot <- dat_chrom_plot[dat_chrom_plot$log2foldchange > 0.58 | dat_chrom_plot$log2foldchange < -0.58,]
dat_chrom_plot <- dat_chrom_plot%>%mutate(threshold = ifelse(dat_chrom_plot$log2foldchange > 0, "pos", "neg"))
dat_chrom_plot_pos <- subset(dat_chrom_plot, threshold %in% "pos")
dat_chrom_plot_neg <- subset(dat_chrom_plot, threshold %in% "neg")

#chromPlot function
#histogram
chromPlot(gaps = hg_gap, 
          annot1 = dat_chrom_plot_pos, 
          annot2 = dat_chrom_plot_neg,
          colAnnot1 = c("#cc5f43"),
          colAnnot2 = c("#8176cc"),
          chrSide=c(-1, 1, 1, 1, 1, 1, 1, 1))

#scatterplot
dat_chrom_plot_neg$log2foldchange <- abs(dat_chrom_plot_neg$log2foldchange)

chromPlot(gaps = hg_gap,
          stat = dat_chrom_plot_pos, 
          stat2 = dat_chrom_plot_neg,
          statCol = "log2foldchange",
          statCol2 = "log2foldchange",
          statName = "log2foldchange",
          statName2 = "log2foldchange",
          statTyp = "p",
          scex = 0.7,
          spty = 20,
          statSumm = "none",
          colStat = c("#cc5f43"),
          colStat2 = c("#8176cc"),
          chrSide=c(-1, 1, 1, 1, 1, 1, -1, 1))


#chromosome comparisons

library(ggplot2)

data <- data.frame(Gene = as.character(dat_chrom$Gene),
                    Chrom = as.character(dat_chrom$chromosome),
                    chromosome = as.character(dat_chrom$chromosome_name),
                    padj = as.numeric(dat_chrom$padj),
                    log2foldchange = as.numeric(dat_chrom$log2FoldChange))

data <- data[data$padj < 0.01, ]
data <- data[data$log2foldchange > 0.58 | data$log2foldchange < -0.58,]
data <- data%>%mutate(threshold = ifelse(data$log2foldchange > 0, "pos", "neg"))
data <- subset(data, chromosome %in% c("1", "19"))
data_pos <- subset(data, threshold %in% "pos")
data_neg <- subset(data, threshold %in% "neg")
write.table(data_pos, "pos.txt", sep = "\t")
write.table(data_neg, "neg.txt", sep = "\t")

data_pos <- read.delim("pos.txt")
data_neg <- read.delim("neg.txt")

data_pos <- data_pos$threshold
data_neg <- data_neg$threshold

pos <- as.data.frame(table(data_pos))
neg <- as.data.frame(table(data_neg))

data <- data.frame(chromosome = as.character(neg$data_neg),
                   frequency = as.numeric(pos$Freq/neg$Freq))

data$chromosome = factor(data$chromosome, 
                     levels = c("1p", "1q", "19p", "19q"), 
                     labels = c("1p", "1q", "19p", "19q"))
data$frequency <- as.numeric(format(round(data$frequency, digits = 2)))

ggplot(data, aes(x = chromosome, y = frequency, fill = chromosome)) +
  geom_bar(color = "Black", stat="identity", 
           width = 0.8, position = position_dodge(width = 1)) +
  geom_text(aes(label=frequency), position = position_dodge(width = 0.9), vjust = -0.25) +
  xlab("Chromosome") + 
  ylab("DE Genes positive/negative") +
  ggtitle("") + 
  scale_fill_manual(values=c("#cc5f43", "#8176cc", "#cc5f43", "#8176cc")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
        axis.line = element_line(colour = "black"),
        legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15)) +
  scale_y_continuous(limits=c(0, 6), expand = c(0,0))






