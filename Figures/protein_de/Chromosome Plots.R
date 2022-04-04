library(biomaRt)
library(chromPlot)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)

### Set Paths
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Proteome/output"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Protein_DE"

#Get Chromosome information and chromPlot files
ensembl = useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
chromosomes <- getBM(attributes=c('uniprot_gn_id', 'start_position', 'end_position', 'chromosome_name', 'band'), mart = ensembl)
data(hg_gap)

### Read in files
II_III <- read.delim(file.path(path, "DE_AstrocytomaVSOligodendroglioma.txt"))
II_IIII <- read.delim(file.path(path, "DE_GlioblastomaVSOligodendroglioma.txt"))
III_IIII <- read.delim(file.path(path, "DE_GlioblastomaVSAstrocytoma.txt"))

#Merge chromosome information and make plot dataframe
chroms <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")

II_III_chrom <- merge(II_III, chromosomes, by.x = "uniprot", by.y = "uniprot_gn_id")
II_III_chrom <- II_III_chrom[II_III_chrom$chromosome_name  %in% chroms,]
II_III_plot <- data.frame(Gene = as.character(II_III_chrom$name),
                        Chrom = as.character(II_III_chrom$chromosome_name),
                        Start = as.integer(II_III_chrom$start_position),
                        End = as.integer(II_III_chrom$end_position),
                        padj = as.numeric(II_III_chrom$adj.pval.deqms),
                        log2foldchange = as.numeric(II_III_chrom$log2FC),
                        band = as.character(substring(II_III_chrom$band, 1, 1)))

II_IIII_chrom <- merge(II_IIII, chromosomes, by.x = "uniprot", by.y = "uniprot_gn_id")
II_IIII_chrom <- II_IIII_chrom[II_IIII_chrom$chromosome_name  %in% chroms,]
II_IIII_plot <- data.frame(Gene = as.character(II_IIII_chrom$name),
                        Chrom = as.character(II_IIII_chrom$chromosome_name),
                        Start = as.integer(II_IIII_chrom$start_position),
                        End = as.integer(II_IIII_chrom$end_position),
                        padj = as.numeric(II_IIII_chrom$adj.pval.deqms),
                        log2foldchange = as.numeric(II_IIII_chrom$log2FC),
                        band = as.character(substring(II_IIII_chrom$band, 1, 1)))

III_IIII_chrom <- merge(III_IIII, chromosomes, by.x = "uniprot", by.y = "uniprot_gn_id")
III_IIII_chrom <- III_IIII_chrom[III_IIII_chrom$chromosome_name  %in% chroms,]
III_IIII_plot <- data.frame(Gene = as.character(III_IIII_chrom$name),
                        Chrom = as.character(III_IIII_chrom$chromosome_name),
                        Start = as.integer(III_IIII_chrom$start_position),
                        End = as.integer(III_IIII_chrom$end_position),
                        padj = as.numeric(III_IIII_chrom$adj.pval.deqms),
                        log2foldchange = as.numeric(III_IIII_chrom$log2FC),
                        band = as.character(substring(III_IIII_chrom$band, 1, 1)))

#Filter out genes for plot and assign placements
II_III_plot <- II_III_plot[II_III_plot$padj < 0.05, ]
II_III_plot$threshold <- ifelse(II_III_plot$log2foldchange > 0, "pos", "neg")
II_III_plot_pos <- subset(II_III_plot, threshold %in% "pos")
II_III_plot_neg <- subset(II_III_plot, threshold %in% "neg")

II_IIII_plot <- II_IIII_plot[II_IIII_plot$padj < 0.05, ]
II_IIII_plot$threshold <- ifelse(II_IIII_plot$log2foldchange > 0, "pos", "neg")
II_IIII_plot_pos <- subset(II_IIII_plot, threshold %in% "pos")
II_IIII_plot_neg <- subset(II_IIII_plot, threshold %in% "neg")

III_IIII_plot <- III_IIII_plot[III_IIII_plot$padj < 0.05, ]
III_IIII_plot$threshold <- ifelse(III_IIII_plot$log2foldchange > 0, "pos", "neg")
III_IIII_plot_pos <- subset(III_IIII_plot, threshold %in% "pos")
III_IIII_plot_neg <- subset(III_IIII_plot, threshold %in% "neg")

#chromPlot function
#histogram
pdf(file.path(outdir, "IvsII_chromPlot.pdf"), width = 3, height = 3.9)
chromPlot(gaps = hg_gap, 
          annot1 = II_III_plot_pos, 
          annot2 = II_III_plot_neg,
          colAnnot1 = c("#cc5f43"),
          colAnnot2 = c("#8176cc"),
          chrSide=c(-1, 1, 1, 1, 1, 1, 1, 1),
          cex = 0.75)
dev.off()

pdf(file.path(outdir, "IvsIII_chromPlot.pdf"), width = 3, height = 3.9)
chromPlot(gaps = hg_gap, 
          annot1 = II_IIII_plot_pos, 
          annot2 = II_IIII_plot_neg,
          colAnnot1 = c("#cc5f43"),
          colAnnot2 = c("#8176cc"),
          chrSide=c(-1, 1, 1, 1, 1, 1, 1, 1),
          cex = 0.75)
dev.off()

pdf(file.path(outdir, "IIvsIII_chromPlot.pdf"), width = 3, height = 3.9)
chromPlot(gaps = hg_gap, 
          annot1 = III_IIII_plot_pos, 
          annot2 = III_IIII_plot_neg,
          colAnnot1 = c("#cc5f43"),
          colAnnot2 = c("#8176cc"),
          chrSide=c(-1, 1, 1, 1, 1, 1, 1, 1),
          cex = 0.75)
dev.off()

#Chromosome 1 and 19 comparisons
Chrom <- c("1p", "1q", "19p", "19q")

II_III_plot <- subset(II_III_plot, Chrom %in% c("1", "19"))
frequency <- c(nrow(II_III_plot[II_III_plot$Chrom == "1" & II_III_plot$band == "p" & II_III_plot$threshold == "pos", ])/
                 nrow(II_III_plot[II_III_plot$Chrom == "1" & II_III_plot$band == "p" & II_III_plot$threshold == "neg", ]),
               nrow(II_III_plot[II_III_plot$Chrom == "1" & II_III_plot$band == "q" & II_III_plot$threshold == "pos", ])/
                 nrow(II_III_plot[II_III_plot$Chrom == "1" & II_III_plot$band == "q" & II_III_plot$threshold == "neg", ]),
               nrow(II_III_plot[II_III_plot$Chrom == "19" & II_III_plot$band == "p" & II_III_plot$threshold == "pos", ])/
                 nrow(II_III_plot[II_III_plot$Chrom == "19" & II_III_plot$band == "p" & II_III_plot$threshold == "neg", ]),
               nrow(II_III_plot[II_III_plot$Chrom == "19" & II_III_plot$band == "q" & II_III_plot$threshold == "pos", ])/1)
frequency <- round(frequency, 2)
II_III_freq <- as.data.frame(cbind(Chrom, frequency))
II_III_freq$Chrom <- factor(II_III_freq$Chrom, levels = c("1p", "1q", "19p", "19q"))
II_III_freq$frequency <- as.numeric(II_III_freq$frequency)

II_IIII_plot <- subset(II_IIII_plot, Chrom %in% c("1", "19"))
frequency <- c(nrow(II_IIII_plot[II_IIII_plot$Chrom == "1" & II_IIII_plot$band == "p" & II_IIII_plot$threshold == "pos", ])/
                 nrow(II_IIII_plot[II_IIII_plot$Chrom == "1" & II_IIII_plot$band == "p" & II_IIII_plot$threshold == "neg", ]),
               nrow(II_IIII_plot[II_IIII_plot$Chrom == "1" & II_IIII_plot$band == "q" & II_IIII_plot$threshold == "pos", ])/
                 nrow(II_IIII_plot[II_IIII_plot$Chrom == "1" & II_IIII_plot$band == "q" & II_IIII_plot$threshold == "neg", ]),
               nrow(II_IIII_plot[II_IIII_plot$Chrom == "19" & II_IIII_plot$band == "p" & II_IIII_plot$threshold == "pos", ])/
                 nrow(II_IIII_plot[II_IIII_plot$Chrom == "19" & II_IIII_plot$band == "p" & II_IIII_plot$threshold == "neg", ]),
               nrow(II_IIII_plot[II_IIII_plot$Chrom == "19" & II_IIII_plot$band == "q" & II_IIII_plot$threshold == "pos", ])/
                 nrow(II_IIII_plot[II_IIII_plot$Chrom == "19" & II_IIII_plot$band == "q" & II_IIII_plot$threshold == "neg", ]))
frequency <- round(frequency, 2)
II_IIII_freq <- as.data.frame(cbind(Chrom, frequency))
II_IIII_freq$Chrom <- factor(II_IIII_freq$Chrom, levels = c("1p", "1q", "19p", "19q"))
II_IIII_freq$frequency <- as.numeric(II_IIII_freq$frequency)

III_IIII_plot <- subset(III_IIII_plot, Chrom %in% c("1", "19"))
frequency <- c(nrow(III_IIII_plot[III_IIII_plot$Chrom == "1" & III_IIII_plot$band == "p" & III_IIII_plot$threshold == "pos", ])/
                 nrow(III_IIII_plot[III_IIII_plot$Chrom == "1" & III_IIII_plot$band == "p" & III_IIII_plot$threshold == "neg", ]),
               nrow(III_IIII_plot[III_IIII_plot$Chrom == "1" & III_IIII_plot$band == "q" & III_IIII_plot$threshold == "pos", ])/
                 nrow(III_IIII_plot[III_IIII_plot$Chrom == "1" & III_IIII_plot$band == "q" & III_IIII_plot$threshold == "neg", ]),
               nrow(III_IIII_plot[III_IIII_plot$Chrom == "19" & III_IIII_plot$band == "p" & III_IIII_plot$threshold == "pos", ])/
                 nrow(III_IIII_plot[III_IIII_plot$Chrom == "19" & III_IIII_plot$band == "p" & III_IIII_plot$threshold == "neg", ]),
               nrow(III_IIII_plot[III_IIII_plot$Chrom == "19" & III_IIII_plot$band == "q" & III_IIII_plot$threshold == "pos", ])/
                 nrow(II_IIII_plot[II_IIII_plot$Chrom == "19" & II_IIII_plot$band == "q" & II_IIII_plot$threshold == "neg", ]))
frequency <- round(frequency, 2)
III_IIII_freq <- as.data.frame(cbind(Chrom, frequency))
III_IIII_freq$Chrom <- factor(III_IIII_freq$Chrom, levels = c("1p", "1q", "19p", "19q"))
III_IIII_freq$frequency <- as.numeric(III_IIII_freq$frequency)

### Set theme
my_theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
                  axis.line = element_line(colour = "black"),
                  legend.position="none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.text = element_text(size = 15),
                  axis.title = element_text(size = 15))

### Plot chromosome frequencies
pdf(file.path(outdir, "IvsII_chromFreq.pdf"), width = 3, height = 5)
ggplot(I_II_freq, aes(x = Chrom, y = frequency, fill = Chrom)) +
  geom_bar(color = "Black", stat="identity", 
           width = 0.8, position = position_dodge(width = 1)) +
  geom_text(aes(label=frequency), position = position_dodge(width = 0.9), vjust = -0.25, size = 6) +
  xlab("Chromosome") + 
  ylab("DE Genes positive/negative") +
  ggtitle("") + 
  scale_fill_manual(values=c("#cc5f43", "#8176cc", "#cc5f43", "#8176cc")) +
  my_theme +
  scale_y_continuous(limits = c(0,12))
dev.off()

pdf(file.path(outdir, "IvsIII_chromFreq.pdf"), width = 3, height = 5)
ggplot(I_III_freq, aes(x = Chrom, y = frequency, fill = Chrom)) +
  geom_bar(color = "Black", stat="identity", 
           width = 0.8, position = position_dodge(width = 1)) +
  geom_text(aes(label=frequency), position = position_dodge(width = 0.9), vjust = -0.25, size = 6) +
  xlab("Chromosome") + 
  ylab("DE Genes positive/negative") +
  ggtitle("") + 
  scale_fill_manual(values=c("#cc5f43", "#8176cc", "#cc5f43", "#8176cc")) +
  my_theme +
  scale_y_continuous(limits = c(0,3))
dev.off()

pdf(file.path(outdir, "IIvsIII_chromFreq.pdf"), width = 3, height = 5)
ggplot(II_III_freq, aes(x = Chrom, y = frequency, fill = Chrom)) +
  geom_bar(color = "Black", stat="identity", 
           width = 0.8, position = position_dodge(width = 1)) +
  geom_text(aes(label=frequency), position = position_dodge(width = 0.9), vjust = -0.25, size = 6) +
  xlab("Chromosome") + 
  ylab("DE Genes positive/negative") +
  ggtitle("") + 
  scale_fill_manual(values=c("#cc5f43", "#8176cc", "#cc5f43", "#8176cc")) +
  my_theme +
  scale_y_continuous(limits = c(0,4))
dev.off()