library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Transcriptome"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/TCGA_heatmap"

dat <- read.delim(file.path(path, "LGG.uncv2.mRNAseq_RSEM_all.txt"))
master <- read.delim(file.path(path, "LGG-GBM Cohort.txt"))

### Format metadata
master <- master[!duplicated(master), ]
row.names(master) <- master$TCGA.ID
master <- master[, 2:7]
master$sample <- row.names(master)

master$IDH <- factor(master$IDH, levels = c("Mut", "WT"),
                     labels = c("Mutant", "Wildtype"))
master$X1p19q <- factor(master$X1p19q, levels = c("LOH", "ROH"),
                        labels = c("Loss of Heterozygosity", "Retained"))
master$CIC <- factor(master$CIC, levels = c("LOF", "Missense", "WT"),
                     labels = c("Loss of Function", "Missense", "Wildtype"))
master$Type <- factor(master$Type, levels = c("I", "II", "III"),
                      labels = c("IDH mutant, 1p/19q LOH", "IDH mutant, 1p/19q ROH", "IDH wildtype"))

### Use this to plot only ODGs
master <- master[master$Type == c("IDH mutant, 1p/19q LOH"), ]

### Format expression matrix
row.names(dat) <- dat$HYBRIDIZATION.R
dat$Gene.ID <- strsplit(as.character(row.names(dat)), split="[|]")
dat$Gene <- sapply(dat$Gene.ID, "[[", 1)
col_idx <- grep("Gene", names(dat))
dat <- dat[, c(col_idx, (1:ncol(dat))[-col_idx])]
dat <- dat[, -c(1)]
dat <- dat[!duplicated(dat$Gene),]
dat <- dat[!(dat$Gene=="?"),]
row.names(dat) <- dat$Gene

dat_heatmap <- dat[,3:ncol(dat)]
dat_heatmap <- dat_heatmap[ ,colnames(dat_heatmap) %in% rownames(master)]
dat_heatmap[dat_heatmap == 0] <- NA
dat_heatmap <- dat_heatmap[complete.cases(dat_heatmap), ]
dat_heatmap <- as.data.frame(dat_heatmap)

master <- master[row.names(master) %in% colnames(dat_heatmap), ]

### Order expression matrix
dat_heatmap <- dat_heatmap[, row.names(master)]

#Log transform data
dat_log <- as.data.frame(lapply(dat_heatmap, function(x) log2(as.numeric(as.character(x))) ))
row.names(dat_log) <- row.names(dat_heatmap)
dat_trans_log <- t(dat_log)

#scale data average = 0 variant =1
sprDat <- t(scale(t(dat_log))) %>% as.data.frame()
str(sprDat, max.level = 0, give.attr = FALSE)

round(data.frame(avgBefore = rowMeans(head(dat_log)), avgAfter = rowMeans(head(sprDat)),
                 varBefore = apply(head(dat_log), 1, var), varAfter = apply(head(sprDat), 1, var)),
      2)
dat_to_plot <- sprDat
dat_clip <- sprDat %>% as.data.frame()

#Take the top 1000 most variable by MAD
mads = apply(dat_clip,1,mad)
dat_top <- dat_clip[rev(order(mads))[1:500],] %>% as.data.frame()

#set >5 to 5, and <-5 to -5 for clipped heatmap
dat_top[] <- lapply(dat_top, function(x) ifelse(x>5,5,x))
dat_top[] <- lapply(dat_top, function(x) ifelse(x< -5,-5,x))

### Convert to matrix
dat_top <- as.matrix(dat_top)

# Set annotations
data_IDH <- as.matrix(master$IDH)
row.names(data_IDH) <- row.names(master)

data_codel <- as.matrix(master$X1p19q)
row.names(data_codel) <- row.names(master)

data_CIC <- as.matrix(master$CIC)
row.names(data_CIC) <- row.names(master)

data_type <- as.matrix(master$Type)
row.names(data_type) <- row.names(master)

#Set Annotation colors
col_fun <- colorRamp2(c(-3, 0, 3), 
                      c("#313695", "white", "#A50026"))
col_IDH <- c(Mutant = "#cc573f", Wildtype = "#b0b043")
col_codel <- c("Loss of Heterozygosity" = "#4fb7a7", "Retained" = "#c06e94")
col_CIC <- c("Loss of Function" = "#cc573f", Missense = "#5c7e3b", Wildtype = "#b0b043")
col_type <- c("IDH mutant, 1p/19q LOH" = "#5eb956", "IDH mutant, 1p/19q ROH" = "#c18743", "IDH wildtype" = "#6f7dcb")

### Set Annotations
top_annotation <- HeatmapAnnotation(IDH = data_IDH,
                                    "1p/19q" = data_codel,
                                    "Molecular Subtype" = data_type,
                                    CIC = data_CIC,
                                    col = list(IDH = col_IDH,
                                               "1p/19q" = col_codel,
                                               "Molecular Subtype" = col_type,
                                               CIC = col_CIC),
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 13),
                                                                   labels_gp = gpar(fontsize = 11)),
                                    annotation_name_gp = gpar(fontsize = 15),
                                    border = TRUE,
                                    simple_anno_size = unit(0.5, "cm"))

### Set Legend Parameters
heatmap_legend_param = list(title = "Log2Expression", 
                            border = TRUE,
                            at = c(-3, 0, 3),
                            title_gp = gpar(fontsize = 13),
                            labels_gp = gpar(fontsize = 11),
                            direction = "horizontal",
                            legend_width = unit(4, "cm"))

### Plot and save heatmap
pdf(file.path(outdir, "ODG_Heatmap.pdf"), width = 8, height = 8)
heatmap <- Heatmap(dat_top,
                   col = col_fun,
                   top_annotation = top_annotation,
                   heatmap_legend_param = heatmap_legend_param,
                   show_column_dend = TRUE,
                   show_column_names = FALSE,
                   show_row_names = FALSE,
                   row_title_rot = 0,
                   row_title_gp = gpar(fontsize = 15),
                   border = TRUE)
draw(heatmap, heatmap_legend_side = "right", merge_legend = TRUE)
dev.off()
