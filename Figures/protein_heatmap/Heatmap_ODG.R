library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Proteome/output"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Protein_heatmap"

dat <- read.delim(file.path(path, "protein_abundance.txt"))
master <- read.delim(file.path(path, "metadata.txt"))

### Format metadata
master <- master[!(master$Type %in% c("Supermix", "PIS")), ]
row.names(master) <- master$set

master$cic_mutation <- factor(master$cic_mutation, levels = c("LOF", "Missense", "WT", ""),
                     labels = c("Loss of Function", "Missense", "Wildtype", "Wildtype"))
master$Type <- factor(master$Type, levels = c("Normal", "Oligodendroglioma", "Astrocytoma", "Glioblastoma"),
                      labels = c("Normal Brain", "Type I", "Type II", "Type III"))

### Use this to plot only ODGs
master <- master[master$Type == c("Type I"), ]

### Format expression matrix
row.names(dat) <- dat$Index
dat_heatmap <- dat[complete.cases(dat), 10:ncol(dat)]

### Order expression matrix
master <- master[row.names(master) %in% colnames(dat_heatmap), ]
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
dat_top <- dat_clip[rev(order(mads))[1:5000],] %>% as.data.frame()

#set >5 to 5, and <-5 to -5 for clipped heatmap
dat_top[] <- lapply(dat_top, function(x) ifelse(x>5,5,x))
dat_top[] <- lapply(dat_top, function(x) ifelse(x< -5,-5,x))

### Convert to matrix
dat_top <- as.matrix(dat_top)

# Set annotations
data_CIC <- as.matrix(master$cic_mutation)
row.names(data_CIC) <- row.names(master)

data_exp <- as.matrix(dat_clip[row.names(dat_clip) %in% c("Q96RK0"), ])
data_exp <- (data_exp-min(data_exp))/(max(data_exp)-min(data_exp))
data_exp <- t(data_exp)
colnames(data_exp) <- "CIC Protein Expression"

data_type <- as.matrix(master$Type)
row.names(data_type) <- row.names(master)

#Set Annotation colors
col_fun <- colorRamp2(c(-3, 0, 3), 
                      c("#313695", "white", "#A50026"))
col_CIC <- c("Loss of Function" = "#cc573f", Missense = "#5c7e3b", Wildtype = "grey")
col_type <- c("Normal Brain" = "grey", "Type I" = "#5eb956", "Type II" = "#c18743", "Type III" = "#6f7dcb")
col_exp <- colorRamp2(c(0, 1), 
                      c("white", "#A50026"))

### Set Annotations
top_annotation <- HeatmapAnnotation("Molecular Subtype" = data_type,
                                    "CIC Mutation Status" = data_CIC,
                                    "CIC Protein Expression" = data_exp,
                                    col = list("Molecular Subtype" = col_type,
                                               "CIC Mutation Status" = col_CIC,
                                               "CIC Protein Expression" = col_exp),
                                    annotation_legend_param = list(title_gp = gpar(fontsize = 13),
                                                                   labels_gp = gpar(fontsize = 11)),
                                    annotation_name_gp = gpar(fontsize = 13),
                                    border = TRUE,
                                    simple_anno_size = unit(0.5, "cm"))

### Set Legend Parameters
heatmap_legend_param = list(title = "Log2Protein\nAbundance", 
                            border = TRUE,
                            at = c(-3, 0, 3),
                            title_gp = gpar(fontsize = 13),
                            labels_gp = gpar(fontsize = 11),
                            direction = "horizontal",
                            legend_width = unit(4, "cm"))

### Plot and save heatmap
pdf(file.path(outdir, "Heatmap_ODG.pdf"), width = 8, height = 6)
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
