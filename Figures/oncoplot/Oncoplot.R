library(ComplexHeatmap)
library(dplyr)
library(matrixStats)
library(circlize)

### Set variables
path <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Genome"
outdir <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Oncoplot"

data_onco <- file.path(path, "mutations.txt")
data_samples <- file.path(path, "subtype.txt")

### Import data
data_onco <- read.delim(data_onco, header = FALSE)
data_samples <- read.delim(data_samples)

### Format data tables
data_onco <- data_onco[, 1:3]
colnames(data_onco) <- c("Sample", "Gene", "Alteration")

names(data_samples)[names(data_samples) == 'X1p'] <- 'Chromosome 1p'
names(data_samples)[names(data_samples) == 'X19q'] <- 'Chromosome 19q'

data_samples$Molecular_Subtype <- factor(data_samples$Molecular_Subtype, levels = c("I", "II", "III"),
                                         labels = c("Type I", "Type II", "Type III"))

data_samples$Diagnosis <- factor(data_samples$Diagnosis, levels = c("Oligodendroglioma", "Astrocytoma", "Glioblastoma"))

### Remove any extraneous spaces that may be within the data.
cols_to_be_rectified <- names(data_onco)[vapply(data_onco, is.character, logical(1))]
data_onco[,cols_to_be_rectified] <- lapply(data_onco[,cols_to_be_rectified], trimws)

### Cast data to wide format
data_cast <- reshape(data_onco, timevar="Gene", idvar="Sample", direction="wide")
colnames(data_cast) <- gsub("Alteration.", "", colnames(data_cast))

### Merge samples with mutations for copy number (1p/19q)
data_samples <- merge(data_samples, data_cast, by = c("Sample"), all = TRUE)
row.names(data_samples) <- data_samples$Sample

### Order samples
data_samples <- data_samples[order(data_samples$`Chromosome 1p`,
                                   data_samples$`Chromosome 19q`,
                                   data_samples$IDH1,
                                   data_samples$IDH2,
                                   data_samples$CIC,
                                   data_samples$TP53,
                                   data_samples$ATRX,
                                   data_samples$CDKN2A,
                                   data_samples$PIK3CA,
                                   data_samples$NOTCH1,
                                   data_samples$MET,
                                   data_samples$BRAF,
                                   data_samples$FUBP1,
                                   data_samples$EGFR,
                                   data_samples$MDM2,
                                   data_samples$PIK3R1,
                                   data_samples$PTEN,
                                   data_samples$RB1,
                                   data_samples$NF1,
                                   data_samples$TERT), ]

### Split into dataframe for oncoplot
data_onco <- data_samples[, -c(1:3)]
data_samples <- data_samples[ , 1:3]

### Set clinical annotations
data_type <- as.matrix(data_samples$Molecular_Subtype)
row.names(data_type) <- data_samples$Sample

data_diag <- as.matrix(data_samples$Diagnosis)
row.names(data_diag) <- data_samples$Diagnosis

### Transpose oncoplot matrix
data_onco <- as.data.frame(t(data_onco))
data_onco <- as.matrix(data_onco)

### Set annotation and oncoplot colours
col <- c(MISSENSE_DRIVER = "#33A02C", MISSENSE = "#B2DF8A",
         TRUNC_DRIVER = "black", INDEL = "#FF7F00",
         HETLOSS = "#A6CEE3", HOMODEL = "#1F78B4", 
         GAIN = "#FB9A99", AMP = "#E31A1C")
col_type <- c("Type I" = "#5eb956", "Type II" = "#c18743", "Type III" = "#6f7dcb")
col_diag <- c(Oligodendroglioma = "#5eb956", Astrocytoma = "#c18743", Glioblastoma = "#6f7dcb")

### Set variables
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = "grey80", col = NA)),
  MISSENSE_DRIVER = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["MISSENSE_DRIVER"], col = NA)),
  MISSENSE = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["MISSENSE"], col = NA)),
  TRUNC_DRIVER = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["TRUNC_DRIVER"], col = NA)),
  INDEL = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["INDEL"], col = NA)),
  HETLOSS = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["HETLOSS"], col = NA)),
  HOMODEL = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["HOMODEL"], col = NA)),
  GAIN = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["GAIN"], col = NA)),
  AMP = function(x, y, w, h)
    grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = col["AMP"], col = NA)))

### Compile annotation layers
top_annotation <- HeatmapAnnotation("Molecular Subtype" = data_type,
                                    "Histology Subtype" = data_diag,
                                    annotation_name_side = "left",
                                    annotation_name_gp = gpar(fontsize = 13),
                                    height = unit(1.5, "cm"),
                                    simple_anno_size = unit(0.5, "cm"),
                                    border = TRUE,
                                    annotation_name_rot = 0,
                                    col = list("Molecular Subtype" = col_type, "Histology Subtype" = col_diag))

### Set labels and legends
heatmap_legend_param = list(title = "Alternations", 
                            at = c("MISSENSE_DRIVER", "MISSENSE", "TRUNC_DRIVER", "INDEL", "HETLOSS", "HOMODEL", "GAIN", "AMP"), 
                            labels = c("Missense Driver", "Missense", "Truncating", "Insertion/Deletion", 
                                       "Loss of Heterozygosity", "Homozygous Deletion", "Gain", "Amplification"),
                            nrow = 4)

### Set column/row order and splits
sample_order <- colnames(data_onco)
gene_order <- c("Chromosome 1p", "Chromosome 19q", "IDH1", "IDH2", "CIC", "TP53", "ATRX", "CDKN2A", "PIK3CA", 
                "NOTCH1", "MET", "BRAF", "FUBP1", "EGFR", "MDM2", "PIK3R1", "PTEN", "RB1", "NF1", "TERT")
column_split <- data_samples$Molecular_Subtype
column_split <- factor(column_split, levels = c("Type I", "Type II", "Type III"))

### Generate oncoprint
pdf(file.path(outdir, "Oncoplot.pdf"), width = 8, height = 6)
oncoPrint <- oncoPrint(data_onco,
                       alter_fun = alter_fun, 
                       col = col,
                       row_names_side = "left", pct_side = "right",
                       top_annotation = top_annotation,
                       column_title = NULL,
                       heatmap_legend_param = heatmap_legend_param,
                       row_order = gene_order,
                       column_order = sample_order,
                       row_names_gp = gpar(fontsize = 13),
                       column_split = column_split,
                       border = TRUE,
                       border_gp = gpar(col = "black"))
draw(oncoPrint, heatmap_legend_side = "bottom", merge_legend = TRUE)
dev.off()
