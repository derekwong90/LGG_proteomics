library(tidyverse)
library(dplyr)
library(ggplot2)
library(nnet)
library(pROC)
library(caret)
library(MASS)

### Set variables
options(scipen = 100)
path <- "/Users/derekwong/Google Drive/Post-Doc/Yip_Proteomics/Protein_model"
proteome <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Proteome/output/protein_abundance.txt"
samples <- "/Users/derekwong/OneDrive - UHN/Post-Doc/Yip_projects/LGG_proteomics/Proteome/output/metadata.txt"

### Read in data
proteome <- read.delim(proteome)
samples <- read.delim(samples)

### Subset samples
samples <- samples[samples$Type %in% c("Oligodendroglioma", "Glioblastoma", "Astrocytoma", "Normal"), ]
samples$Type <- factor(samples$Type, levels = c("Normal", "Oligodendroglioma", "Astrocytoma", "Glioblastoma"),
                       labels = c("Normal", "Type I", "Type II", "Type III"))

### Format proteome data
proteome <- proteome[complete.cases(proteome), ]
proteome <- proteome[!(duplicated(proteome$Gene)), ]
row.names(proteome) <- proteome$Gene
proteome <- proteome[, c(10:ncol(proteome))]
proteome <- as.data.frame(t(proteome))

### Center and normalize data
scaled_data <- as.data.frame(scale(proteome))

### Merge with metadata
scaled_data <- merge(scaled_data, samples, by.x = "row.names", by.y = c("set"))

### Keep model proteins
keep <- c("FLNA","SLC4A1","CNN3","HAPLN1","IQGAP1", "MUCL1", "VIM", "H2BU1", "SPOCK3", "MFSD10", "ALDH2", "NES", "NTM", "GNPTG", "BTBD17", "CNRIP1",   
          "VAT1L", "FBLL1", "MACROH2A2", "AK3", "PHYHIPL", "ABAT", "PURA", "PURG", "SEZ6L", "ALDH7A1", "IDI1", "SLC15A2", "KCNN3", "S100A1", "FXYD7", "LY6H",     
          "TSPAN7", "SVIP", "PRNP", "GNG3", "BASP1", "C1QC", "IGHV3-7", "ALB", "IGHV1-46", "ITIH2", "IGHG3", "TIMP1", "IGHG2", "VKORC1", "IGKV3-20", "IGHV3-49", 
          "SERPINA3", "IGHV4-28", "IGHV3-64D", "IGHV1-69", "CHI3L1", "BCAT1", "Type")

scaled_data <- scaled_data[scaled_data$Type %in% c("Type I", "Type II", "Type III"), colnames(scaled_data) %in% keep]
scaled_data$Type <- factor(scaled_data$Type, levels = c("Type I", "Type II", "Type III"),
                           labels = c("Type_I", "Type_II", "Type_III"))

### Make the model
train_control <- trainControl(method='repeatedcv', number=10, repeats = 10, verboseIter = TRUE, classProbs = TRUE, savePredictions = TRUE)

model_multi <- caret::train(Type~., data = scaled_data, trControl = train_control, method='multinom')
print(model_multi)

ClassPredicted <- predict(model_multi, newdata = scaled_data, "raw")
ClassProb <- predict(model_multi, newdata = scaled_data, "prob")

### Get p-values
z <- summary(model_multi)$coefficients/summary(model_multi)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1))*2 # we are using two-tailed z test

### Get AUC values
response <- model_multi$pred$obs
response <- as.factor(response)
predictor <- cbind(model_multi$pred$Type_I, model_multi$pred$Type_II, model_multi$pred$Type_III)
colnames(predictor) <- c("Type_I", "Type_II", "Type_III")

auc <- multiclass.roc(response, predictor)

### Get curve values
# Type I vs Type II
length <- ifelse(length(auc[["rocs"]][["Type_I/Type_II"]][[1]][["sensitivities"]]) > length(auc[["rocs"]][["Type_I/Type_II"]][[2]][["sensitivities"]]),
                 length(auc[["rocs"]][["Type_I/Type_II"]][[2]][["sensitivities"]]), length(auc[["rocs"]][["Type_I/Type_II"]][[1]][["sensitivities"]]))
Type_I_sens <- cbind(auc[["rocs"]][["Type_I/Type_II"]][[1]][["sensitivities"]], auc[["rocs"]][["Type_I/Type_II"]][[2]][["sensitivities"]])
Type_I_sens <- Type_I_sens[1:length, ]
Type_I_sensitivity <- rowMeans(Type_I_sens)
Type_I_sensitivity <- 1-Type_I_sensitivity

Type_I_spec <- cbind(auc[["rocs"]][["Type_I/Type_II"]][[1]][["specificities"]], auc[["rocs"]][["Type_I/Type_II"]][[2]][["specificities"]])
Type_I_spec <- Type_I_spec[1:length, ]
Type_I_specificity <- rowMeans(Type_I_spec)

Type_I <- as.data.frame(cbind(Type_I_sensitivity, Type_I_specificity))
colnames(Type_I) <- c("sensitivity", "specificity")

auc_I <- round(sum(Type_I_specificity*diff(c(0, Type_I_sensitivity))),3)

# Type I vs Type III
length <- ifelse(length(auc[["rocs"]][["Type_I/Type_III"]][[1]][["sensitivities"]]) > length(auc[["rocs"]][["Type_I/Type_III"]][[2]][["sensitivities"]]),
                 length(auc[["rocs"]][["Type_I/Type_III"]][[2]][["sensitivities"]]), length(auc[["rocs"]][["Type_I/Type_III"]][[1]][["sensitivities"]]))
Type_II_sens <- cbind(auc[["rocs"]][["Type_I/Type_III"]][[1]][["sensitivities"]], auc[["rocs"]][["Type_I/Type_III"]][[2]][["sensitivities"]])
Type_II_sens <- Type_II_sens[1:length, ]
Type_II_sensitivity <- rowMeans(Type_II_sens)
Type_II_sensitivity <- 1-Type_II_sensitivity

Type_II_spec <- cbind(auc[["rocs"]][["Type_I/Type_III"]][[1]][["specificities"]], auc[["rocs"]][["Type_I/Type_III"]][[2]][["specificities"]])
Type_II_spec <- Type_II_spec[1:length, ]
Type_II_specificity <- rowMeans(Type_II_spec)

Type_II <- as.data.frame(cbind(Type_II_sensitivity, Type_II_specificity))
colnames(Type_II) <- c("sensitivity", "specificity")

auc_II <- round(sum(Type_II_specificity*diff(c(0, Type_II_sensitivity))),3)

# Type II vs Type III
length <- ifelse(length(auc[["rocs"]][["Type_II/Type_III"]][[1]][["sensitivities"]]) > length(auc[["rocs"]][["Type_II/Type_III"]][[2]][["sensitivities"]]),
                 length(auc[["rocs"]][["Type_II/Type_III"]][[2]][["sensitivities"]]), length(auc[["rocs"]][["Type_II/Type_III"]][[1]][["sensitivities"]]))
Type_III_sens <- cbind(auc[["rocs"]][["Type_II/Type_III"]][[1]][["sensitivities"]], auc[["rocs"]][["Type_II/Type_III"]][[2]][["sensitivities"]])
Type_III_sens <- Type_III_sens[1:length, ]
Type_III_sensitivity <- rowMeans(Type_III_sens)
Type_III_sensitivity <- 1-Type_III_sensitivity

Type_III_spec <- cbind(auc[["rocs"]][["Type_II/Type_III"]][[1]][["specificities"]], auc[["rocs"]][["Type_II/Type_III"]][[2]][["specificities"]])
Type_III_spec <- Type_III_spec[1:length, ]
Type_III_specificity <- rowMeans(Type_III_spec)

Type_III <- as.data.frame(cbind(Type_III_sensitivity, Type_III_specificity))
colnames(Type_III) <- c("sensitivity", "specificity")

auc_III <- round(sum(Type_III_specificity*diff(c(0, Type_III_sensitivity))),3)

### Mean AUC
auc_mean <- round((auc_I + auc_II + auc_III)/3, 3)

### Plot curves 
AUC_plot<- ggplot() +
  geom_line(data = Type_I, aes(x = sensitivity, y = specificity, color = "I")) +
  geom_line(data = Type_II, aes(x = sensitivity, y = specificity, color = "II")) +
  geom_line(data = Type_III, aes(x = sensitivity, y = specificity, color = "III")) +
  geom_abline(intercept = 0, linetype = "dashed", alpha = 0.5) +
  xlab("False Positive Rate") + 
  ylab("True Positive Rate") +
  labs(color = "") + 
  ggtitle("") + 
  scale_color_manual(labels = c(paste0("Type I vs Type II: AUC = ", auc_I),
                                paste0("Type I vs Type III: AUC = ", auc_II),
                                paste0("Type II vs Type III: AUC = ", auc_III)), 
                     values = c("I" = "#359B73", "II" = "#D55E00", "III" = "#2271B2")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size = 15),
      legend.position = c(0.65, 0.2),
      legend.background = element_blank(),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15)) +
  annotate("text", x = 0.65, y = 0.3, label = paste0("Model AUC: ", auc_mean), size = 6)
AUC_plot

ggsave(file.path(path, "protein_classifier.pdf"), width = 5.5, height = 5)
