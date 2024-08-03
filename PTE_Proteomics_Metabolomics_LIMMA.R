# This is to perform single-omics analysis for the proteomics and metabolomics data using limma

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading Libraries
library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel)

#----Proteomics----
# Loading in dataframe
dt <- read.csv("PTE_Proteomics_Pre-processed_LIMMA.csv", sep = ",", header = T, 
               stringsAsFactors = F, row.names = 1)
dt <- dt %>% select(order(colnames(.))) #Re-ordering
# Loading experimental design file
exp_des <- read.csv("PTE_OriginalMatrix_ExpDes.csv", stringsAsFactors = F, sep = ",", header = 1)
exp_des <- exp_des %>%
  filter(!grepl("Ref",label)) %>% # Removing reference channels
  unite("new_label", 2:3, sep = "_", remove = F) %>% # Setting new sample labels
  select(!label) %>%
  rename(label = new_label) %>%
  arrange(label)
exp_des$condition  <- gsub("Sha","Sham", exp_des$condition) # Changing 'Sha' to 'Sham'

# Removing outliers from whole dataset - Naive 1, PTE 7, Sham 10
dt_noOut <- dt[,-c(1,18,22)] 
exp_noOut <- exp_des[-c(1,18,22),] 

# Design matrix and contrast
exp_noOut$condition <- as.factor(exp_noOut$condition)
exp_noOut$batch <- as.factor(exp_noOut$batch)

dt_des <- model.matrix(~0 + condition + batch, exp_noOut)
colnames(dt_des) <- gsub("condition","", colnames(dt_des))
colnames(dt_des) <- gsub(" ","", colnames(dt_des))

cont_anova <- makeContrasts(
  # Contrasts of interest
  `Sham vs Naïve` = "Sham - Naïve", # effects of craniotomy 
  `TBI-PTE vs Naïve` = "TBI - Naïve", # protein changes in TBI
  `TBI+PTE vs Naïve` = "PTE - Naïve", # protein changes in PTE
  `TBI+PTE vs TBI-PTE` = "PTE - TBI", # protein changes in PTE post-TBI
  
  # Levels
  levels = dt_des
)

# Run Limma
# Fit the expression matrix to a linear model
fit <- lmFit(dt_noOut, dt_des)
# Compute contrast
fit_contrast <- contrasts.fit(fit, cont_anova)
# Bayes statistics of differential expression
fit_bayes <- eBayes(fit_contrast, trend = T)
# Plot the mean-variance relationship to see if there is a trend in the expression of the features
pdf("PTE_Proteomics_MeanVar_NoOutliers.pdf")
plotSA(fit_bayes) 
title("Mean-Variance Relationship")
dev.off()

# Summary of number of differentially expressed genes (across all contrasts)
result <- decideTests(fit_bayes, method = "separate", adjust.method = "BH", 
                      p.value = 0.05, lfc = 0.585)
summary(result)
res_df <- as.data.frame(summary(result))
colnames(res_df) <- c("Direction", "Pair", "Number")
write.csv(res_df, "PTE_Proteomics_NoOutliers_FC1.5_DE Result Summary.csv", row.names = F)

# Getting a topTable for each individual contrast
top_list <- vector('list', length(colnames(cont_anova)))

for (i in 1:length(colnames(cont_anova))) {
  top_list[[i]] <- topTable(fit_bayes, n=Inf, adjust.method = "BH", coef = i)
  names(top_list) <- colnames(result)
}
# Adding information about up- or down-regulation
top_mod <- map(top_list,
               ~ mutate(.x, "Sig DE" = case_when(logFC >= log2(1.5) & adj.P.Val <= 0.05 ~ "Up",
                                                 logFC <= -log2(1.5) & adj.P.Val <= 0.05 ~ "Down",
                                                 TRUE ~ "Not Sig")))
# Generating file for each individual table
iwalk(top_mod, ~write.csv(.x, file.path("PTE Proteomics LIMMA", paste0(.y, "_FC1.5_NoOutliers_TBI.B.csv"))))

# Make volcano plots for each of the individual comparisons
pdf("PTE_Proteomics_Volcano Plots_NoOutliers_FC1.5.pdf")
for (i in seq(length(top_mod))) {
  colors <- c(Up = "firebrick3", Down = "dodgerblue3", `Not Sig` = "gray50")
  p <- ggplot(top_mod[[i]], aes(logFC, -log(P.Value,10))) +
    geom_point(aes(color = `Sig DE`), size = 1) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"(P-Value)")) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    geom_text_repel(data = top_mod[[i]],
                    aes(logFC, -log(P.Value,10), 
                        label = ifelse(top_mod[[i]]$`Sig DE` !="Not Sig", 
                                       rownames(top_mod[[i]]), "")),
                    size = 3,
                    max.overlaps = 50) +
    ggtitle(names(top_mod)[[i]])
  print(p)
}
dev.off() 

#----Metabolomics----
# Reading in files
dt <- read.csv("PTE_Metabolomics_Pre-processed_SampImp_LIMMA.csv", sep = ",", header = T, 
               stringsAsFactors = F, row.names = 1)
colnames(dt)  <- gsub("\\."," ", colnames(dt))
# Loading experimental design file
exp_des <- read.csv("PTE_Metabolomics_ExpDes.csv", header = T, stringsAsFactors = F)
exp_des <- exp_des %>%
  select(!sample) %>% # Removing sample group
  filter(group != "PQC") %>% # Removing quality control samples
  arrange(label, colnames(dt))

# Removing outliers - Naive 1, Sham 5 (SampImp)
dt_noOut <- dt[,-c(1,28)]
exp_noOut <- exp_des[-c(1,28),]

# Design matrix and contrast
exp_noOut$group <- as.factor(exp_noOut$group)
dt_des <- model.matrix(~0 + group, exp_noOut)
colnames(dt_des) <- gsub("group","", colnames(dt_des))

cont_anova <- makeContrasts(
  # Contrasts of interest
  `Sham vs Naïve` = "Sham - Naïve", # effects of craniotomy 
  `TBI-PTE_A vs Naïve` = "TBI - Naïve", # protein changes in TBI
  `TBI+PTE vs Naïve` = "PTE - Naïve", # protein changes in PTE
  `TBI+PTE vs TBI-PTE_A` = "PTE - TBI", # protein changes in PTE post-TBI
  `TBI-PTE_A vs TBI-PTE_B` = "TBI - TBI.B",
  
  # Levels
  levels = dt_des
)

# Run Limma
# Fit the expression matrix to a linear model
fit <- lmFit(dt_noOut, dt_des)
# Compute contrast
fit_contrast <- contrasts.fit(fit, cont_anova)
# Bayes statistics of differential expression
fit_bayes <- eBayes(fit_contrast, trend = T)
# Plot the mean-variance relationship to see if there is a trend in the expression of the features
pdf("PTE_Metabolomics_MeanVar_NoOutliers.pdf")
plotSA(fit_bayes) 
title("Mean-Variance Relationship")
dev.off()

# Summary of number of differentially expressed genes
result <- decideTests(fit_bayes, method = "separate", adjust.method = "BH", 
                      p.value = 0.05, lfc = 0.585)
summary(result)
res_df <- as.data.frame(summary(result))
colnames(res_df) <- c("Direction", "Pair", "Number")
write.csv(res_df, "PTE_Metabolomics_NoOutliers_FC1.5_DE Result Summary.csv", row.names = F)

# Getting a topTable for each individual contrast
top_list <- vector('list', length(colnames(cont_anova)))

for (i in 1:length(colnames(cont_anova))) {
  top_list[[i]] <- topTable(fit_bayes, n=Inf, adjust.method = "BH", coef = i)
  names(top_list) <- colnames(result)
}

# Adding information about up- or down-regulation
top_mod <- map(top_list,
               ~ mutate(.x, "Sig DE" = case_when(logFC >= log2(1.5) & adj.P.Val <= 0.05 ~ "Up",
                                                 logFC <= -log2(1.5) & adj.P.Val <= 0.05 ~ "Down",
                                                 TRUE ~ "Not Sig")))
# Generating file for each individual table
iwalk(top_mod, ~write.csv(.x, file.path("PTE Metabolomics LIMMA", paste0(.y, "_FC1.5_NoOutliers_TBI.B.csv"))))

# Volcano plots for each of the individual comparisons
pdf("PTE_Metabolomics_Volcano Plots_NoOutliers_FC1.5.pdf")
for (i in seq(length(top_mod))) {
  colors <- c(Up = "firebrick3", Down = "dodgerblue3", `Not Sig` = "gray50")
  p <- ggplot(top_mod[[i]], aes(logFC, -log(P.Value,10))) +
    geom_point(aes(color = `Sig DE`), size = 1) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"(P-Value)")) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    geom_text_repel(data = top_mod[[i]],
                    aes(logFC, -log(P.Value,10), 
                        label = ifelse(top_mod[[i]]$`Sig DE` !="Not Sig", 
                                       rownames(top_mod[[i]]), "")),
                    size = 3,
                    max.overlaps = Inf) +
    ggtitle(names(top_mod)[[i]])
  print(p)
}
dev.off() 
```