# This is to run the pre-processing and clustering analysis of Metabolomics data.
# Again, this will be performed in a similar fashion to the GAERS experiment.

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

#----Preprocessing----
# Samples should be in columns and features in rows

# Loading libraries
library(qmtools)
library(tidyverse)
library(ggplot2)
library(naniar)
library(DEP)

# Loading file
df_init <- read.csv("P22_0549_Exp_2_data_for_MA.csv", row.names = 1) # For analysis
exp_des <- read.csv("PTE_Metabolomics_ExpDes.csv", header = T, stringsAsFactors = F) # Experimental Design

# Changing the names of the initial file
rownames(df_init) <- exp_des$label
# Setting up the file for pre-processing
miss_check <- df_init %>%
  select(!label) %>%
  arrange(rownames(.)) %>%
  t(.) %>%
  as.data.frame(.) %>%
  select(!starts_with("PQC"))

# Making SummarisedExperiment class
data <- SummarizedExperiment(assays = as.matrix(log2(miss_check)),
                             rowData = data.frame(metabolites = rownames(miss_check)),
                             colData = arrange(exp_des[-c(11:18),c(2:3)], label))
assayNames(data) <- "raw"
colData(data) # Check that column label matches samples
data
grp_info <- factor(colData(data)$group) # Sample group information that we want to plot

# Missing value plot
visdat::vis_miss(miss_check) + 
  theme(axis.text.x = element_text(angle = 90)) # add ggplot2 aesthetics to the plot 
# Density plot for missing values
p <- plot_detect(data) 

# Check filtering dimensions
dim(removeFeatures(data, i = "raw", method = "missing", cut = 0.5, group = colData(data)$group))
se <- removeFeatures(data, i = "raw",
                     method = "missing", cut = 0.5, 
                     group = colData(data)$group)
# Imputation
se <- imputeIntensity(se, i = "raw",
                      method = "QRILC",
                      sigma = 1.8,
                      MARGIN = 2L,
                      name = "QLC") 
plotBox(se, i = "QLC", group = grp_info, log2 = F) + ggtitle("QLC")
# Unlogging for VSN
dt_imp <- as.data.frame(assay(se, "QLC"))
dt_imp <- 2^dt_imp
dt <- SummarizedExperiment(assays = as.matrix(dt_imp),
                           rowData = data.frame(metabolites = rownames(dt_imp)),
                           colData = arrange(exp_des[-c(11:18),c(2:3)], label))
assayNames(dt) <- "QLC"
plotBox(dt, i = "QLC", group = grp_info, log2 = T) + ggtitle("QLC")
# Normalisation
dt <- normalizeIntensity(dt, i = "QLC", 
                         method = "vsn",
                         name = "QLC_VSN")
plotBox(dt, i = "QLC_VSN", group = grp_info) + ggtitle("QLC_VSN")
# Pareto scaling
dt <- normalizeIntensity(dt, i = "QLC_VSN",
                         method = "feature.scale",
                         type = "pareto",
                         name = "QLC_VSN_PS") 
plotBox(dt, i = "QLC_VSN_PS", group = grp_info) + ggtitle("QLC_VSN_PS")

# Clustering
m_pca <- reduceFeatures(dt, i = "QLC_VSN", method = "pca", ncomp = 5)
summary(m_pca)
plotReduced(m_pca, group = grp_info, label = T, ellipse = T) + ggtitle("QLC_VSN")
# With Pareto Scaling
pca_PS <- reduceFeatures(dt, i = "QLC_VSN_PS", method = "pca", ncomp = 5)
summary(pca_PS)
plotReduced(pca_PS, group = grp_info, label = T, ellipse = T) + ggtitle("QLC_VSN_PS")

# Saving as an RData file
saveRDS(dt, "PTE_Metabolomics_QLC_VSN_PS_SampImp_20240101.rds")
# dt <- readRDS("PTE_Metabolomics_QLC_VSN_PS_SampImp_20240101.rds") #<-- to load

# Exporting for LIMMA
dt_limma <- as.data.frame(assay(dt, "QLC_VSN_PS"))
# Transposing for PCA
data_pca <- dt_limma %>%
  t() %>%
  as.data.frame() %>%
  add_column(label = exp_des[-c(11:18),3], .before = 1)
# Saving
write.csv(dt_limma, "PTE_Metabolomics_QLC_VSN_PS_LIMMA.csv", row.names = T)
write.csv(data_pca, "PTE_Metabolomics_QLC_VSN_PS_PCA.csv", row.names = T)

# Saving plots
pdf("PTE_Metabolomics_Pre-processing_QLC_VSN_PS_.pdf")
# Initial sample distribution
plotReduced((reduceFeatures(se, i = "raw", method = "pca")), group = grp_info) + ggtitle("Raw Intensities")
# Missing values
visdat::vis_miss(miss_check) + 
  theme(axis.text.x = element_text(angle = 90))
plot(p)
# Normalised
plotBox(dt, i = "QLC", group = grp_info, log2 = T) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_discrete(name = "Group", labels = c("Naive", "PTE", "Sham", "TBI")) +
  ggtitle("QLC Feature Imputation")
# after normalization
plotBox(dt, i = "QLC_VSN", group = grp_info, log2 = F)  +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_discrete(name = "Group", labels = c("Naive", "PTE", "Sham", "TBI")) +
  ggtitle("QLC VSN")
# Scaling
plotBox(dt, i = "QLC_VSN_PS", group = grp_info, log2 = F) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_discrete(name = "Group", labels = c("Naive", "PTE", "Sham", "TBI")) +
  ggtitle("Pareto Scaled")
# Final sample clustering
plotReduced((reduceFeatures(dt, i = "QLC_VSN_PS", method = "pca")), group = grp_info, label = F, ellipse = F) + ggtitle("QLC_VSN_PS")
plotReduced((reduceFeatures(dt, i = "QLC_VSN_PS", method = "pca")), group = grp_info, label = T, ellipse = F) + ggtitle("QLC_VSN_PS")
plotReduced((reduceFeatures(dt, i = "QLC_VSN_PS", method = "pca")), group = grp_info, label = T, ellipse = T) + ggtitle("QLC_VSN_PS")
plotReduced((reduceFeatures(dt, i = "QLC_VSN", method = "pca")), group = grp_info, label = T, ellipse = F) + ggtitle("QLC_VSN")
dev.off()

#---- Clustering----
library(ggfortify)
library(plotly)
library(ggplot2)
library(tidyverse)
library(ggrepel)

dt <- read.csv("PTE_Metabolomics_Pre-processed_SampImp_PCA.csv", stringsAsFactors = F, header = T)

# Removing Naives and Sham
dt_PTE_TBI <- dt %>%
  filter(label == "PTE" | label == "TBI")
# Removing outliers
dt_noOut <- dt[-c(1,28),]

# PCA
pca_res <- prcomp(dt_noOut[,-1])
pdf("PTE_Metabolomics_QLC_SampImp_VSN_PS_Clustering.pdf")
# Plot PCA
autoplot(pca_res, data = dt) +
  geom_point(aes(color = label)) + 
  stat_ellipse(aes(color = label)) +
  geom_text_repel(label = rownames(dt), max.overlaps = 35, force = 10) +
  ggtitle("QLC_SampImp_VSN_PS")
# Scree Plot
# Finding the total variance explained by each PC
var_exp <- pca_res$sdev^2/sum(pca_res$sdev^2)
# Creating Scree Plot
qplot(c(1:(length(var_exp))), var_exp) +
  geom_line() +
  xlab("Principal Component") + 
  ylab("Variance Explained") +
  ggtitle("") +
  ylim(0, 0.20)
# Finding exact percentage total variance explained by each principal component
print(var_exp)
dev.off()

# MDS
# Setting as matrix
data_raw<-as.matrix(dt[,-1])
# Clustering Analysis
distance <- dist(data_raw) 
mds1<-cmdscale(distance, k=2) 
plot(mds1, type='n')
text(mds1, labels=rownames(data_raw), cex=0.6, adj=0.5)

# Dendrogram
dist <- dist(as.matrix(dt[,-1]), method = "euclidian", diag = T) # Calculating the euclidian distance between samples
hc <- hclust(dist) # Performing hierarchical clustering
plot(hc)