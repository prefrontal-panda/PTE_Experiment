# This is to perform the pre-processing and clustering of the proteomics data
# This will be performed in a somewhat similar fashion to the GAERS experiment with the exception of batch effect handling.
# In this experiment, the batch effects will be accounted for during limma analysis rather than before (the only exception being to make PCAs)

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

#----Preprocessing----
# Loading libraries
library(DEP)
library(tidyverse)
library(limma)
library(ggplot2)
library(ggfortify)

# Load in dataframe
prot_init <- read.csv("Imputed_matrix.csv", stringsAsFactors = F, sep = ",") # dataset for analysis
prot_id <- read.csv("Gene_ProtID.csv", stringsAsFactors = F, sep = ",") # Gene and protein ID data
prot_id <- prot_id[order(match(prot_id$Gene.Name, prot_init$ProteinID)),] # rearranging in the order of the protein table

# Load experimental design file
exp_des <- read.csv("PTE_OriginalMatrix_ExpDes.csv", stringsAsFactors = F, sep = ",", header = 1)
exp_des$condition <- gsub("Sha","Sham", exp_des$condition) # Changing 'Sha' to 'Sham'
# Changing the protein table column names
new_labels <- unite(exp_des, "new_label", 2:4, sep = "_", remove = T)
colnames(prot_init)[-1] <- new_labels$new_label
# Changing the initial column names. Include ID, replicate no, and batch
new_exp_des <- exp_des %>%
  unite("new_label", 2:4, sep = "_", remove = F) %>%
  select(!label) %>%
  rename(label = new_label) %>%
  filter(!grepl("Ref", label))

# Adding the Gene Name and Protein ID columns
prot_init <- prot_init %>%
  add_column(`Gene Name` = prot_id$Gene.Name, .before = 1) %>%
  add_column(`Protein ID`= prot_id$Protein.IDs, .before = 2) %>%
  select(!ProteinID)

# IRS normalisation
# Separation out reference channels
refs <- select(prot_init, starts_with("Ref"))
# Finding geometric means while also ignoring NAs
refs$geomean <- apply(refs, 1, function(x) exp(mean(log(x))))
# Scaling factors
refs$GRef1 <- refs$geomean / refs$Ref_1_A
refs$GRef2 <- refs$geomean / refs$Ref_2_B
refs$GRef3 <- refs$geomean / refs$Ref_3_C
# Applying scaling factors
prot_irs <- select(prot_init, ends_with("_A")) * refs$GRef1
prot_irs <- cbind(prot_irs, select(prot_init, ends_with("_B")) * refs$GRef2)
prot_irs <- cbind(prot_irs, select(prot_init, ends_with("_C")) * refs$GRef3)
# Adding the Protein ID column
prot_irs <- prot_irs %>% 
  add_column(`Gene Name` = prot_init$`Gene Name`, .before = 1) %>%
  add_column(`Protein ID`= prot_init$`Protein ID`, .before = 2) %>%
  select(-starts_with("Ref"))

# Getting ready for DEP
# Making unique names for proteins with no Gene names - use the protein ID
data_unique <- make_unique(prot_init, "Gene Name", "Protein ID", delim = ";")
# Getting the column number of the samples
samp_cols <- c(grep("Sham", colnames(data_unique)),
               grep("Naïve", colnames(data_unique)),
               grep("PTE", colnames(data_unique)),
               grep("TBI", colnames(data_unique)))
# Making the final SummarisedExperiment dataframe
data_se <- make_se(data_unique, samp_cols, new_exp_des)
#Plotting boxplot
plot_normalization(data_se)
plot_pca(data_se, x = 1, y = 2, indicate =  "condition") + ggtitle("IRS Normalisation")
plot_pca(data_se, x = 1, y = 2, indicate =  "batch") + ggtitle("IRS Normalisation")
# VSN
data_vsn <- normalize_vsn(data_se)
plot_normalization(data_vsn) + ggtitle("VSN")
plot_pca(data_vsn, x = 1, y = 2, indicate =  "condition") + ggtitle("VSN")
plot_pca(data_vsn, x = 1, y = 2, indicate =  "batch") + ggtitle("VSN")

# Convert it back to a dataframe for future analysis use.
data_prepro <- data_vsn %>%
  get_df_wide(.) %>%
  select(!Gene.Name:ID)

# Batch Correction for visualisation
# Limma's batch correction (This should only be used for visualisation purposes)
lim_batch <- as.data.frame(removeBatchEffect(x = as.matrix(data_prepro[2:46]), 
                                             batch = new_exp_des$batch,
                                             group = new_exp_des$condition))

# Final data frame
data_final <- lim_batch %>%
  add_column(`Gene Name` = prot_id$Gene.Name, .before = 1) %>%
  add_column(`Protein ID`= prot_id$Protein.IDs, .before = 2)

# Final PCA
data_fin_pca <- data_final %>%
  column_to_rownames("Protein ID") %>%
  select(-c("Gene Name")) %>%
  t(.) %>%
  as.data.frame(.) %>%
  add_column(condition = new_exp_des$condition, .before = 1) %>%
  add_column(batch = new_exp_des$batch, .after = 1)
pca_res <- prcomp(data_fin_pca[,-c(1:2)])
autoplot(pca_res, data = data_fin_pca, colour = "condition", size = 3) + 
  ggtitle("LIMMA Corrected")
boxplot(data_final[-c(1:2)])

# Save an .RData object for future use
save(data_se, data_imp, data_vsn, lim_batch, combat_batch, new_exp_des, 
     data_fin_pca, data_fin_pca_comb, pca_res, pca_res_comb,
     file = "PTE_Proteomics_Pre-processing Files_PreImpute Data.RData")
load("PTE_Proteomics_Pre-processing Files.RData")
write.csv(data_fin_pca, "PTE_Proteomics_Pre-processed_PCA.csv", row.names = T)
write.csv(data_final, "PTE_Proteomics_Pre-processed_MOFA.csv", row.names = F)
write.csv(data_prepro, "PTE_Proteomics_Pre-processed_LIMMA.csv", row.names = F)

# Saving the plots
pdf("PTE_Proteomics_Pre-processing Plots_PreImpute Data.pdf")
plot_frequency(data_se)
plot_numbers(data_se)
plot_normalization(data_se) + ggtitle("IRS Normalisation")
plot_pca(data_se, x = 1, y = 2, indicate =  "condition") + ggtitle("IRS Normalisation")
plot_pca(data_se, x = 1, y = 2, indicate =  "batch") + ggtitle("IRS Normalisation")
plot_normalization(data_vsn) + ggtitle("VSN")
plot_pca(data_vsn, x = 1, y = 2, indicate =  "condition") + ggtitle("VSN")
plot_pca(data_vsn, x = 1, y = 2, indicate =  "batch") + ggtitle("VSN")
autoplot(pca_res, data = data_fin_pca, colour = "condition", size = 3) + ggtitle("LIMMA Batch Correction")
autoplot(pca_res, data = data_fin_pca, colour = "batch", size = 3) + ggtitle("LIMMA Batch Correction")
autoplot(pca_res_comb, data = data_fin_pca_comb, colour = "condition", size = 3) + ggtitle("ComBat Corrected")
autoplot(pca_res_comb, data = data_fin_pca_comb, colour = "batch", size = 3) + ggtitle("ComBat Corrected")
dev.off()

#----Clustering Analysis----
# Loading library
library(ggfortify)
library(plotly)
library(ggplot2)
library(tidyverse)
library(ggrepel)

# Loading data
dt <- read.csv("PTE_Proteomics_Pre-processed_PCA.csv", stringsAsFactors = F, header = T,
               row.names = 1)

# Separating out the groups of interest
dt_noOut <- dt %>%
  filter(condition == "TBI" | condition == "PTE")
# No outliers
out <- c("Naïve_1", "Sham_10", "PTE_7")
dt_noOut <- dt[!row.names(dt) %in% out,]

# PCA
pca_res <- prcomp(dt_noOut[,-c(1:2)])
# Opening a pdf file to save the plot
pdf("PTE_Proteomics.pdf")
autoplot(pca_res, data = dt_noOut) + 
  geom_point(aes(color = condition)) + 
  stat_ellipse(aes(color = condition)) +
  geom_text_repel(label = rownames(dt_noOut), max.overlaps = 35, force = 10) +
  ggtitle("Clustering")
# Scree Plot
var_exp <- pca_res$sdev^2/sum(pca_res$sdev^2) # Finding the total variance explained by each PC
# Creating Scree Plot
qplot(c(1:(length(var_exp))), var_exp) + # Change the numbers in c() to the length of var_exp
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
data_raw<-as.matrix(dt[,-c(1:2)])
# Clustering Analysis
distance <- dist(data_raw) 
mds1<-cmdscale(distance, k=2) 
plot(mds1, type='n')
text(mds1, labels=rownames(data_raw), cex=0.6, adj=0.5)

# Dendrogram
dist <- dist(as.matrix(dt[,-c(1:2)]), method = "euclidian", diag = T) # Calculating the euclidian distance between samples
hc <- hclust(dist) # Performing hierarchical clustering
plot(hc)