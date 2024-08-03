# This script is to run MOFA to integrate and analyse the transcriptomics, proteomics, and metabolomics data from the PTE experiment.
# To run MOFA, a MOFA object must be created with four dimensions (samples, features, views (feature group), and groups (for samples))

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(MOFA2)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(GGally)

#----Creating Long Data Frame----
# Loading data files
# Transcriptomics
trans <- read.csv("PTE_Transcriptomics_MOFA.csv", stringsAsFactors = F, sep = ",")
colnames(trans)[1] <- "sample"
trans <- trans %>%
  add_column(group = c(rep("Naive", 9), rep("PTE",10), rep("Sham",11), rep("TBI", 13)), .after = "sample")
trans$sample <- gsub("Naïve", "Naive", trans$sample)
# Proteomics
prot <- read.csv("PTE_Proteomics_Pre-processed_MOFA.csv", stringsAsFactors = F, sep = ",", row.names = 2)
prot <- prot %>%
  select(!Gene.Name) %>%
  select(order(colnames(.))) %>% 
  t() %>%
  as.data.frame() %>%
  add_column(group = c(rep("Naive",10), rep("PTE",10), rep("Sham",12), rep("TBI",13)), 
             .before = 1) %>%
  add_column(sample = rownames(.), .before = 1) %>%
  remove_rownames()
prot$sample <- gsub("_"," ", prot$sample)
prot$sample <- gsub("Naïve", "Naive", prot$sample)
# Metabolomics
met <- read.csv("PTE_Metabolomics_Pre-processed_SampImp_PCA.csv", stringsAsFactors = F, sep = ",")
colnames(met)[1] <- "sample"
colnames(met)[2] <- "group"
met[,c(1:2)] <- lapply(met[,c(1:2)], gsub, pattern = "Naïve", replacement = "Naive")

# Removing outliers from the individual datasets
prot_noOut  <- prot[-c(1,18,22),] 
met_noOut <- met[-c(1,28),]
# Keeping PTE and TBI samples
prot_noOut <- prot_noOut[prot_noOut$group == "PTE" | prot_noOut$group == "TBI",] 
met_noOut <- met_noOut[met_noOut$group == "PTE" | met_noOut$group == "TBI",] #[c(11:20,33:45),]
trans_noOut <- trans[trans$group == "PTE" | trans$group == "TBI",]  #[c(10:19,31:43),]

# Getting the common samples left
prot_common <- prot_noOut[prot_noOut$sample %in% met_noOut$sample,]
met_common <- met_noOut[met_noOut$sample %in% prot_noOut$sample,]
trans_common <- trans_noOut[-8,]  #trans[-c(1, 17, 20, 26),]
# Remove Naive 7, Sham 1 from prot and met
prot_common <- prot_common[-c(7,19),]
met_common <- met_common[-c(7,19),]

# Unit variance scaling the data (make samples rows and features columns).
prot_scale <- 
  data.frame(scale(prot_common[,-c(1:2)]), center = T, scale = T) %>%
  add_column(., prot_common[,1:2], .before = 1) %>%
  select(!c(center, scale))

met_scale <- 
  data.frame(scale(met_common[,-c(1:2)]), center = T, scale = T) %>%
  add_column(., met_common[,1:2], .before = 1) %>%
  select(!c(center, scale))

trans_scale <- 
  data.frame(scale(trans_common[,-c(1:2)], center = T, scale = T)) %>%
  add_column(., trans_common[,1:2], .before = 1) 

# Turning into long data frame format
# Proteomics
prot_long <- prot_scale %>%
  pivot_longer(!c("sample", "group"), # all columns except strain
               names_to = "feature", # proteins to go into column; 'features'
               values_to =  "value") %>% # values to go into column; 'value'
  arrange(feature) %>% # arranging by proteins (ascending order) 
  mutate(view = "Proteomic", .after = feature) # add 'view' column (or in this case, omic)
# Metabolomics
met_long <- met_scale %>%
  pivot_longer(!c("sample", "group"),
               names_to = "feature",
               values_to = "value") %>%
  arrange(feature) %>%
  mutate(view = "Metabolomic", .after = feature)
# transcriptomics
trans_long <- trans_scale %>%
  pivot_longer(!c("sample", "group"),
               names_to = "feature",
               values_to = "value") %>%
  arrange(feature) %>%
  mutate(view = "Transcriptomic", .after = feature)
# Concatenating both data frames
join <- rbind(trans_long,prot_long, met_long)

# Saving for future use
write.csv(join, "PTE_MOFA_TBIPTEOnly_20240509.csv", row.names = F)

#---- Training the Model----
# Load in dataframe
df <- read.csv("PTE_MOFA_NoOut_New_20240517.csv", stringsAsFactors = F, sep = ",", row.names = 1)

# Keeping top 10,000 genes
top_ord <- read.csv("PTE_Transcriptomics_CPMord.csv", stringsAsFactors = F, row.names = 1)
top_ord <- as.data.frame(t(top_ord))
#Removing 
remove_ids <- top_ord %>% slice_tail(n = 5186)
df_filt <- df[!df$feature %in% rownames(remove_ids),]

# Creating a MOFA object (no multi-group intferance)
df_nogrp <- df_filt[,-2]
# Making a new MOFAobject
MOFAobj <- create_mofa(df_nogrp)
print(MOFAobj)
plot_data_overview(MOFAobj)

# Data options
data_opts <- get_default_data_options(MOFAobj)
# If you want to change something:
data_opts$scale_views <- F
head(data_opts)

# Model options
model_opts <- get_default_model_options(MOFAobj)
# Changing options
model_opts$num_factors <- 5
model_opts$spikeslab_weights <- TRUE # introduces sparsity
head(model_opts)

# Training options
train_opts <- get_default_training_options(MOFAobj)
head(train_opts)
# Modifying training options
train_opts$maxiter <- 1500
train_opts$convergence_mode <- "medium"

# Now we build the final MOFA object
model_train <- prepare_mofa(object = MOFAobj,
                            data_options = data_opts,
                            model_options = model_opts,
                            training_options = train_opts)
# Training the object. We may need to tweak the above options according to the final training result
# Ideally we will select the model with the highest Evidence Lower Bound (ELBO)
outfile = file.path(getwd(),"PTE_MOFA_ProtTransOnly_TBIPTEOnly.hdf5")
MOFAobj.trained <- run_mofa(model_train, outfile, save_data = T, use_basilisk = T)

#---- Downstream Analysis----
# Loading in the file
filepath <- file.path(getwd(), "PTE_MOFA_TBIPTEOnly_20240509.hdf5")
trained_model <- load_model(filepath)
plot_data_overview(trained_model) # looking at the data overview

# Adding metadata
exp_des <- read.csv("PTE_Metadata.csv", header = T, stringsAsFactors = F)
new_lab <- c("sample", "PTE","OFT\nTime", "EPM\nTime", "Sucrose %", 
             "MWM\nAcqusition", "MWM\nReversal", "Beam Slips\nand Falls" )
exp_des <- exp_des %>%
  select(!Debbie.s.ID) %>%
  arrange(Sample.ID) %>%
  filter(Sample.ID %in% samples_metadata(trained_model)$sample)
colnames(exp_des) <- new_lab
# Change condition label
exp_des$condition <- c(rep("Naive",8), rep("PTE",9), rep("Sham",9), rep("TBI",13)) 
exp_des$PTE <- c(rep(1, 9), rep(0, 13))  
# adding it to the model
samples_metadata(trained_model) <- exp_des
samples_metadata(trained_model) # Checking

# Checking correlation between factors
plot_factor_cor(trained_model,
                col=colorRampPalette(c("blue", "white", "red"))(200))

# Now we can quantify the amount of variance explained by individual data modules
# We can look at either the total variance per view or per group
head(trained_model@cache$variance_explained$r2_total[[1]])
head(trained_model@cache$variance_explained$r2_per_factor[[1]])
# Plotting the amount of variance explained
plot_variance_explained(trained_model, x="group", y="factor", plot_total = T)[[2]] # Total variance
plot_variance_explained(trained_model, x="view", y="factor") # Per factor

# Plotting Factors
# Visualisation of single factors
plot_factor(trained_model,
            factors = 1,
            color_by = "condition",
            legend = T)
# Visualisation of multiple factors
plot_factors(trained_model, factors = 1:5, color_by = "condition", legend = T)

# Adding labels to the factor plots
plts <- vector('list',length(colnames(trained_model@expectations[["Z"]][["single_group"]])))
for (factor in 1:length(colnames(trained_model@expectations[["Z"]][["single_group"]]))) {
  plts[[factor]] <- plot_factor(trained_model,
                                factors = factor,
                                color_by = "condition",
                                legend = T)
  names(plts) <- colnames(trained_model@expectations[["Z"]][["single_group"]])
}
plts_df <- vector('list', length(plts))
for (factor in seq(length(plts))) {
  plts_df[[factor]] <- plts[[factor]][["plot_env"]][["df"]]
  names(plts_df) <- colnames(trained_model@expectations[["Z"]][["single_group"]])
}
names(plts_df) <- lapply(names(plts_df), gsub, pattern = "([a-z])([0-9])", replacement = "\\1 \\2")
# Plotting using ggplot2
for (factor in seq(length(plts_df))) {
  p <- ggplot(data = plts_df[[factor]], mapping = aes(y = value, x = group_by)) + 
    geom_jitter(position = position_jitter(seed = 1), aes(colour = color_by)) + 
    labs(color = "Condition") +
    theme_classic() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          legend.title = element_text()) +
    geom_text_repel(aes(label = sample), position = position_jitter(seed = 1), max.overlaps = Inf) + 
    geom_hline(yintercept = 0, linetype = "dashed", linewidth=0.2, alpha=0.5) +
    ggtitle(names(plts_df)[[factor]])
  print(p)
}

# Saving Plots
pdf("PTE_TBIPTEOnly_TransProtOnly_FactorPlts.pdf")
plot_data_overview(trained_model)
plot_factor_cor(trained_model,
                col=colorRampPalette(c("blue", "white", "red"))(200))
# Total variance
plot_variance_explained(trained_model, x="group", y="factor", plot_total = T)[[2]]
# Per factor
plot_variance_explained(trained_model, x="view", y="factor")
plot_factors(trained_model, factors = 1:5, color_by = "condition", legend = T)
for (factor in seq(length(plts_df))) {
  p <- ggplot(data = plts_df[[factor]], mapping = aes(y = value, x = group_by)) + 
    geom_jitter(position = position_jitter(seed = 1), aes(colour = color_by)) + 
    labs(color = "Condition") +
    theme_classic() + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          legend.title = element_text()) +
    geom_text_repel(aes(label = sample), position = position_jitter(seed = 1), max.overlaps = Inf) + 
    geom_hline(yintercept = 0, linetype = "dashed", linewidth=0.2, alpha=0.5) +
    ggtitle(names(plts_df)[[factor]])
  print(p)
}
# Association analysis with metadata
# Correlation plot
correlate_factors_with_covariates(object = trained_model,
                                  covariates = c("condition","OFT\nTime", "EPM\nTime", "Sucrose %",
                                                 "MWM\nAcqusition", "MWM\nReversal", 
                                                 "Beam Slips\nand Falls"),
                                  plot = "r",
                                  col=colorRampPalette(c("blue", "white", "red"))(200),
                                  mar = c(0,0,1,0),
                                  win.asp = 0.5)
# Heatmap
correlate_factors_with_covariates(object = trained_model,
                                  covariates = c("condition", "OFT\nTime", "EPM\nTime", "Sucrose %",
                                                 "MWM\nAcqusition", "MWM\nReversal", 
                                                 "Beam Slips\nand Falls"),
                                  plot = "log_pval",
                                  alpha = 0.05,
                                  win.asp = 0.5)
dev.off()


#----Weight Plots----
# Changing the labels of the molecules
feat_name <- features_names(trained_model) # Getting the name of features for the model

# Proteomics
p_lab <- read.csv("Gene_ProtID.csv", stringsAsFactors = F)
p_trained <- feat_name[["Proteomic"]] # Getting proteomics data from the trained_model
p_labMatch <- p_lab[p_lab$Protein.IDs %in% p_trained,] # Subset the main label dataframe 
p_labMatch <- p_labMatch[order(match(p_labMatch$Protein.IDs, p_trained)),] # Ordering against trained_model
# Metabolomics
m_lab <- read.csv("Metabolie_NewLabels.csv", stringsAsFactors = F, header = 1)
m_trained <- feat_name[["Metabolomic"]]
m_labMatch <- m_lab[m_lab$old_label %in% m_trained,] # Checking what is missing
m_labMatch <- m_labMatch[order(match(m_labMatch$old_label, m_trained)),] 

# Making a new list and changing
feat_list <- list(Metabolomic = m_labMatch$new_label,
                  Proteomic = p_labMatch$Gene.Name,
                  Transcriptomic = feat_name$Transcriptomic)
# Changing feature name
features_names(trained_model) <- feat_list

# Plotting weights
pdf("PTE_MOFA_NoOuts_Weights.pdf")
for (i in 1:5) {
  p <- plot_top_weights(trained_model,
                        view = "Proteomic",
                        factors = i,
                        nfeatures = 20,
                        scale = T)
  m <- plot_top_weights(trained_model,
                        view = "Metabolomic",
                        factors = i,
                        nfeatures = 20,
                        scale = T)
  t <- plot_top_weights(trained_model,
                        view = "Transcriptomic",
                        factors = i,
                        nfeatures = 20,
                        scale = T)
  plot(t)
  plot(p)
  plot(m)
}
dev.off()