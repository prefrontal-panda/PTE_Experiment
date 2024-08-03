# This is to run transcriptomics differential expression analysis using edgeR

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(edgeR)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

# Reading in the files
count <- read.csv("PTE Transcriptomics/PTE_salmon_counts_length_scaled.csv")
colnames(count) <- gsub("^X","", colnames(count))
# Reading in new name file
newLab <- read.csv("PTE_Transcript_Labels.csv")

# Making DEGList
count_list <- count_list[,newLab$ID] # Rearranging
colnames(count_list) <- newLab$New_ID
rownames(count_list) <- count$gene_id
# Grouping information
lab <- c(rep("Naive",10), rep("PTE",10), rep("Sham",12), rep("TBI",13))
# Sample information
samps <- colnames(count_list)
# Setting DGEList
list <- DGEList(counts = count_list, group = lab, genes = count[,c(1:2)],
                samples = samps)
list$samples # to see library size of the samples

# Filtering
keep <- filterByExpr(list) # Getting list of genes to remove
table(keep)
list_filt <- list[keep, , keep.lib.sizes = F] # performing filtering

# Plotting density plot
# Setting up
nsamples <- ncol(list)
col.density <- colorRampPalette(brewer.pal(12, name = "Paired"))(nsamples)
# Un-filtered data
lcpm_dens <- cpm(list, log = T)
plot(density(lcpm_dens[, 1], bw = 0.4), col = col.density[1], lwd = 2, ylim = c(0, 1),
     las = 2, main = "", xlab = "")
title(main = "Raw data", xlab = "Log-cpm")
for (i in 2:nsamples) {
  den <- density(lcpm_dens[, i], bw = 0.4)
  lines(den$x, den$y, col = col.density[i], lwd = 2)
}
legend("topright", samps, text.col = col.density, bty = "n", ncol = 5)
Raw_dens <- recordPlot()
# Filtered Data
lcpm_filt <- cpm(list_filt, log = T)
plot(density(lcpm_filt[, 1], bw = 0.4), col = col.density[1], lwd = 2, ylim = c(0, 1),
     las = 2, main = "", xlab = "")
title(main = "Filtered data", xlab = "Log-cpm")
for (i in 2:nsamples) {
  den <- density(lcpm_filt[, i], bw = 0.4)
  lines(den$x, den$y, col = col.density[i], lwd = 2)
}
legend("topright", samps, text.col = col.density, bty = "n", ncol = 5)
Filt_dens <- recordPlot()

# Normalisation
list_filt <- normLibSizes(list_filt, method = "TMM")
list_filt$samples

# Boxplot
# Before normalisation
lcpm <- cpm(list, log = T)
par(mar=c(8,4,4,1))
boxplot(lcpm, las=2)
title(main = "Unormalised data", ylab = "Log-cpm")
Unnorm_boxplt <- recordPlot()
# After normalisation
lcpm2 <- cpm(list_filt, log = T)
par(mar=c(8,4,4,1))
boxplot(lcpm2, las=2)
title(main = "Normalised data", ylab = "Log-cpm")
norm_boxplt <- recordPlot()

# Removing outliers
par(mar=c(5,4,4,2))
plotMDS(list_filt)
  # Can also export data for PCA analysis
  # pca <- cpm(list_filt, log = T)
  # write.csv(pca, "PTE_Transcriptomics_PCA.csv", row.names = T)
outliers <- c("Sham 1", "Naïve 7")
lt_noOut <- list_filt[ ,which(!list_filt$samples$samples %in% outliers)]
# Check
lt_noOut$samples
head(lt_noOut[["counts"]])
# New MDS
plotMDS(lt_noOut)

#----Exporting for MOFA-----
lt_noOut_cpm <- cpm(lt_noOut, log = T) # Getting normalised dataframe for MOFA

# Checking data distribution
nsamp_new <- length(colnames(lt_noOut_cpm))
plot(density(lt_noOut_cpm[, 1], bw = 0.4), col = col.density[1], lwd = 2, ylim = c(0, 1),
     las = 2, main = "", xlab = "")
title(main = "Filtered + Normalised data", xlab = "Log-cpm")
for (i in 2:nsamp_new) {
  den <- density(lt_noOut_cpm[, i], bw = 0.4)
  lines(den$x, den$y, col = col.density[i], lwd = 2)
}
legend("topright", colnames(lt_noOut_cpm), text.col = col.density, bty = "n", ncol = 5)
fin_dens <- recordPlot()
# Changing rownames so that it is gene names
rownames(lt_noOut_cpm) <- lt_noOut$genes$gene_name

# Reorder by variance  to filter out lowly variable transcripts (do this before turning lt_noOut_cpm to dataframe)
library(matrixStats)
o <- order(rowVars(lt_noOut_cpm), decreasing = TRUE)
cpm_ord <- lt_noOut_cpm[o,]
cpm_ord <- as.data.frame(t(cpm_ord))
write.csv(cpm_ord, "PTE_Transcriptomics_CPMord.csv")
# Re-structuring
lt_noOut_cpm <- lt_noOut_cpm %>%
  as.data.frame(.) %>%
  select(order(colnames(.)))
lt_noOut_cpm <- as.data.frame(t(lt_noOut_cpm))
# Writing
write.csv(lt_noOut_cpm, "PTE_Transcriptomics_MOFA.csv",  row.names = T)

#----Saving Plots----
pdf("PTE_Transcriptomics_Pre-processing Plots.pdf")
Raw_dens
Filt_dens
Unnorm_boxplt
norm_boxplt
par(mar=c(5,4,4,2))
plotMDS(list_filt)
plotMDS(lt_noOut)
fin_dens
dev.off()

# Saving as RData
save.image(file = "PTE_Transcriptomics_Preprocessing.RData")
# Loading
load("PTE_Transcriptomics_Preprocessing.RData") 

#----Differential expression testing----
# Creating labels and batch informations
lab_noOut <- as.factor(c(rep("Naïve",9), rep("PTE",10), rep("Sham",11), rep("TBI",13)))
batch <- as.factor(c(1,3,rep(1,3), rep(2,2),3,3, 
                     1,3,rep(1,2), rep(2,3), rep(3,3),
                     rep(3,3), rep(1,3), rep(2,4),3,
                     1, rep(3,4), rep(1,3), rep(2,5)))
# Creating design and contrast matrix
des <- model.matrix(~0 + lab_noOut + batch)
rownames(des) <- colnames(lt_noOut)
colnames(des) <- gsub("lab_noOut","", colnames(des))

cont_anova <- makeContrasts(
  # Contrasts of interest
  `Sham vs Naïve` = "Sham - Naïve", # effects of craniotomy 
  `TBI vs Naïve` = "TBI - Naïve", # changes in TBI
  `PTE vs Naïve` = "PTE - Naïve", # changes in PTE
  `PTE vs TBI` = "PTE - TBI", # changes in PTE post-TBI
  #`PTE vs Sham` = "PTE - Sham",
  
  # Levels
  levels = des
)

# Dispersion Estimate
lt_noOut <- estimateDisp(lt_noOut, des, robust = T)
lt_noOut$common.dispersion
plotBCV(lt_noOut)

# Running the individual contrasts
# Setting variables
top_n <- 10
colors <- c(Up = "firebrick3", Down = "dodgerblue3", `Not Sig` = "gray50")
# Fitting model
fit <- glmQLFit(lt_noOut, des, robust = T)
# For loop
pdf("PTE_Trans_DEAnalysisPlots_Batch.pdf")
for (comp in colnames(cont_anova)) {
  qlf <- glmQLFTest(fit,contrast=cont_anova[,comp])
  
  # 2 FC
  FC2 <- decideTests(qlf, adjust.method = "BH", p.value = 0.05, lfc = 1)
  sum_FC2 <- as.data.frame(summary(FC2))
  colnames(sum_FC2) <- c("Direction", "Pair", "Number")
  write.csv(sum_FC2, file = paste0("PTE_Trans_", comp, "_DESummary_FC2.csv"), row.names = F)
  
  # 1.5 FC
  FC2 <- decideTests(qlf, adjust.method = "BH", p.value = 0.05, lfc = 0.585)
  sum_FC2 <- as.data.frame(summary(FC2))
  colnames(sum_FC2) <- c("Direction", "Pair", "Number")
  write.csv(sum_FC2, file = paste0("PTE_Trans_", comp, "_DESummary_FC2.csv"), row.names = F) 
  
  # Top transcripts
  tab <- topTags(qlf, adjust.method = "BH", sort.by = "PValue", n = Inf)
  
  # Plotting MD
  is.de <- decideTests(qlf)
  plotMD(qlf, status = is.de, values=c(1,-1), col=c("red","blue"),
         legend = "topright", main = comp)
  
  # Volcano Plots
  # FC 1.5
  # Adding direction of differentiation
  tab_mod <- tab$table %>%
    mutate("Sig DE" = case_when(logFC >= log2(1.5) & FDR <= 0.05 ~ "Up",
                                logFC <= -log2(1.5) & FDR <= 0.05 ~ "Down",
                                TRUE ~ "Not Sig"))
  # Only taking top 10
  top_genes_1.5 <- rbind(
    tab_mod %>%
      filter(`Sig DE` == "Up") %>%
      arrange(FDR, desc(abs(logFC))) %>%
      head(top_n),
    tab_mod %>%
      filter(`Sig DE` == "Down") %>%
      arrange(FDR, desc(abs(logFC))) %>%
      head(top_n)
  )
  p_1.5 <- ggplot(tab_mod, aes(logFC, -log(PValue,10))) + 
    geom_point(aes(color = `Sig DE`), size = 1) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"(P-Value)")) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    geom_text_repel(data = top_genes_1.5,
                    aes(logFC, -log(PValue,10), 
                        label = gene_name),
                    size = 3,
                    max.overlaps = 50) +
    ggtitle(paste0(comp, ", FC = 1.5"))
  print(p_1.5)
  # Saving
  write.csv(tab_mod, file = paste0("PTE_Trans_", comp, "_DEAll_FC1.5.csv"), row.names = T)
  # FC 2
  tab_mod_FC2 <- tab$table %>%
    mutate("Sig DE" = case_when(logFC >= log2(2) & FDR <= 0.05 ~ "Up",
                                logFC <= -log2(2) & FDR <= 0.05 ~ "Down",
                                TRUE ~ "Not Sig"))
  # Plotting
  top_genes_2 <- rbind(
    tab_mod_FC2 %>%
      filter(`Sig DE` == "Up") %>%
      arrange(FDR, desc(abs(logFC))) %>%
      head(top_n),
    tab_mod_FC2 %>%
      filter(`Sig DE` == "Down") %>%
      arrange(FDR, desc(abs(logFC))) %>%
      head(top_n)
  )
  p_2 <- ggplot(tab_mod_FC2, aes(logFC, -log(PValue,10))) +
    geom_point(aes(color = `Sig DE`), size = 1) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"(P-Value)")) +
    scale_color_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) +
    geom_text_repel(data = top_genes_2,
                    aes(logFC, -log(PValue,10), 
                        label = gene_name),
                    size = 3,
                    max.overlaps = 50) +
    ggtitle(paste0(comp, ", FC = 2"))
  print(p_2)
  # Saving
  write.csv(tab_mod_FC2, file = paste0("PTE_Trans_", comp, "_DEAll_FC2.csv"), row.names = T)
}
dev.off()

# Saving as RData
save.image(file = "PTE_Transcriptomics_DEAnalysis.RData")
load("PTE Transcriptomics/PTE_Transcriptomics_DEAnalysis.RData")