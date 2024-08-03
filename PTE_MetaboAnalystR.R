# This is to run MetaboAnalystR for the PTE metabolomics analysis
# We will run metabolite set enrichment analysis (MSEA).
# This is based on : https://rdrr.io/github/simscr/metaboanalyst/f/vignettes/Introduction_to_MetaboAnalystR.Rmd 

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(tidyverse)
library(MetaboAnalystR)

# Preparing for the analysis
m_KEGG <- read.csv("PTE Enrichment Analysis/PTE_TBIPTEOnly_MA_KEGG.csv", stringsAsFactors = F)
m_noNA <- m_KEGG %>%
  select(-contains("NA")) # Removing metabolites without KEGG ID
write.csv(m_noNA, "PTE Enrichment Analysis/PTE_TBIPTEOnly_MA_noNA.csv", row.names = F)

# Running the analysis
# This creates a dataSet object for storing processed data (in the form of an R list), a analSet object, imgSet object, and msg.vec
mSet <- InitDataObjects("conc", "msetqea", FALSE) # Define the dataType and analysis type

# Read in the data and fill in the dataSet list
# Specify the file path, data format and label type
mSet <- Read.TextData(mSet, "PTE_TBIPTEOnly_MA_noNA.csv", 
                      "rowu", # Here, samples are in rows and unpaired (so "rowu") 
                      "disc") # Discrete data
# To view messages from the data import and processing
mSet$msgSet$read.msg

# Data integrity check
# Ensure data is valid and suitable for subsequent analysis
mSet<-SanityCheckData(mSet)

# Pre-processing data
# Imputation
# Default is 1/5 of the min positive values of the valriables
mSet <- ImputeMissingVar(mSet, method="min") # Default imputation method
mSet$msgSet$replace.msg
# Matching compound name to KEGG IDs
mSet<-CrossReferencing(mSet, "kegg");

# Preparing for downstream analysis
mSet<-CreateMappingResultTable(mSet)
# Normalisation
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 300, width=NA) # feature-wise view of the data normalization
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 300, width=NA) # sample-wise view of the normalization
# Filtering
mSet<-SetMetabolomeFilter(mSet, F);

# Enrichment analysis
# Selecting library
mSet<-SetCurrentMsetLib(mSet, "RaMP_pathway", 2); # Give library name and min no. of compounds within the set
mSet<-CalculateGlobalTestScore(mSet) # perform QEA with globaltest
# Plotting
mSet<-PlotQEA.Overview(mSet, "qea_0_", "net", "png", 300, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "qea", "qea_dot_0_", "png", 300, width=NA)

# Create Biomarker Sweave report 
PreparePDFReport(mSet, "PTE_MAR_Enrichment")

# To save all files created during your session
SaveTransformedData(mSet)
# Saving as RData session
save.image(file = "MA_EnrichmentAnalysis_PTETBIOnly.RData")
load("MA_EnrichmentAnalysis_PTETBIOnly.RData")