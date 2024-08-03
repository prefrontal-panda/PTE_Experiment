# This is to run ReactomGSA to perform PADOG enrichment analysis on the transcriptomic and proteomic datasets for the PTE experiment.
# We use this as there are several additional functionalities not found online
# Based on: https://bioconductor.org/packages/release/bioc/vignettes/ReactomeGSA/inst/doc/using-reactomegsa.html 

# Setting working directory
setwd("//ad.monash.edu/home/User028/dcho0009/Desktop/r")
getwd()

# Loading libraries
library(ReactomeGSA)
library(tidyverse)

# Reading in the files required
# If using a data.frame for the counts/intensities, ensure samples are in columns
trans <- read.csv("PTE_Trans_Raw_Reactome_TBIPTEOnly.csv", 
                  stringsAsFactors = F, row.names = 1)
prot <- read.csv("PTE_Proteomics_UniProtID_Reactome.csv",
                 stringsAsFactors = F, row.names = 1)
# Metadata (row = sample, column = metadata property)
trans_meta <- data.frame(samples = colnames(trans),
                         group = c("TBI", "TBI", "PTE", "PTE", "TBI", "TBI",
                                   "PTE", "PTE", "TBI", "PTE", "TBI", "TBI",
                                   "TBI", "PTE", "TBI", "TBI", "TBI", "PTE",
                                   "TBI", "PTE", "PTE", "TBI", "PTE"))
prot_meta <- data.frame(samples = colnames(prot),
                        group = c("PTE", "TBI", "PTE", "TBI", "PTE", "TBI",
                                  "TBI",  "TBI", "PTE", "TBI", "PTE", "TBI",
                                  "PTE", "TBI", "TBI", "PTE", "PTE", "TBI",
                                  "PTE", "TBI", "PTE", "TBI", "TBI"))

# Checking the data types that can be analysed
get_reactome_data_types()

# Checking which methods are available
available_methods <- get_reactome_methods(print_methods = FALSE, return_result = TRUE)
available_methods$name # only show the names of the available methods
# Checking the parameters
params <- available_methods$parameters[available_methods$name == "PADOG"][[1]] # show the parameter names for the method
paste0(params$name, " (", params$type, ", ", params$default, ")") # print

# To perform a GSEA
# First, start an analysis request
gsea_req <- ReactomeAnalysisRequest(method = "PADOG")
# Setting the parameters (as required)
gsea_req <- set_parameters(request = gsea_req, max_missing_value = 0.5)

# Adding dataset
# Transcriptome
gsea_req <- add_dataset(request = gsea_req,
                        expression_values = trans, # dataset
                        name = "Transcriptomics",
                        type = "rnaseq_counts",
                        comparison_factor = "group", # Value here should match column name in metadata
                        comparison_group_1 = "PTE",
                        comparison_group_2 = "TBI",
                        sample_data = trans_meta, # metadata
                        additional_factors = NULL,
                        overwrite = F,
                        discrete_norm_function = "TMM")
gsea_req # Checking
# Proteomics
gsea_req <- add_dataset(request = gsea_req,
                        expression_values = prot,
                        name = "Proteomics",
                        type = "proteomics_int",
                        comparison_factor = "group", # Value here should match column name in metadata
                        comparison_group_1 = "PTE",
                        comparison_group_2 = "TBI",
                        additional_factors = NULL,
                        sample_data = prot_meta,
                        overwrite = F)
gsea_req # Checking

# Running the analysis
res <- perform_reactome_analysis(request = gsea_req, compress = F)

# Saving
save.image(file = "PTE_TBIPTEOnly_ReactomeGSA_PADOG.RData")
load("PTE_TBIPTEOnly_ReactomeGSA_PADOG.RData")

#----Result Investigation----
names(res) # Names of dataset
result_types(res) #result types available
# Getting the fold change result
prot_fc <- get_result(res, type = "fold_changes", name = "Proteomics")
trans_fc <- get_result(res, type = "fold_changes", name = "Transcriptomics")
# Getting pathway results
path <- pathways(res)
common_path <- path %>%
  filter(sig.Transcriptomics == T & sig.Proteomics == T)
up_path <- path %>%
  filter(Direction.Transcriptomics == "Up" & Direction.Proteomics == "Up")
down_path <- path %>%
  filter(Direction.Transcriptomics == "Down" & Direction.Proteomics == "Down")

# Visualising Results
plot_volcano(res, 2) # Volcano plot
plot_correlations(res) # Pathway-level similarity
plot_heatmap(res, fdr = 0.05, break_long_names = T) # Plotting direction of change within the two datasets

# Getting heatmap of overall top 20 pathways
comm_path <- plot_heatmap(res, fdr = 0.05, max_pathways =  20, break_long_names = T) +
  ggplot2::theme(axis.text.x = element_text(size = 9),
                 axis.title.y = element_blank(),
                 axis.text.y = element_text(size = 8),
                 title = element_text(size = 10)) + 
  ggtitle("FDR = 0.05")
# Saving
ggsave(filename = "PTETBIOnly_ReactomeGSA_CommonPath.tiff", plot = comm_path, 
       width = 20, height = 20, units = "cm", device='tiff', dpi=300)

# Getting individual omics datasets
trans_path <- path %>%
  select(c("Name",contains("Transcriptomics")))
prot_path <- path %>%
  select(c("Name",contains("Proteomics")))
# Saving dataframes
write.csv(prot_path, "PTE_PTETBIOnly_ReactomeGSA_Protein.csv", row.names = T)