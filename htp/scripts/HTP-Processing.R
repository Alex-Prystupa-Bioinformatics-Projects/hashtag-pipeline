# HTP-Processing.R

library(Seurat)
library(dplyr)
library(ggplot2)
library(glue)
library(gridExtra)
library(grid)
library(tidyr)
library(SCpubr)
library(harmony)

source("htp/scripts/helper-functions/utils.R")
source("htp/scripts/helper-functions/processing-utils.R")

# Create Save Directories
create.dir("htp-outs/Downstream-Analysis")
create.dir("htp-outs/Downstream-Analysis/plots")
create.dir("htp-outs/Downstream-Analysis/plots/UMAP")
create.dir("htp-outs/Downstream-Analysis/plots/UMAP/specific-resolutions")
create.dir("htp-outs/Downstream-Analysis/markers")
create.dir("htp-outs/Downstream-Analysis/RDS")

# Load Data
seu_obj <- readRDS("htp-outs/QC-Output/RNA-QC/RDS/seu-merged-post-qc.RDS")

# Run Pipeline
seu_obj <- multi_sample_pipeline(seu_obj)
saveRDS(seu_obj, "htp-outs/Downstream-Analysis/RDS/final-htp-seu-obj.RDS")

# Find Markers
markers_results <- run.presto.specific.columns(seu_obj, RNA_cols = T)
markers_results$`Raw-Markers` %>% writexl::write_xlsx("htp-outs/Downstream-Analysis/markers/Raw-Markers.xlsx")
markers_results$`Filtered-Markers` %>% writexl::write_xlsx("htp-outs/Downstream-Analysis/markers/Filtered-Markers.xlsx")
markers_results$`Top-Markers` %>% writexl::write_xlsx("htp-outs/Downstream-Analysis/markers/Top-Markers.xlsx")

# Plots
plot.umaps(seu_obj, save_path = "htp-outs/Downstream-Analysis/plots/UMAP")




