# RNA-QC.R

library(Seurat)
library(dplyr)
library(ggplot2)
library(glue)
library(gridExtra)
library(grid)
library(tidyr)

source("htp/scripts/helper-functions/utils.R")
source("htp/scripts/helper-functions/rna-utils.R")

# Create Save Directories
create.dir("htp-outs/QC-Output/RNA-QC")
create.dir("htp-outs/QC-Output/RNA-QC/plots")
create.dir("htp-outs/QC-Output/RNA-QC/plots/VlnPlots")
create.dir("htp-outs/QC-Output/RNA-QC/plots/FeatureScatterPlots")
create.dir("htp-outs/QC-Output/RNA-QC/qc-stats")
create.dir("htp-outs/QC-Output/RNA-QC/RDS")

# Load Data
seu_obj_list <- readRDS("htp-outs/QC-Output/HTO-QC/RDS/seu-list-singlets.RDS")

# Add percent.mt
seu_obj_list <- lapply(seu_obj_list, FUN = function(seu_obj) {
  seu_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu_obj, pattern = "^mt-") # for mouse
  
  return(seu_obj)
})

# QC Plots
hto_run_names <- names(seu_obj_list)

## Vln Plot
for (name in hto_run_names) {
  pdf(glue("htp-outs/QC-Output/RNA-QC/plots/VlnPlots/{name}-RNA-QC-VlnPlot.pdf"), height = 12, width = 16)
  plot.rna.qc.vln.plots(seu_obj_list[[name]])
  dev.off()
}

## Scatter Plot
for (name in hto_run_names) {
  pdf(glue("htp-outs/QC-Output/RNA-QC/plots/FeatureScatterPlots/{name}-RNA-QC-FeatureScatterPlot.pdf"), height = 12, width = 16)
  plot.rna.qc.scatter.plot(seu_obj_list[[name]])
  dev.off()
}

# Get Summary Stats
for (name in hto_run_names) {
  get.rna.qc.summary(seu_obj_list[[name]]) %>%
    write.csv(glue("htp-outs/QC-Output/RNA-QC/qc-stats/{name}-rna-qc-summary.csv"))
}

# Merge Seurat Object
merged_seu_obj <- merge(seu_obj_list[[1]], y=c(seu_obj_list[-1]))

# Filter Seurat Object **** MOST IMPORTANT STEP ***
merged_seu_obj <- subset(merged_seu_obj, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 10)

# Save Seurat Objects
saveRDS(seu_obj_list, "htp-outs/QC-Output/RNA-QC/RDS/seu-list-post-qc.RDS")
saveRDS(merged_seu_obj, "htp-outs/QC-Output/RNA-QC/RDS/seu-merged-post-qc.RDS")




