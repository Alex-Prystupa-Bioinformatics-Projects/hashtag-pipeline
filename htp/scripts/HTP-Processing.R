# HTP-Processing.R

library(argparse)
library(Seurat)
library(dplyr)
library(ggplot2)
library(glue)
library(gridExtra)
library(grid)
library(tidyr)
library(SCpubr)
library(harmony)
options(future.globals.maxSize = 4 * 1024^3)  # Harmony failing without this

source("htp/scripts/helper-functions/utils.R")
source("htp/scripts/helper-functions/processing-utils.R")

# Parse command line arguments
parser <- ArgumentParser(description = "HTP Processing pipeline")
parser$add_argument("--res_cols", default = "RNA", choices = c("RNA", "SCT"),
                    help = "Which assay resolution columns to use (RNA or SCT). Default is RNA.")
args <- parser$parse_args()

res_cols <- args$res_cols

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
# Make sure multi_sample_pipeline is updated to accept res_cols param and use it (e.g. assay.use = res_cols)
seu_obj <- multi_sample_pipeline(seu_obj, res_cols = res_cols)

saveRDS(seu_obj, "htp-outs/Downstream-Analysis/RDS/final-htp-seu-obj.RDS")

# Find Markers
markers_results <- run.presto.specific.columns(seu_obj, res_cols = res_cols)

markers_results$`Raw-Markers` %>% writexl::write_xlsx("htp-outs/Downstream-Analysis/markers/Raw-Markers.xlsx")
markers_results$`Filtered-Markers` %>% writexl::write_xlsx("htp-outs/Downstream-Analysis/markers/Filtered-Markers.xlsx")
markers_results$`Top-Markers` %>% writexl::write_xlsx("htp-outs/Downstream-Analysis/markers/Top-Markers.xlsx")

# Plots
plot.umaps(seu_obj, save_path = "htp-outs/Downstream-Analysis/plots/UMAP", res_cols = res_cols)
plot.umaps.labeled(seu_obj, save_path = "htp-outs/Downstream-Analysis/plots/UMAP", res_cols = res_cols)
