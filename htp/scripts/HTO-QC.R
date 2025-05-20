library(Seurat)
library(dplyr)
library(ggplot2)
library(glue)
library(gridExtra)
library(grid)

source("htp/scripts/helper-functions/utils.R")
source("htp/scripts/helper-functions/hto-utils.R")

# Create Save Directories
create.dir("htp-outs/QC-Output/HTO-QC")
create.dir("htp-outs/QC-Output/HTO-QC/plots")
create.dir("htp-outs/QC-Output/HTO-QC/qc-stats")
create.dir("htp-outs/QC-Output/HTO-QC/RDS")

# Load Samples Sheet
samples_sheet <- readxl::read_xlsx("htp/raw-data/hto-samples-sheet.xlsx")

# Load Samples Groups
samples_groups <- load.xlsx.multi.sheet("htp/raw-data/hto-samples-groups.xlsx")

# Run HTO Preprocessing Workflow
seu_obj_list <- list()

for (i in 1:nrow(samples_sheet)) {
  path <- samples_sheet$path[i]
  project_name <- samples_sheet$hto_run_name[i]
  samples_groups_df <- samples_groups[[project_name]]
  
  seu_obj_list[[i]] <- run.hto.workflow(path = path, project_name = project_name, samples_groups_df = samples_groups_df)
}

names(seu_obj_list) <- samples_sheet$hto_run_name

# Get HTO Classification Stats Pre-Singlet Filtering
hto_demux_stats <- get.hto.demux.stats(seu_obj_list)
hto_demux_stats %>% write.csv("htp-outs/QC-Output/HTO-QC/qc-stats/demux-stats.csv")

# Keep Singlets
seu_obj_list_singlets <- lapply(seu_obj_list, function(seu_obj) {return(seu_obj %>% subset(subset = HTO_classification.global == "Singlet"))})

# Get Filtering Stats
for (name in names(seu_obj_list)) {
  filter.stats(seu_obj = seu_obj_list[[name]], seu_obj_singlet = seu_obj_list_singlets[[name]], hto_run = name) %>%
    write.csv(glue("htp-outs/QC-Output/HTO-QC/qc-stats/{name}-filter-stats.csv"))
}

# Feature Scatter Plots
for (name in names(seu_obj_list)) {
  pdf(glue("htp-outs/QC-Output/HTO-QC/plots/{name}-Feature-Scatter.pdf", height = 10, width = 12))
  plot.feature_scatter(seu_obj = seu_obj_list[[name]]) %>% print()
  dev.off()
  
  pdf(glue("htp-outs/QC-Output/HTO-QC/plots/{name}-Singlet-Feature-Scatter.pdf", height = 10, width = 12))
  plot.feature_scatter(seu_obj = seu_obj_list_singlets[[name]]) %>% print()
  dev.off()
}

# Save RDS Files
saveRDS(seu_obj_list, "htp-outs/QC-Output/HTO-QC/RDS/seu-list.RDS")
saveRDS(seu_obj_list_singlets, "htp-outs/QC-Output/HTO-QC/RDS/seu-list-singlets.RDS")



