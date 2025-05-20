# processing-utils.R

# 1. Multi Sample Pipeline
multi_sample_pipeline <- function(merged_seu_obj, split_col="orig.ident", harmony=FALSE) {
  
  merged_seu_obj <- NormalizeData(merged_seu_obj)
  merged_seu_obj <- FindVariableFeatures(merged_seu_obj)
  merged_seu_obj <- ScaleData(merged_seu_obj)
  merged_seu_obj <- RunPCA(merged_seu_obj)
  merged_seu_obj <- RunHarmony(merged_seu_obj, assay.use="RNA", group.by.vars = "orig.ident")
  merged_seu_obj <- RunUMAP(merged_seu_obj, reduction = "harmony", dims = 1:40)
  merged_seu_obj <- FindNeighbors(merged_seu_obj, reduction = "harmony", dims = 1:40) %>% FindClusters(resolution = c(1:10 / 10))
  
  return(merged_seu_obj)
}

# 2. Find Markers
run.presto.specific.columns <- function(seu_obj, col_name_list = NULL, RNA_cols = F) {
  
  DefaultAssay(seu_obj) <- "RNA"
  seu_obj <- JoinLayers(seu_obj)
  
  if (RNA_cols) {
    col_name_list <- colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), "RNA_snn")]
  }
  
  # Get Wilcox AUC results
  wilcox.auc.list <- list()
  
  i <- 1
  for (col in col_name_list) {
    wilcox.auc.list[[i]] <- presto::wilcoxauc(seu_obj, group_by = col_name_list[i])
    i <- i + 1
  }
  
  names(wilcox.auc.list) <- col_name_list
  
  # Get Filered Wilcox AUC results
  filtered.markers.list <- lapply(wilcox.auc.list, FUN = function(x) {
    x <- x %>%
      dplyr::filter(auc > 0.5) %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::filter(pct_in > 50)
  })
  
  # Presto Top Markers
  top.markers.list <- lapply(wilcox.auc.list, FUN = function(x) {
    x <- presto::top_markers(x, n = 100, auc_min = 0.5, padj_max = 0.05, pct_in_min = 50)
  })
  
  # Final Results List
  final_results_list <- list(wilcox.auc.list, filtered.markers.list, top.markers.list)
  names(final_results_list) <- c("Raw-Markers", "Filtered-Markers", "Top-Markers")
  
  return(final_results_list)
}

# 2. Plot UMAP
plot.umaps <- function(seu_obj, save_path) {
  
  # Plot All Together in One File
  pdf(glue("{save_path}/all_resolutions_UMAPs.pdf"), height = 12, width = 16)
  for (res in colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), "RNA_snn")]) {
    res_num <- stringr::str_replace(string = res, pattern = "RNA_snn_res.", replacement = "")
    a <- SCpubr::do_DimPlot(seu_obj, group.by = res, plot.axes = T, plot.title = glue("UMAP Resolution {res_num}")) + theme(plot.title = element_text(hjust=0.5))
    a  %>% print()
  }
  dev.off()
  
  # Plot All Resolutions in Separate Files
  for (res in colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), "RNA_snn")]) {
    res_num <- stringr::str_replace(string = res, pattern = "RNA_snn_res.", replacement = "")
    
    pdf(glue("{save_path}/specific-resolutions/UMAP-res-{res_num}.pdf"), height = 12, width = 16)
    a <- SCpubr::do_DimPlot(seu_obj, group.by = res, plot.axes = T, plot.title = glue("UMAP Resolution {res_num}")) + theme(plot.title = element_text(hjust=0.5))
    a %>% print()
    dev.off()
  }
}

plot.umaps.labeled <- function(seu_obj, save_path) {
  
  # Plot All Together in One File
  pdf(glue("{save_path}/all_resolutions_UMAPs-Labeled.pdf"), height = 12, width = 16)
  for (res in colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), "RNA_snn")]) {
    res_num <- stringr::str_replace(string = res, pattern = "RNA_snn_res.", replacement = "")
    a <- SCpubr::do_DimPlot(seu_obj, group.by = res, plot.axes = T, plot.title = glue("UMAP Resolution {res_num}"), label = T) + theme(plot.title = element_text(hjust=0.5))
    a  %>% print()
  }
  dev.off()
  
  # Plot All Resolutions in Separate Files
  for (res in colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), "RNA_snn")]) {
    res_num <- stringr::str_replace(string = res, pattern = "RNA_snn_res.", replacement = "")
    
    pdf(glue("{save_path}/specific-resolutions/UMAP-res-{res_num}-Labeled.pdf"), height = 12, width = 16)
    a <- SCpubr::do_DimPlot(seu_obj, group.by = res, plot.axes = T, plot.title = glue("UMAP Resolution {res_num}"), label = T) + theme(plot.title = element_text(hjust=0.5))
    a %>% print()
    dev.off()
  }
}




