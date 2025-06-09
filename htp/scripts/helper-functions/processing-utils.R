# processing-utils.R

# 1. Multi Sample Pipeline
multi_sample_pipeline <- function(merged_seu_obj, split_col = "orig.ident", res_cols = "RNA") {
  
  if (res_cols == "SCT") {
    merged_seu_obj <- SCTransform(merged_seu_obj, verbose = FALSE)
  } else {
    merged_seu_obj <- NormalizeData(merged_seu_obj)
    merged_seu_obj <- FindVariableFeatures(merged_seu_obj)
    merged_seu_obj <- ScaleData(merged_seu_obj)
  }
  
  merged_seu_obj <- RunPCA(merged_seu_obj)
  merged_seu_obj <- RunHarmony(merged_seu_obj, assay.use = res_cols, group.by.vars = split_col)
  merged_seu_obj <- RunUMAP(merged_seu_obj, reduction = "harmony", dims = 1:40)
  merged_seu_obj <- FindNeighbors(merged_seu_obj, reduction = "harmony", dims = 1:40)
  merged_seu_obj <- FindClusters(merged_seu_obj, resolution = c(1:10 / 10))
  
  return(merged_seu_obj)
}

# 2. Find Markers
run.presto.specific.columns <- function(seu_obj, col_name_list = NULL, res_cols = "RNA") {
  
  DefaultAssay(seu_obj) <- "RNA"
  
  if (res_cols == "RNA") {
    seu_obj <- JoinLayers(seu_obj)
    col_name_list <- colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), "RNA_snn")]
  }

  if (res_cols == "SCT") {

    # RNA Preprocessing (I like using RNA for DE Genes)
    seu_obj <- NormalizeData(seu_obj)
    seu_obj <- FindVariableFeatures(seu_obj)
    seu_obj <- ScaleData(seu_obj)
    seu_obj <- JoinLayers(seu_obj)

    col_name_list <- colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), "SCT_snn")]
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
plot.umaps <- function(seu_obj, save_path, res_cols = "RNA") {
  
  prefix <- if (res_cols == "SCT") "SCT_snn" else "RNA_snn"
  pattern <- if (res_cols == "SCT") "SCT_snn_res." else "RNA_snn_res."
  
  # Plot All Together in One File
  pdf(glue("{save_path}/all_resolutions_UMAPs.pdf"), height = 12, width = 16)
  for (res in colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), prefix)]) {
    res_num <- stringr::str_replace(string = res, pattern = pattern, replacement = "")
    a <- SCpubr::do_DimPlot(seu_obj, group.by = res, plot.axes = TRUE, plot.title = glue("UMAP Resolution {res_num}")) + 
      theme(plot.title = element_text(hjust = 0.5))
    print(a)
  }
  dev.off()
  
  # Plot All Resolutions in Separate Files
  for (res in colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), prefix)]) {
    res_num <- stringr::str_replace(string = res, pattern = pattern, replacement = "")
    
    pdf(glue("{save_path}/specific-resolutions/UMAP-res-{res_num}.pdf"), height = 12, width = 16)
    a <- SCpubr::do_DimPlot(seu_obj, group.by = res, plot.axes = TRUE, plot.title = glue("UMAP Resolution {res_num}")) + 
      theme(plot.title = element_text(hjust = 0.5))
    print(a)
    dev.off()
  }
}

plot.umaps.labeled <- function(seu_obj, save_path, res_cols = "RNA") {
  
  prefix <- if (res_cols == "SCT") "SCT_snn" else "RNA_snn"
  pattern <- if (res_cols == "SCT") "SCT_snn_res." else "RNA_snn_res."
  
  # Plot All Together in One File
  pdf(glue("{save_path}/all_resolutions_UMAPs-Labeled.pdf"), height = 12, width = 16)
  for (res in colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), prefix)]) {
    res_num <- stringr::str_replace(string = res, pattern = pattern, replacement = "")
    a <- SCpubr::do_DimPlot(seu_obj, group.by = res, plot.axes = TRUE, plot.title = glue("UMAP Resolution {res_num}"), label = TRUE) + 
      theme(plot.title = element_text(hjust = 0.5))
    print(a)
  }
  dev.off()
  
  # Plot All Resolutions in Separate Files
  for (res in colnames(seu_obj@meta.data)[startsWith(colnames(seu_obj@meta.data), prefix)]) {
    res_num <- stringr::str_replace(string = res, pattern = pattern, replacement = "")
    
    pdf(glue("{save_path}/specific-resolutions/UMAP-res-{res_num}-Labeled.pdf"), height = 12, width = 16)
    a <- SCpubr::do_DimPlot(seu_obj, group.by = res, plot.axes = TRUE, plot.title = glue("UMAP Resolution {res_num}"), label = TRUE) + 
      theme(plot.title = element_text(hjust = 0.5))
    print(a)
    dev.off()
  }
}




