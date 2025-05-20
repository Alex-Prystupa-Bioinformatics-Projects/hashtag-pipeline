# rna-utils.R

# 1. Plot QC Vln Plots
plot.rna.qc.vln.plots <- function(seu_obj) {
  Seurat::VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) %>% print()
}

# 2. Plot QC Scatter Plot
plot.rna.qc.scatter.plot <- function(seu_obj) {
  
  plot1 <- Seurat::FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- Seurat::FeatureScatter(seu_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  plot <- plot1 + plot2
  plot %>% print()
}

# 3. Get RNA QC Stats
get.rna.qc.summary <- function(seu_obj) {
  
  # 1. nFeature low cutoff
  df_low <- seu_obj@meta.data %>% select(nFeature_RNA, HTO_maxID) %>% 
    mutate(nFeature_RNA_low_pass = ifelse(nFeature_RNA > 500, "nFeature_RNA_low_pass", "nFeature_RNA_low_fail"))
  
  df_low_summary <- df_low %>%
    group_by(HTO_maxID, nFeature_RNA_low_pass) %>%
    summarize(count = n()) %>%
    spread(key = nFeature_RNA_low_pass, value = count, fill = 0)
  
  # 2. nFeature high cutoff
  df_high <- seu_obj@meta.data %>% select(nFeature_RNA, HTO_maxID) %>% 
    mutate(nFeature_RNA_high_pass = ifelse(nFeature_RNA < 7000, "nFeature_RNA_high_pass", "nFeature_RNA_high_fail"))
  
  df_high_summary <- df_high %>%
    group_by(HTO_maxID, nFeature_RNA_high_pass) %>%
    summarize(count = n()) %>%
    spread(key = nFeature_RNA_high_pass, value = count, fill = 0)
  
  # 3. percent.mt cutoff
  df_pct_mt <- seu_obj@meta.data %>% select(percent.mt, HTO_maxID) %>% 
    mutate(pct_mt_pass = ifelse(percent.mt < 20, "pct_mt_pass", "pct_mt_fail"))
  
  df_pct_mt_summary <- df_pct_mt %>%
    group_by(HTO_maxID, pct_mt_pass) %>%
    summarize(count = n()) %>%
    spread(key = pct_mt_pass, value = count, fill = 0) 
  
  # 4. Total Counts
  df <- seu_obj@meta.data %>% select(nFeature_RNA, percent.mt, HTO_maxID) %>% 
    mutate(final_pass = ifelse(nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 20, 
                               "qc_pass_count", "qc_fail_count"))
  
  df_summary <- df %>%
    group_by(HTO_maxID, final_pass) %>%
    summarize(count = n()) %>% 
    spread(key = final_pass, value = count, fill = 0)
  
  df <- dplyr::left_join(df_low_summary, df_high_summary, by = "HTO_maxID") %>% 
    dplyr::left_join(y = df_pct_mt_summary, by = "HTO_maxID") %>%
    dplyr::left_join(y = df_summary, by = "HTO_maxID")
  
  return(df)
}



