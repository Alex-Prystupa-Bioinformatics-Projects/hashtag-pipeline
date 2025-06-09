# hto-utils.R

# 1. Load Data
load.hto.data <- function(path, project_name) {

  seu_data <- Read10X(data.dir = path)
  seu_obj <- CreateSeuratObject(counts = seu_data$`Gene Expression`, project = project_name, min.cells = 3, min.features = 200)
  
  print("Object Created")
  # 2. Subset Counts to Contain only barcodes in both Gene Expression & HTO
  joint_bcs <- intersect(rownames(seu_obj@meta.data), colnames(seu_data$`Cell Hashing`))
  
  seu_obj <- seu_obj[ ,joint_bcs]
  seu_data$`Cell Hashing` <- seu_data$`Cell Hashing`[ ,joint_bcs]
  
  # 3. Add HTO Assay
  seu_obj[["HTO"]] <- CreateAssayObject(counts = seu_data$`Cell Hashing`)

  return(seu_obj)
}

# 2. Downstream Analysis + HTO Demux
normalize.hto.demux <- function(seu_obj) {
  
  # Normalize RNA
  seu_obj <- NormalizeData(seu_obj)
  
  # Find and scale variable features
  seu_obj <- FindVariableFeatures(seu_obj, selection.method = "mean.var.plot")
  seu_obj <- ScaleData(seu_obj, features = VariableFeatures(seu_obj))
  
  # Normalize HTO
  seu_obj <- NormalizeData(seu_obj, assay = "HTO", normalization.method = "CLR")
  
  # HTO Demux
  seu_obj <- HTODemux(seu_obj, assay = "HTO", positive.quantile = 0.99)
  
  return(seu_obj)
}

# 3. Add Meta Data
seu.add.meta.data <- function(seu_obj, samples_groups_df) {
  
  # 1. Subset by Hashtags in Samples groups
  seu_obj <- subset(seu_obj, subset = HTO_maxID %in% samples_groups_df$Hashtag)
  
  # 2. Join with Meta Data
  seu_obj@meta.data <- seu_obj@meta.data %>% tibble::rownames_to_column("barcodes")
  seu_obj@meta.data <- left_join(seu_obj@meta.data, samples_groups_df, by = c("HTO_maxID" = "Hashtag"))
  seu_obj@meta.data <- seu_obj@meta.data %>% tibble::column_to_rownames("barcodes")
  
  return(seu_obj)
}

# 4. Full HTO Workflow
run.hto.workflow <- function(path, project_name, samples_groups_df) {
  
  # 1. Load Data
  seu_obj <- load.hto.data(path, project_name)
  
  # 2. Downstream Analysis + HTO Demux
  seu_obj <- normalize.hto.demux(seu_obj)
  
  # 3. Add Meta Data
  seu_obj <- seu.add.meta.data(seu_obj, samples_groups_df)
  
  return(seu_obj)
}

# 5. Get HTO Demux Stats
get.hto.demux.stats <- function(seu_obj_list) {
  
  hto_stats_list <- list(); i <- 1
  
  for (name in names(seu_obj_list)) {
    hto_stats_list[[i]] <- table(seu_obj_list[[name]]@meta.data$HTO_classification.global); i <- i + 1
  }
  
  names(hto_stats_list) <- names(seu_obj_list)
  
  hto_demux_stats <- data.frame(do.call(rbind, hto_stats_list))
  
  return(hto_demux_stats)
}

# 6. Get Filtering Stats per Hashtag
filter.stats <- function(seu_obj, seu_obj_singlet, hto_run) {
  
  pre_filter_table <- table(seu_obj_list[[hto_run]]$HTO_maxID, seu_obj_list[[hto_run]]$orig.ident)
  post_filter_table <- table(seu_obj_list_singlets[[hto_run]]$HTO_maxID, seu_obj_list_singlets[[hto_run]]$orig.ident)
  
  stats_df <- cbind(pre_filter_table, post_filter_table) %>% 
    magrittr::set_colnames(c(glue("{hto_run}_all_Freq"), glue("{hto_run}_singlets_Freq")))
  
  return(stats_df)
}

# 7. Plot Feature Scatter Plots
plot.feature_scatter <- function(seu_obj) {
  
  hashtag.combs <- combn(unique(seu_obj$HTO_maxID), m=2)
  
  for (i in 1:ncol(hashtag.combs)) {
    FeatureScatter(seu_obj, feature1 = hashtag.combs[1,i], feature2 = hashtag.combs[2,i]) %>% print()
  }
}


