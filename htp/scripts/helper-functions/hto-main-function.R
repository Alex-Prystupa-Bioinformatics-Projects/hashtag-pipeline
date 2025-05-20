source("helper-functions/hto-utils.R")
source("helper-functions/utils.R")

run.hto.workflow <- function(path, project_name, samples_groups_df) {
  
  # 1. Load Data
  seu_obj <- load.hto.data(path, project_name, samples_groups_df)
  
  # 2. Downstream Analysis + HTO Demux
  seu_obj <- normalize.hto.demux(seu_obj)
}
