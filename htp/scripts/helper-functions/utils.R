# utils.R

# 1. Create New Save Directories
create.dir <- function(path) {
  if (!dir.exists(path)) {dir.create(path, recursive = T)}
}

# 2. Load Multi Sheets Excel
load.xlsx.multi.sheet <- function(path) {
  sheet_names <- readxl::excel_sheets(path)
  sheets <- lapply(sheet_names, function(x) {readxl::read_xlsx(path, sheet = x)})
  names(sheets) <- sheet_names
  
  return(sheets)
}
