data_tay_conversion <- function(fpath){
  loc <- strsplit(fpath, split = "//data-tay/TAYLOR-LAB/")[[1]][2]
  fpath_new <- paste("/Volumes/TAYLOR-LAB/", loc, sep = "")
  
  return(fpath_new)
}

Table_1 <- lapply(TablePaths, data_tay_conversion)
