library(pacman)

pacman::p_load(ggplot2, data.table, dplyr, ggfx, viridis, RColorBrewer, ggpubr, lemon, ggbreak, tidyfast, ggbeeswarm, R.utils)
filter <- dplyr::filter

#running script on mac?
mac = TRUE

#the file structure differs between mac and windows. This function converts the windows structure to the mac structure
data_tay_conversion <- function(fpath){
  loc <- strsplit(fpath, split = "//data-tay/TAYLOR-LAB/")[[1]][2]
  fpath_new <- paste("/Volumes/TAYLOR-LAB/", loc, sep = "")
  
  return(fpath_new)
}

#Set path of analyzed Images
TablePaths <- c(
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88 TRAF6/20230405 1.5nM_cl069_TRAF6_MyD88 001/Essential.csv.gz",
  #"//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88 TRAF6/20230405 1.5nM_cl069_TRAF6_MyD88 006/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88 TRAF6/20230413 2.5nM_cl069_TRAF6_MyD88 001/Essential.csv.gz",
  
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TRAF6-BD TRAF6/20220516 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TRAF6-BD TRAF6/20220610 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TRAF6-BD TRAF6/20220615 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz",
  
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20220610 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 007/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20230405 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20230413 2.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz"
)

if(mac == TRUE){
  TablePaths <- lapply(TablePaths, data_tay_conversion)
}
  
Table <- lapply(TablePaths, 
                  fread)
  
#fill=FALSE since all images are two color
Table <- rbindlist(Table, fill = FALSE)
  
#Create the same order for all the following plots
Table$COHORT <- factor(
  Table$COHORT, levels = c("MyD88 TRAF6",
                           "MyD88-TRAF6-BD TRAF6",
                           "MyD88-TIR-TRAF6-BD TRAF6"
  )
)

#set a color palette for the plots
color_violin<-c(
  "MyD88-TRAF6-BD TRAF6" = "#117733", 
  "MyD88 TRAF6" = "#44AA99",
  "MyD88-TIR-TRAF6-BD TRAF6" = "#332288"
)

#Standard error of the mean function
sem <- function(x) sd(x)/sqrt(length(x))

#create a list for comparisons
my_comparison <- 
  list(
    c("MyD88 TRAF6", "MyD88-TRAF6-BD TRAF6"),
    c("MyD88 TRAF6", "MyD88-TIR-TRAF6-BD TRAF6"),
    c("MyD88-TRAF6-BD TRAF6", "MyD88-TIR-TRAF6-BD TRAF6")
)
