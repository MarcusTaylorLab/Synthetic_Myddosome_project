library(pacman)

pacman::p_load(ggplot2, data.table, dplyr, ggfx, viridis, RColorBrewer, ggpubr, lemon, ggbreak, tidyfast, ggbeeswarm, R.utils, signal, ggh4x)
filter <- dplyr::filter

#Set path of analyzed Images
TablePaths <- c(
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88 TRAF6/20230405 1.5nM_cl069_TRAF6_MyD88 001/Essential.csv.gz",
  #"/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88 TRAF6/20230405 1.5nM_cl069_TRAF6_MyD88 006/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88 TRAF6/20230413 2.5nM_cl069_TRAF6_MyD88 001/Essential.csv.gz",
  
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TRAF6-BD TRAF6/20220516 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TRAF6-BD TRAF6/20220610 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TRAF6-BD TRAF6/20220615 1.5nM_cl232_TRAF6_MyD88-TRAF6-BD-GFP 002/Essential.csv.gz",
  
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20220610 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20220615 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 007/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20230405 1.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-TIR-TRAF6-BD TRAF6/20230413 2.5nM_cl240_TRAF6_MyD88-TIR-TRAF6-BD-GFP 001/Essential.csv.gz",
  
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-DHF91-TRAF6-BD TRAF6/20220615 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-DHF91-TRAF6-BD TRAF6/20221207 4nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-DHF91-TRAF6-BD TRAF6/20230405 1.5nM_cl236_TRAF6_MyD88-DHF91-TRAF6-BD-GFP 003/Essential.csv.gz",
  
  "/Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/20230602_batch/Output/BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230602 6nM_cl321-BDLD57H_TRAF6_MyD88_002/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/20230616_batch/Output/BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230616 3nM_cl321-BDLD57H_TRAF6_MyD88 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/20230623_batch/Output/BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6/20230623 3nM_cl321-BDLD57H_TRAF6_MyD88 002/Essential.csv.gz"
)

Table <- lapply(TablePaths, 
                fread)

#fill=FALSE since all images are two color
Table <- rbindlist(Table, fill = FALSE)

#Create the same order for all the following plots
Table$COHORT <- factor(
  Table$COHORT, levels = rev(c("MyD88 TRAF6",
                           "MyD88-TRAF6-BD TRAF6",
                           "MyD88-TIR-TRAF6-BD TRAF6",
                           "BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6",
                           "MyD88-DHF91-TRAF6-BD TRAF6"
                           
  )
))

#set a color palette for the plots
color_violin<-c(
  "MyD88-TRAF6-BD TRAF6" = "#117733", 
  "MyD88 TRAF6" = "#333333",
  "MyD88-TIR-TRAF6-BD TRAF6" = "#332288",
  "MyD88-DHF91-TRAF6-BD TRAF6" = "#44AA99",
  "BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6" = "#AA4499"
)

#create a list for comparisons
my_comparison <- list(
  c("MyD88 TRAF6", "MyD88-TRAF6-BD TRAF6"),
  c("MyD88 TRAF6", "MyD88-TIR-TRAF6-BD TRAF6"),
  c("MyD88 TRAF6", "MyD88-DHF91-TRAF6-BD TRAF6"),
  c("MyD88 TRAF6", "BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6"),
  c("MyD88-TRAF6-BD TRAF6", "MyD88-TIR-TRAF6-BD TRAF6"),
  c("MyD88-TRAF6-BD TRAF6", "MyD88-DHF91-TRAF6-BD TRAF6"),
  c("MyD88-TRAF6-BD TRAF6", "BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6"),
  c("MyD88-TIR-TRAF6-BD TRAF6", "MyD88-DHF91-TRAF6-BD TRAF6"),
  c("MyD88-TIR-TRAF6-BD TRAF6", "BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6"),
  c("MyD88-DHF91-TRAF6-BD TRAF6", "BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6")
)

#Standard error of the mean function
sem <- function(x) sd(x)/sqrt(length(x))

