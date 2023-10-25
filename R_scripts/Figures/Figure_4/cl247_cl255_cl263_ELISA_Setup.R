library(pacman)
pacman::p_load(data.table, tidyr, dplyr, ggplot2, interleave, ggpubr, stringr, ggbreak)
filter <- dplyr::filter

#Specify which data to analyse
Table_path <- 
  c("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20220609_Elisa/Output/ELISA_Table.csv",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20220623_Elisa/Output/ELISA_Table.csv",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20221006_Elisa_doseresponse/Output/ELISA_Table.csv",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20221124_Elisa_doseresponse/Output/ELISA_Table.csv",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20221208_Elisa_doseresponse/Output/ELISA_Table.csv"
    )

All_days_data <- data.frame()

for (Table in Table_path){
  
  All_plates_data <- fread(Table)
  All_plates_data$Date <- (strsplit(strsplit(Table, "/")[[1]][7], "_"))[[1]][1]
  All_days_data <- rbind(All_days_data, All_plates_data)
  
  rm(
    All_plates_data
  )
  
}

#Standard error of the mean function
sem <- function(x) sd(x)/sqrt(length(x))