library(pacman)
pacman::p_load(data.table, tidyr, dplyr, ggplot2, interleave, ggpubr, stringr)
filter <- dplyr::filter

#provide location of template script
TEMPLATE_SCRIPT <- "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/ELISA analysis_template.R"

#Specify which data to analyse
Input_Directory_List <- 
  c("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20220609_Elisa",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20220623_Elisa",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20220701_Elisa"
  )

# Run the setup -----------------------------------------------------------
All_days_data <- data.frame()

for (Input_Directory in Input_Directory_List){
  
  print(Input_Directory)
  source(TEMPLATE_SCRIPT, local = T)
  All_plates_data$Date <- (strsplit(strsplit(Input_Directory, "/")[[1]][7], "_"))[[1]][1]
  All_days_data <- rbind(All_days_data, All_plates_data)
  
  rm(
    All_plates_data
  )
  
}

All_days_data$Stimulation_Condition <- factor(All_days_data$Stimulation_Condition, levels=c("Unstimulated", "Stimulated"))

#Standard error of the mean function
sem <- function(x) sd(x)/sqrt(length(x))
