library(pacman)
pacman::p_load(data.table, tidyr, dplyr, ggplot2, interleave, ggpubr)

#Specify Dilution factor of supernatant
Dilution_Factor <- 5

#Creating Output Folder
Output_Directory <- file.path(Input_Directory, "Output")
if(!file.exists(Output_Directory)){
  dir.create(Output_Directory)
}

#Calculating Number of plates
for (x in 1:10){
  if(file.exists(paste0(Input_Directory, "/Plate_", x))){
    Number_of_plates <- x
  }
}

ELISA_Fx <- function(plate_number){
  
  # Reading Data ------------------------------------------------------------
  #Getting Path of treatment conditions and values
  Input_plate <- paste0(Input_Directory, "/Plate_", plate_number)
  
  #Reading Plate Treatement 
  Values_Measured <- fread(paste0(Input_plate, "/Values_Measured.csv"), header = F)
  Cohort <- fread(paste0(Input_plate, "/Cohort.csv"), header = F)
  Stimulation_Condition <- fread(paste0(Input_plate, "/Stimulation_Condition.csv"), header = F)
  Sample_Day <- fread(paste0(Input_plate, "/Sample_Day.csv"), header = F)
  
  
  #Converting tables into vector for to make a single table
  Values_Measured <- as.vector(as.matrix(Values_Measured))
  Cohort <- as.vector(as.matrix(Cohort))
  Stimulation_Condition <- as.vector(as.matrix(Stimulation_Condition))
  Sample_Day <- as.vector(as.matrix(Sample_Day))
  
  #Creating Table containing all plate Information
  Plate <- NULL
  Plate$Values_Measured <- Values_Measured
  Plate$Cohort <- Cohort
  Plate$Stimulation_Condition <- Stimulation_Condition
  Plate$Sample_Day <- Sample_Day
  
  rm(
    Values_Measured,
    Cohort,
    Stimulation_Condition,
    Sample_Day
  )
  
  Plate <- Plate %>% as.data.table()
  
  #Removing Empty Wells
  Plate <- Plate %>% 
    filter(
      Cohort != "Blank"
    ) %>% as.data.table()
  
  
  # Standard Curve ----------------------------------------------------------
  #Creating a Standard Curve
  Plate_Standards <- Plate %>% 
    filter(
      Stimulation_Condition == "Calibration"
    ) %>% 
    group_by(
      Cohort
    ) %>% 
    summarise(
      Values_Measured_mean = mean(Values_Measured)
      # Values_Measured_median = median(Values_Measured)
    ) %>%  
    mutate(
      Cohort = as.numeric(Cohort)
    ) %>% 
    arrange(
      Cohort
    )
  
  Fit <- lm(Cohort ~ Values_Measured_mean -1, data = Plate_Standards) #linear model of the Standard curve. -1 omits the intercept
  
  R <- summary(Fit)$r.squared
  
  Rsquare <- signif(R, digits = 4)
  
  rm(
    R
  )
  
  print(paste0("IL2-Amount = slope*Intensity"))
  print(paste0("IL2-Amount = ", Fit$coefficients[1],"*Intensity"))
  
  Plate_Standards <- 
    Plate_Standards %>% 
    mutate(
      Fit_Test = (Fit$coefficients[1]*Values_Measured_mean)
    )
  
  ggplot(
    data = Plate_Standards,
  ) +
    geom_point(
      aes(
        x = Values_Measured_mean,
        y = Cohort
      ),
      size = 5
    ) +
    geom_line(
      aes(
        x = Values_Measured_mean,
        y = Fit_Test
      ),
      linetype = "dashed"
    ) +
    annotate(
      'text',
      x = 0.15,
      y = 700,
      label = paste0("R^2 = ",Rsquare),
      size = 10
    ) +
    annotate(
      'text',
      x = max(Plate_Standards$Values_Measured_mean) - (0.25*max(Plate_Standards$Values_Measured_mean)),
      y = 150,
      label = paste0("IL-Amount = \n", signif(Fit$coefficients[1], digits = 4),"*Intensity")
    ) +
    labs(
      x = "Measured Values",
      y = "IL-Concentration (pg/mL)"
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 20)
    )
  
  Save_Name <- paste0("Plate_Number_" , plate_number, "_Standard_Curve.pdf")
  Save_Name <- file.path(Output_Directory, Save_Name)
  
  ggsave(
    Save_Name,
    plot = last_plot(),
    height = 3*3,
    width = 5*4
  )
  
  rm(
    Save_Name,
    Plate_Standards,
    Rsquare
  )
  
  
  # Fitting Data To Standarad Curve -----------------------------------------
  Plate <- Plate %>% 
    filter(
      Stimulation_Condition != "Calibration"
    ) %>% 
    mutate(
      Values_Measured = as.numeric(Values_Measured),
      IL2_concentration = (Fit$coefficients[1]*Values_Measured),
      IL2_concentration_Dilution_Factor = IL2_concentration*Dilution_Factor
    )
  
  rm(
    Fit
  )
  
  return(Plate)
}

All_plates_data = data.frame()

for (plate_number in 1:Number_of_plates){
  

  data_temp <- ELISA_Fx(plate_number)
  data_temp$Plate <- plate_number
  
  All_plates_data <- rbind(All_plates_data, data_temp)
  
  rm(data_temp)
}