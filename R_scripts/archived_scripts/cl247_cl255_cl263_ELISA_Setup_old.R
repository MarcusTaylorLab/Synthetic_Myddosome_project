library(pacman)
pacman::p_load(data.table, tidyr, dplyr, ggplot2, interleave, ggpubr, stringr, ggbreak)
filter <- dplyr::filter

#provide location of template script
TEMPLATE_SCRIPT <- "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/ELISA analysis_template.R"

#Specify which data to analyse
Input_Directory_List <- 
  c("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20220609_Elisa",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20220623_Elisa",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20221006_Elisa_doseresponse",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20221124_Elisa_doseresponse",
    "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/ELISA data/ELISA analysis in R/20221208_Elisa_doseresponse")

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

#Define the plotting function
Plot_Fx <- function(chimeric_MyD88_plot, chimeric_MyD88_stats, chimeric_MyD88_test){  
  
  color_elisa <- c("Unstimulated" = "white",
                   "Stimulated" = "grey")
  
  ggplot(
    data = chimeric_MyD88_stats,
    aes(
      x = Cohort,
      y = IL2_concentration_Dilution_Factor_mean,
      fill = Stimulation_Condition
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.5),
      color = "black",
      width = 0.5
    ) +
    geom_errorbar(
      data = chimeric_MyD88_stats,
      aes(
        x = Cohort,
        ymin = IL2_concentration_Dilution_Factor_mean - IL2_concentration_Dilution_Factor_sd,
        ymax = IL2_concentration_Dilution_Factor_mean + IL2_concentration_Dilution_Factor_sd
      ),
      linewidth = .75,
      position = position_dodge(width = 0.5),
      width = 0.25
    ) +
    geom_point(
      data = chimeric_MyD88_plot,
      aes(
        x = Cohort,
        y = IL2_concentration_Dilution_Factor_mean
      ),
      size = 2,
      position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4)
    )+
    stat_compare_means(
      data = chimeric_MyD88_test,
      method = "wilcox.test",
      label = "p.signif"
    )+
    scale_y_continuous(
      breaks = scales::breaks_width(1000)
    )+
    fill_palette(
      palette = color_elisa
    ) +
    labs(
      y = "IL-2 (pg/mL)"
    ) +
    theme_classic(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.title.x = element_blank(),
      axis.text.y = element_text(color = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(color = "black")
    )
}
