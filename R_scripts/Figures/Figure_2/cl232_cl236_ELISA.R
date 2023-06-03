# #Figure 2: DHF91-TIR-T6BM, WT, 3x KO-------------------------------------------------------------------------

DHF91_TIR_T6BM_data <- 
  All_days_data %>% 
  filter(
    Cohort %in% c("DHF91-TIR-T6BM low", "MyD88-T6BM"),
    Date == "20220623"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-T6BM" ~ "MyD88-T6BM",
      Cohort == "DHF91-TIR-T6BM low" ~ "DHF91-TIR-T6BM"))

Triple_KO <-
  All_days_data %>% 
  filter(
    Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-" ~ "MyD88-/-/ IRAK4-/-/\n IRAK1-/-"))

DHF91_TIR_T6BM_data <- rbind(DHF91_TIR_T6BM_data, Triple_KO)


DHF91_TIR_T6BM_data$Cohort <-
  factor(DHF91_TIR_T6BM_data$Cohort,
         levels=c( "MyD88-/-/ IRAK4-/-/\n IRAK1-/-", 
                   "MyD88-T6BM", 
                   "DHF91-TIR-T6BM"
         ))

DHF91_TIR_T6BM_test <-
  DHF91_TIR_T6BM_data %>% 
  mutate(
    IL2_concentration_Dilution_Factor_mean = IL2_concentration_Dilution_Factor
  ) %>% 
  as.data.table()

DHF91_TIR_T6BM_plot <- 
  DHF91_TIR_T6BM_data %>% 
  group_by(
    Cohort,
    Stimulation_Condition,
    Sample_Day
  ) %>% 
  summarise(
    IL2_concentration_Dilution_Factor_mean = mean(IL2_concentration_Dilution_Factor),
    IL2_concentration_Dilution_Factor_median = median(IL2_concentration_Dilution_Factor)
  ) %>% 
  as.data.table()

DHF91_TIR_T6BM_stats <-
  DHF91_TIR_T6BM_data %>% 
  group_by(
    Cohort,
    Stimulation_Condition
  ) %>% 
  summarise(
    IL2_concentration_Dilution_Factor_mean = mean(IL2_concentration_Dilution_Factor),
    IL2_concentration_Dilution_Factor_median = median(IL2_concentration_Dilution_Factor),
    IL2_concentration_Dilution_Factor_sd = sd(IL2_concentration_Dilution_Factor)
  ) %>% 
  as.data.table()

Plot_Fx(DHF91_TIR_T6BM_plot, DHF91_TIR_T6BM_stats, DHF91_TIR_T6BM_test)

setwd("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 2")

ggsave(
  "cl232_cl236_ELISA.pdf",
  plot = last_plot(),
  scale = 3,
  units = "mm",
  height = 45,
  width = 50
)
