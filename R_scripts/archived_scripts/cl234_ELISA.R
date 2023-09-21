# #Figure 1D: MyD88-T6BM-3xA, WT, 3x KO-------------------------------------------------------------------------
chimeric_MyD88_3xA_data <- 
  All_days_data %>% 
  filter(
    Cohort %in% c("MyD88-T6BM-3xA", "MyD88-GFP/ TRAF6-mScarlet"),
    Date == "20220701"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-GFP/ TRAF6-mScarlet" ~ "MyD88-GFP/\n TRAF6-mScarlet",
      Cohort == "MyD88-T6BM-3xA" ~ "MyD88-T6BM-3xA"))

Triple_KO <-
  All_days_data %>% 
  filter(
    Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-" ~ "MyD88-/-/ IRAK4-/-/\n IRAK1-/-"))

chimeric_MyD88_3xA_data <- rbind(chimeric_MyD88_3xA_data, Triple_KO)


chimeric_MyD88_3xA_data$Cohort <-
  factor(chimeric_MyD88_3xA_data$Cohort,
         levels=c( "MyD88-GFP/\n TRAF6-mScarlet", 
                   "MyD88-/-/ IRAK4-/-/\n IRAK1-/-", 
                   "MyD88-T6BM-3xA"
         ))

chimeric_MyD88_3xA_test <-
  chimeric_MyD88_3xA_data %>% 
  mutate(
    IL2_concentration_Dilution_Factor_mean = IL2_concentration_Dilution_Factor
  ) %>% 
  as.data.table()

chimeric_MyD88_3xA_plot <- 
  chimeric_MyD88_3xA_data %>% 
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

chimeric_MyD88_3xA_stats <-
  chimeric_MyD88_3xA_data %>% 
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

Plot_Fx(chimeric_MyD88_3xA_plot, chimeric_MyD88_3xA_stats, chimeric_MyD88_3xA_test)

setwd("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 1")

ggsave(
  "cl234_ELISA.pdf",
  plot = last_plot(),
  scale = 3,
  units = "mm",
  height = 45,
  width = 50
)
