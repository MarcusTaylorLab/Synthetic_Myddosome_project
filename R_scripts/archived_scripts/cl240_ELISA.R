# #Figure 1C: MyD88-TIR-T6BM, WT, 3x KO-------------------------------------------------------------------------
MyD88_TIR_data  <- 
  All_days_data %>% 
  filter(
    Cohort %in% c("MyD88-TIR-T6BM low", "MyD88-TIR-T6BM high", "MyD88-GFP/ TRAF6-mScarlet"),
    Date == "20220623"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-GFP/ TRAF6-mScarlet" ~ "MyD88-GFP/\n TRAF6-mScarlet",
      Cohort == "MyD88-TIR-T6BM low" ~ "TIR-T6BM low",
      Cohort == "MyD88-TIR-T6BM high" ~ "TIR-T6BM high"))

MyD88_TIR_data <- rbind(MyD88_TIR_data, Triple_KO)

MyD88_TIR_data$Cohort <-
  factor(MyD88_TIR_data$Cohort,
         levels=c( "MyD88-GFP/\n TRAF6-mScarlet", 
                   "MyD88-/-/ IRAK4-/-/\n IRAK1-/-", 
                   "TIR-T6BM low",
                   "TIR-T6BM high"
         ))

MyD88_TIR_test <-
  MyD88_TIR_data %>% 
  mutate(
    IL2_concentration_Dilution_Factor_mean = IL2_concentration_Dilution_Factor
  ) %>% 
  as.data.table()

MyD88_TIR_plot <- 
  MyD88_TIR_data %>% 
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

MyD88_TIR_stats <-
  MyD88_TIR_data %>% 
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

Plot_Fx(MyD88_TIR_plot, MyD88_TIR_stats, MyD88_TIR_test)

setwd("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 1")

ggsave(
  "cl240_ELISA.pdf",
  plot = last_plot(),
  scale = 3,
  units = "mm",
  height = 45,
  width = 50
)