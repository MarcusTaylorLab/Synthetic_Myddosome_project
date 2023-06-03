# #Figure 3A: WT, 3xKO, MyD88-GFP-synTRAF6-BD-3x, MyD88-GFP-synTRAF6-BD-3xA--------------------------------------------------------------------------
synTRAF6_BD_3x_data <- 
  All_days_data %>% 
  filter(
    Cohort %in% c("MyD88-GFP-synTRAF6-BD-3x"),
    Date == "20221208",
    Stimulation_Condition %in% c("0e+00", "1e+01")
  ) %>% 
  mutate(
    Stimulation_Condition = case_when(
      Stimulation_Condition == "0e+00" ~ "Unstimulated",
      Stimulation_Condition == "1e+01" ~ "Stimulated"
    ),
    Cohort = case_when(
      Cohort == "MyD88-GFP-synTRAF6-BD-3x" ~ "MyD88-sT6BM-3x"
    )
  )

chimeric_MyD88_data <- 
  All_days_data %>% 
  filter(
    Cohort %in% c("MyD88-T6BM"),
    Date == "20220623"
  )

synTRAF6_BD_3xA_data <- 
  All_days_data %>% 
  filter(
    Cohort %in% c("MyD88-GFP-synTRAF6-BD-3xA"),
    Date == "20221124"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-GFP-synTRAF6-BD-3xA" ~ "MyD88-sT6BM-3xA"
    )
  )

Triple_KO <-
  All_days_data %>% 
  filter(
    Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-" ~ "MyD88-/-/ IRAK4-/-/\n IRAK1-/-"))

synTRAF6_BD_3x_data <- rbind(synTRAF6_BD_3x_data, synTRAF6_BD_3xA_data, Triple_KO, chimeric_MyD88_data)

synTRAF6_BD_3x_data$Stimulation_Condition <- factor(synTRAF6_BD_3x_data$Stimulation_Condition, levels = c("Unstimulated", "Stimulated"))

synTRAF6_BD_3x_data$Cohort <-
  factor(synTRAF6_BD_3x_data$Cohort,
         levels=c( "MyD88-T6BM", 
                   "MyD88-/-/ IRAK4-/-/\n IRAK1-/-",
                   "MyD88-sT6BM-3x", 
                   "MyD88-sT6BM-3xA"
         ))

synTRAF6_BD_3x_test <- 
  synTRAF6_BD_3x_data %>%
  mutate(
    IL2_concentration_Dilution_Factor_mean = IL2_concentration_Dilution_Factor
  ) %>% 
  as.data.table()

synTRAF6_BD_3x_plot <- 
  synTRAF6_BD_3x_data %>% 
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

synTRAF6_BD_3x_stats <-
  synTRAF6_BD_3x_data %>% 
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

Plot_Fx(synTRAF6_BD_3x_plot, synTRAF6_BD_3x_stats, synTRAF6_BD_3x_test
        )

setwd("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 3")

ggsave(
  "cl232_cl255_cl257_ELISA.pdf",
  plot = last_plot(),
  scale = 3,
  units = "mm",
  height = 45,
  width = 50
)
