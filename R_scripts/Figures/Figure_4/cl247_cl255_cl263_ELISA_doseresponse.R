# #Figure 3A: WT, 3xKO, MyD88-GFP-synTRAF6-BD-3x, MyD88-GFP-synTRAF6-BD-3xA--------------------------------------------------------------------------
Dose_data <- 
  All_days_data %>% 
  filter(
    Cohort %in% c("MyD88-GFP-synTRAF6-BD-1x", "MyD88-GFP-synTRAF6-BD-3x", "MyD88-GFP-synTRAF6-BD-5x"),
    Date == "20221208"
  ) %>% 
  group_by(
    Sample_Day
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-GFP-synTRAF6-BD-1x" ~ "MyD88-sT6BM-1x",
      Cohort == "MyD88-GFP-synTRAF6-BD-3x" ~ "MyD88-sT6BM-3x",
      Cohort == "MyD88-GFP-synTRAF6-BD-5x" ~ "MyD88-sT6BM-5x"
    ),
    Relative_IL2_concentration = (Values_Measured-min(Values_Measured))/(max(Values_Measured)-min(Values_Measured)),
    Stimulation_Condition = as.numeric(Stimulation_Condition)
  )
  
Dose_plot <- 
  Dose_data %>% 
  group_by(
    Cohort,
    Stimulation_Condition,
    Sample_Day
  ) %>% 
  summarise(
    Relative_IL2_concentration_mean = mean(Relative_IL2_concentration),
    IL2_concentration_Dilution_Factor_mean = mean(IL2_concentration_Dilution_Factor)
  ) %>% 
  as.data.table()

Dose_stats <-
  Dose_data %>% 
  group_by(
    Cohort,
    Stimulation_Condition
  ) %>% 
  summarise(
    Relative_IL2_concentration_mean = mean(Relative_IL2_concentration),
    Relative_IL2_concentration_sd = sd(Relative_IL2_concentration),
    IL2_concentration_Dilution_Factor_mean = mean(IL2_concentration_Dilution_Factor),
    IL2_concentration_Dilution_Factor_sd = sd(IL2_concentration_Dilution_Factor)
  ) %>% 
  as.data.table()

color_doseresponse<-c("MyD88-sT6BM-1x" = "#88CCEE", 
                      "MyD88-sT6BM-3x" = "#999933",
                      "MyD88-sT6BM-5x" = "#CC6677"
)

ggplot(
  data = Dose_stats,
  aes(
    x = Stimulation_Condition,
    y = Relative_IL2_concentration_mean,
    group = Cohort,
    color = Cohort
  )
)+
  geom_path(
    size = 1.5
  )+
  color_palette(
    palette = color_doseresponse
  ) +
  geom_errorbar(
    data = Dose_stats,
    aes(
      x = Stimulation_Condition,
      y = Relative_IL2_concentration_mean,
      ymin = Relative_IL2_concentration_mean - Relative_IL2_concentration_sd,
      ymax = Relative_IL2_concentration_mean + Relative_IL2_concentration_sd
    ),
    linewidth = .75,
    width = 0.25,
    color = "black"
  )+
  geom_point(
    data = Dose_plot,
    aes(
      x = Stimulation_Condition,
      y = Relative_IL2_concentration_mean
    ),
    size = 2
  )+
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100),
    labels = scales::label_comma() #to remove scientific notation
  )+
  labs(
    y = "relative IL-2 secretion",
    x = "IL-1 (ng/mL)"
  ) +
  theme_classic(base_size = 22)+
  theme(
    legend.position = "0",
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
    )


setwd("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 3")

ggsave(
  "cl247_cl255_cl263_ELISA_doseresponse.pdf",
  plot = last_plot(),
  scale = 3,
  units = "mm",
  height = 47,
  width = 88
)
