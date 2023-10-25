WT_data <-
  All_days_data %>% 
  filter(
    Cohort == "WT",
    Date == "20221124"
  ) %>% 
  mutate(
    Stimulation_Condition = as.numeric(Stimulation_Condition),
    Stimulation_Condition = case_when(
      Stimulation_Condition == 0 & Cohort == "WT" ~ "Unstimulated",
      Stimulation_Condition == 10 & Cohort == "WT" ~ "Stimulated"
    )
  ) %>% 
  filter(
    !is.na(Stimulation_Condition )
  )

syn_data <- 
  All_days_data %>% 
  filter(
    Cohort %in% c("MyD88-GFP-synTRAF6-BD-1xA", "MyD88-GFP-synTRAF6-BD-3xA", "MyD88-GFP-synTRAF6-BD-5xA"),
    Date == "20221124"
  ) %>% 
  group_by(
    Sample_Day
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-GFP-synTRAF6-BD-1xA" ~ "1xA",
      Cohort == "MyD88-GFP-synTRAF6-BD-3xA" ~ "3xA",
      Cohort == "MyD88-GFP-synTRAF6-BD-5xA" ~ "5xA"
    ))

syn_T6BM_data <- rbind(syn_data, WT_data) %>% 
  group_by(
    Sample_Day
  ) %>% 
  mutate(
    Relative_Intensity = (Values_Measured-min(Values_Measured))/(max(Values_Measured)-min(Values_Measured))
  )

syn_T6BM_data$Cohort <-
  factor(syn_T6BM_data$Cohort,
         levels=c( "WT",
                   "1xA",
                   "3xA",
                   "5xA"
         ))

syn_T6BM_data$Stimulation_Condition <-
  factor(syn_T6BM_data$Stimulation_Condition,
         levels=c("Unstimulated",
                  "Stimulated"
         ))

syn_T6BM_test <-
  syn_T6BM_data %>% 
  mutate(
    IL2_concentration_Dilution_Factor_mean = IL2_concentration_Dilution_Factor,
    Relative_Intensity_mean = Relative_Intensity
  ) %>% 
  as.data.table()

syn_T6BM_plot <- 
  syn_T6BM_data %>% 
  group_by(
    Cohort,
    Stimulation_Condition,
    Sample_Day
  ) %>% 
  summarise(
    IL2_concentration_Dilution_Factor_mean = mean(IL2_concentration_Dilution_Factor),
    Relative_Intensity_mean = mean(Relative_Intensity)
  ) %>% 
  as.data.table()

syn_T6BM_stats <-
  syn_T6BM_data %>% 
  group_by(
    Cohort,
    Stimulation_Condition
  ) %>% 
  summarise(
    IL2_concentration_Dilution_Factor_mean = mean(IL2_concentration_Dilution_Factor),
    IL2_concentration_Dilution_Factor_median = median(IL2_concentration_Dilution_Factor),
    IL2_concentration_Dilution_Factor_sd = sd(IL2_concentration_Dilution_Factor),
    Relative_Intensity_mean = mean(Relative_Intensity),
    Relative_Intensity_sem = sem(Relative_Intensity)
  ) %>% 
  as.data.table()

color_elisa <- c("Unstimulated" = "white",
                 "Stimulated" = "grey")

ggplot(
  data = syn_T6BM_stats,
  aes(
    x = Cohort,
    y = Relative_Intensity_mean,
    fill = Stimulation_Condition
  )
) +
  geom_col(
    position = position_dodge(width = 0.5),
    color = "black",
    width = 0.5
  ) +
  geom_errorbar(
    data = syn_T6BM_stats,
    aes(
      x = Cohort,
      ymin = Relative_Intensity_mean - Relative_Intensity_sem,
      ymax = Relative_Intensity_mean + Relative_Intensity_sem
    ),
    linewidth = .75,
    position = position_dodge(width = 0.5),
    width = 0.25
  ) +
  geom_point(
    data = syn_T6BM_plot,
    aes(
      x = Cohort,
      y = Relative_Intensity_mean
    ),
    size = 0.75,
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4)
  )+
  stat_compare_means(
    data = syn_T6BM_test,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  )+
  scale_y_continuous(
    breaks = seq(from = 0, to = 1, by = 0.2)
  )+
  scale_fill_manual(
    labels = c("- IL-1", "+ IL-1"),
    values = c("white", "grey50")
  )+
  labs(
    y = "IL-2 release (relative)"
  )+
  guides(
    color = "none",
    fill = guide_legend(nrow = 2)
  )+
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 8, 
                               colour = "black",
                               angle = 0,
                               vjust = 0.6
    ),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 9),
    legend.position = c(0.9,0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 9, color = "black"),
    legend.key.size = unit(3, "mm"),
    #legend.key = element_rect(fill = "transparent")  # Set the legend box background to transparent
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure S4")

ggsave(
  "WT_cl249_cl257_cl265_ELISA.pdf",
  plot = last_plot(),
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 35,
  width = 92
)
