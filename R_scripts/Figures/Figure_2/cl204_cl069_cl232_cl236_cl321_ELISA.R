All_days_data <-
  All_days_data %>% 
  group_by(
    Sample_Day,
    Date
  ) %>% 
  mutate(
    Cohort = fcase(
      Cohort == "3E10_GFP", "MyD88-GFP",
      Cohort != "3E10_GFP", Cohort)
  ) %>% 
  group_by(
    Sample_Day,
    Date
  ) %>% 
  mutate(REL_IL2 = Values_Measured / Values_Measured[Cohort == "MyD88-GFP" & Stimulation_Condition == "Stimulated"])

Factor <- 
  All_days_data %>% 
  ungroup() %>% 
  filter(
    Cohort == "MyD88-GFP/ TRAF6-mScarlet",
    Stimulation_Condition == "Stimulated"
  ) %>% 
  summarise(
    REL_IL2 = mean(REL_IL2)
  ) %>% 
  as.numeric()

All_days_data <-
  All_days_data %>% 
  mutate(
    Relative_Intensity = REL_IL2*(1/Factor)
  )

# #Figure 2: DHF91-TIR-T6BM, WT, 3x KO-------------------------------------------------------------------------
DHF91_TIR_T6BM_data <- 
  All_days_data %>% 
  filter(
    Cohort %in% c("DHF91-TIR-T6BM low", "MyD88-T6BM"),
    Date == "20220623"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-T6BM" ~ "chMyD88",
      Cohort == "DHF91-TIR-T6BM low" ~ "chMyD88-DHF91"))

Triple_KO <-
  All_days_data %>% 
  filter(
    Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-" ~ "3xKO"))

bDLD_data <-
  All_days_data %>% 
  filter(
    Cohort == "BDLD_57_H"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "BDLD_57_H" ~ "chMyD88-bDLD"))

DHF91_TIR_T6BM_data <- rbind(DHF91_TIR_T6BM_data, Triple_KO, bDLD_data)


DHF91_TIR_T6BM_data$Cohort <-
  factor(DHF91_TIR_T6BM_data$Cohort,
         levels=c( "3xKO",
                   "chMyD88",
                   "chMyD88-bDLD",
                   "chMyD88-DHF91"
         ))

DHF91_TIR_T6BM_test <-
  DHF91_TIR_T6BM_data %>% 
  mutate(
    IL2_concentration_Dilution_Factor_mean = IL2_concentration_Dilution_Factor,
    Relative_Intensity_mean = Relative_Intensity
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
    Relative_Intensity_mean = mean(Relative_Intensity)
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
    IL2_concentration_Dilution_Factor_sd = sd(IL2_concentration_Dilution_Factor),
    Relative_Intensity_mean = mean(Relative_Intensity),
    Relative_Intensity_sem = sem(Relative_Intensity)
  ) %>% 
  as.data.table()

color_elisa <- c("Unstimulated" = "white",
                 "Stimulated" = "grey")

ggplot(
  data = DHF91_TIR_T6BM_stats,
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
    data = DHF91_TIR_T6BM_stats,
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
    data = DHF91_TIR_T6BM_plot,
    aes(
      x = Cohort,
      y = Relative_Intensity_mean
    ),
    size = 0.75,
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4)
  )+
  stat_compare_means(
    data = DHF91_TIR_T6BM_test,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
    )+
  scale_y_continuous(
    breaks = seq(from = 0, to = 1, by = 0.2)
  )+
  scale_x_discrete(
    labels = c("3xKO", "cMyD88", bquote(cMyD88^bDLD), bquote(cMyD88^DHF91))
    
  )+
  scale_fill_manual(
    labels = c("- IL-1", "+ IL-1"),
    values = c("white", "grey50")
  )+
  labs(
    y = "relative IL-2 secretion"
  )+
  guides(
    color = "none"
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
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 9, color = "black"),
    legend.key.size = unit(2, "mm"),
    #legend.key = element_rect(fill = "transparent")  # Set the legend box background to transparent
    )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 2")

ggsave(
  "cl204_cl069_cl232_cl236_cl321_ELISA.pdf",
  plot = last_plot(),
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 35,
  width = 86
)
