# All_days_data <-
#   All_days_data %>% 
#   group_by(
#     Sample_Day,
#     Date
#   ) %>% 
#   mutate(REL_IL2 = 
#            case_when(
#              Date == "20220609" ~ Values_Measured / Values_Measured[Cohort == "WT" & Stimulation_Condition == "Stimulated"],
#              Cohort != "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-" ~ Values_Measured / Values_Measured[Cohort == "MyD88-GFP/ TRAF6-mScarlet" & Stimulation_Condition == "Stimulated"]))
# 
# Factor <- 
#   All_days_data %>% 
#   ungroup() %>% 
#   filter(
#     Cohort == "MyD88-GFP/ TRAF6-mScarlet",
#     Stimulation_Condition == "Stimulated"
#   ) %>% 
#   summarise(
#     REL_IL2 = mean(REL_IL2)
#   ) %>% 
#   as.numeric()
# 
# All_days_data <-
#   All_days_data %>% 
#   mutate(
#     Relative_Intensity = REL_IL2*(1/Factor)
#   )

chimeric_MyD88 <- 
  All_days_data %>% 
  filter(
    Date %in% c("20220623", "20220701")
  ) %>% 
  group_by(
    Sample_Day,
    Date
  ) %>% 
  mutate(
    Relative_Intensity = Values_Measured / Values_Measured[Cohort == "MyD88-GFP/ TRAF6-mScarlet" & Stimulation_Condition == "Stimulated"])

# #Figure 1:-------------------------------------------------------------------------
chimeric_MyD88_data <- 
  chimeric_MyD88 %>% 
  filter(
    Cohort %in% c("MyD88-T6BM", "MyD88-TIR-T6BM low", "MyD88-GFP/ TRAF6-mScarlet"),
    Date == "20220623"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-GFP/ TRAF6-mScarlet" ~ "WT",
      Cohort == "MyD88-T6BM" ~ "cMyD88",
      Cohort == "MyD88-TIR-T6BM low" ~ "cMyD88-TIR")
    )

Triple_KO <-
  All_days_data %>% 
  filter(
    Date == "20220609",
  ) %>% 
  group_by(
    Sample_Day,
    Date
  ) %>% 
  mutate(
    Relative_Intensity = Values_Measured / Values_Measured[Cohort == "WT" & Stimulation_Condition == "Stimulated"]
    ) %>% 
  filter(
    Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-/-/ IRAK4-/-/\\n IRAK1-/-" ~ "3xKO"))

chimeric_MyD88_3xA_data <- 
  chimeric_MyD88 %>% 
  filter(
    Cohort %in% c("MyD88-T6BM-3xA"),
    Date == "20220701"
  ) %>% 
  mutate(
    Cohort = case_when(
      Cohort == "MyD88-T6BM-3xA" ~ "cMyD88-3xA"))

chimeric_MyD88_data <- rbind(chimeric_MyD88_data, Triple_KO, chimeric_MyD88_3xA_data)

chimeric_MyD88_data$Cohort <-
  factor(chimeric_MyD88_data$Cohort,
         levels=c( "WT", 
                   "3xKO", 
                   "cMyD88",
                   "cMyD88-TIR",
                   "cMyD88-3xA"
         ))

chimeric_MyD88_test <-
  chimeric_MyD88_data %>% 
  mutate(
    IL2_concentration_Dilution_Factor_mean = IL2_concentration_Dilution_Factor,
    Relative_Intensity_mean = Relative_Intensity
  ) %>% 
  as.data.table()

chimeric_MyD88_plot <- 
  chimeric_MyD88_data %>% 
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

chimeric_MyD88_stats <-
  chimeric_MyD88_data %>% 
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

# chimeric_MyD88_test <-
#   chimeric_MyD88_data %>% 
#   mutate(
#     IL2_concentration_Dilution_Factor_mean = IL2_concentration_Dilution_Factor
#   ) %>% 
#   as.data.table()
# 
# chimeric_MyD88_plot <- 
#   chimeric_MyD88_data %>% 
#   group_by(
#     Cohort,
#     Stimulation_Condition,
#     Sample_Day
#   ) %>% 
#   summarise(
#     IL2_concentration_Dilution_Factor_mean = mean(IL2_concentration_Dilution_Factor),
#     IL2_concentration_Dilution_Factor_median = median(IL2_concentration_Dilution_Factor)
#   ) %>% 
#   as.data.table()
# 
# chimeric_MyD88_stats <-
#   chimeric_MyD88_data %>% 
#   group_by(
#     Cohort,
#     Stimulation_Condition
#   ) %>% 
#   summarise(
#     IL2_concentration_Dilution_Factor_mean = mean(IL2_concentration_Dilution_Factor),
#     IL2_concentration_Dilution_Factor_median = median(IL2_concentration_Dilution_Factor),
#     IL2_concentration_Dilution_Factor_sd = sd(IL2_concentration_Dilution_Factor)
#   ) %>% 
#   as.data.table()
# 
# Plot_Fx(chimeric_MyD88_plot, chimeric_MyD88_stats, chimeric_MyD88_test)

color_elisa <- c("Unstimulated" = "white",
                 "Stimulated" = "grey")

ggplot(
  data = chimeric_MyD88_stats,
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
    data = chimeric_MyD88_stats,
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
    data = chimeric_MyD88_plot,
    aes(
      x = Cohort,
      y = Relative_Intensity_mean
    ),
    size = 0.75,
    position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.4)
  )+
  stat_compare_means(
    data = chimeric_MyD88_test,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  )+
  scale_y_continuous(
    breaks = seq(from = 0, to = 1, by = 0.2)
  )+
  scale_x_discrete(
    labels = c("WT", "3xKO", "cMyD88", bquote(cMyD88^TIR), bquote(cMyD88^"3xA"))
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

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 1")

ggsave(
  "cl069_cl232_cl234_cl240_ELISA.pdf",
  plot = last_plot(),
  scale = 1,
  units = "mm",
  height = 50,
  width = 80
)
