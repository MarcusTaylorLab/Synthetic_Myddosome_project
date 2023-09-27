# Recruitment of TRAF6 over Normalized Intensity --------------------------
Recruitment_List <- 
  Table %>%
  filter(
    NORMALIZED_INTENSITY >= 1, #filter out the noise
    FRAMES_ADJUSTED <=100, #only take the frames without excessive bleaching
    FRAMES_SINCE_LANDING <= 200,
    PROTEIN == "MyD88"
  ) %>% 
  mutate(
    RECRUITMENT = fcase(
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1.5, "1",
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 < 1.5, "0" #If TRAF6 is colocalized RECRUITMENT will be 1, otherwise 0
    ),
    NORMALIZED_INTENSITY = 2*round(NORMALIZED_INTENSITY/2) #so we round to even integers
  ) %>% 
  group_by(
    COHORT,
    NORMALIZED_INTENSITY,
    IMAGE
  ) %>% 
  filter(
    n() >= 10 #so we only look at intensities where there are at least 10 events
  ) %>% 
  summarise(
    NORMALIZED_RECRUITMENT = (sum(RECRUITMENT == 1)/n()),
    SEM_NORMALIZED_RECRUITMENT = sem(RECRUITMENT),
    SD_NORMALIZED_RECRUITMENT = sd(RECRUITMENT)
  ) %>% 
  as.data.table()

Mean_of_Means <-
  Recruitment_List %>% 
  group_by(
    COHORT,
    NORMALIZED_INTENSITY
  ) %>% 
  filter(
    n() >= 2 #so we only look at intensities where there are at least 2 replicates
  ) %>%
  summarise(
    SEM_NORMALIZED_RECRUITMENT = sem(NORMALIZED_RECRUITMENT),
    SD_NORMALIZED_RECRUITMENT = sd(NORMALIZED_RECRUITMENT),
    NORMALIZED_RECRUITMENT = mean(NORMALIZED_RECRUITMENT)
  )

#plot the percentage of TRAF6 recruitment over normalized Intensity
ggplot(
  data = Mean_of_Means
)+
  geom_path(
    aes(
      x = NORMALIZED_INTENSITY,
      y = NORMALIZED_RECRUITMENT*100,
      group = COHORT,
      color = COHORT
    ),
    alpha = 1,
    linewidth = 0.5
  )+
  geom_ribbon(
    aes(
      x = NORMALIZED_INTENSITY,
      ymin = (NORMALIZED_RECRUITMENT-SEM_NORMALIZED_RECRUITMENT)*100,
      ymax = (NORMALIZED_RECRUITMENT+SEM_NORMALIZED_RECRUITMENT)*100,
      group = COHORT,
      fill = COHORT
    ),
    alpha = 0.3,
    linewidth = 0
  )+
  color_palette(
    palette = color_violin
  )+
  fill_palette(
    palette = color_violin
  )+
  scale_x_continuous(
    limits = c(0, 100),
    breaks = scales::breaks_width(20)
  )+
  labs(
    x = "Size of chimeric oligomer",
    y = "% of TRAF6 \n pos. oligomers"
  )+
  theme_classic(base_size = 9)+
  theme(
    legend.position = "0",
    axis.text = element_text(color = "black",
                             size = 7),
    legend.title = element_blank()
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 4")

ggsave(
  "cl247_cl255_cl263_NORM-INT_PCT-TRAF6-REC_path.pdf",
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 30,
  width = 90
)
