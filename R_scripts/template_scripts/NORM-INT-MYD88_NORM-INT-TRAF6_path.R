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
    NORMALIZED_INTENSITY = 2*round(NORMALIZED_INTENSITY/2) #so we round to even integers
  ) %>% 
  group_by(
    COHORT,
    NORMALIZED_INTENSITY
  ) %>% 
  filter(
    n() >= 25 #so we only look at intensities where there are at least 5 events
  ) %>% 
  summarise(
    MEAN_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = mean(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>%
  group_by(
    COHORT
  ) %>% 
  # mutate(
  #   SMOOTH_NORMALIZED_RECRUITMENT = signal::sgolayfilt(NORMALIZED_RECRUITMENT, p = 1, n = 25) #apply a smoothing filter to add to the plot
  # ) %>% 
  as.data.table()

#plot the percentage of TRAF6 recruitment over normalized Intensity
ggplot(
  data = Recruitment_List
)+
  geom_path(
    aes(
      x = NORMALIZED_INTENSITY,
      y = MEAN_COMPLEMENTARY_NORMALIZED_INTENSITY_1,
      group = COHORT,
      color = COHORT
    ),
    alpha = 0.5,
    linewidth = 1.5
  )+
  # geom_path(
  #   aes(
  #     x = NORMALIZED_INTENSITY,
  #     y = SMOOTH_NORMALIZED_RECRUITMENT*100,
  #     group = COHORT,
  #     color = COHORT
  #   ),
  #   linewidth = 1
  # )+
  color_palette(
    palette = color_violin
  )+
  scale_x_continuous(
    limits = c(0, 100),
    breaks = scales::breaks_width(40)
  )+
  # scale_y_continuous(
  #   trans = "log1p",
  #   limits = c(0, 100),
  #   breaks = c(0,5,10,25,50,75,100)
  # )+ 
  labs(
    x = "Normalized Intensity of chimeric MyD88",
    y = "Normalized Intensity of TRAF6"
  )+
  theme_classic(base_size = 22)+
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text = element_text(color = "black")
  )+
  guides(
    color = guide_legend(nrow = 4)
  )
