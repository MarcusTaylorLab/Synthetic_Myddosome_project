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
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 2, "1",
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 <2, "0" #If TRAF6 is colocalized RECRUITMENT will be 1, otherwise 0
    ),
    NORMALIZED_INTENSITY = round(NORMALIZED_INTENSITY) #so we round to even integers
  ) %>% 
  group_by(
    COHORT,
    NORMALIZED_INTENSITY
  ) %>% 
  filter(
    n() >= 25 #so we only look at intensities where there are at least 5 events
  ) %>% 
  summarise(
    NORMALIZED_RECRUITMENT = (sum(RECRUITMENT == 1)/n())
  # ) %>%
  # group_by(
  #   COHORT
  # ) %>%
  # mutate(
  #   SMOOTH_NORMALIZED_RECRUITMENT = signal::sgolayfilt(NORMALIZED_RECRUITMENT, p = 3, n = 25) #apply a smoothing filter to add to the plot
  ) %>% 
  as.data.table()

#plot the percentage of TRAF6 recruitment over normalized Intensity
ggplot(
  data = Recruitment_List
)+
  geom_path(
    aes(
      x = NORMALIZED_INTENSITY,
      y = NORMALIZED_RECRUITMENT*100,
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
    breaks = scales::breaks_width(20)
  )+
  # scale_y_continuous(
  #   trans = "log1p",
  #   limits = c(0, 100),
  #   breaks = c(0,5,10,25,50,75,100)
  # )+ 
  labs(
    x = "Normalized Intensity of chimeric MyD88",
    y = "% colocalization with TRAF6"
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

setwd("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 2")

ggsave(
  "cl069_cl232_cl240_NORM-INT_PCT-TRAF6-REC_path.pdf",
  scale = 3,
  units = "mm",
  height = 60,
  width = 50
)
