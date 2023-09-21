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
    NORMALIZED_INTENSITY = round(NORMALIZED_INTENSITY) #so we round to even integers
  ) %>% 
  group_by(
    COHORT,
    NORMALIZED_INTENSITY
  ) %>% 
  filter(
    n() >= 20 #so we only look at intensities where there are at least 5 events
  ) %>% 
  summarise(
    count = n(),
    NORMALIZED_RECRUITMENT = (sum(RECRUITMENT == 1)/n()),
    SEM_NORMALIZED_RECRUITMENT = sem(RECRUITMENT),
    SD_NORMALIZED_RECRUITMENT = sd(RECRUITMENT)
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
    alpha = 1,
    linewidth = 0.5
  )+
  color_palette(
    palette = color_violin
  )+
  fill_palette(
    palette = color_violin
  )+
  scale_x_continuous(
    limits = c(0, 180),
    breaks = scales::breaks_width(20)
  )+
  labs(
    x = "Size of chimeric oligomer",
    y = "oligomers \n colocalizing with TRAF6"
  )+
  theme_classic(base_size = 9)+
  theme(
    legend.position = "0",
    axis.text = element_text(color = "black",
                             size = 7),
    legend.title = element_blank()
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure S4")

ggsave(
  "cl249_cl257_cl265_NORM-INT_PCT-TRAF6-REC_path.pdf",
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 30,
  width = 90
)
