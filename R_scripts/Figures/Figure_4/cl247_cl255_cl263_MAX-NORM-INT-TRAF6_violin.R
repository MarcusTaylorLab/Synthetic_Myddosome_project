# Comparing Max Normalized Intensities of TRAF6 ------------------------------------

#start grouping and summarizing the data
Tracks_TRAF6<-
  Table %>% 
  filter(
    NORMALIZED_INTENSITY >= 1, #filter out the noise
    FRAMES_ADJUSTED <=100, #only take the frames without excessive bleaching
    FRAMES_SINCE_LANDING <= 200,
    PROTEIN == "MyD88"
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  mutate(
    LIFETIME = max(TIME_ADJUSTED),
    LIFETIME_FRAMES = max(FRAMES_ADJUSTED),
    MAX_NORMALIZED_INTENSITY = max(NORMALIZED_INTENSITY),
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = max(COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  filter(
    LIFETIME_FRAMES >= 3, #only tracks with at least three spots
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1.5,
    FRAMES_ADJUSTED == min(FRAMES_ADJUSTED)
  ) %>% 
  as.data.table()

#Add stats for the images taken to clean up violin
Replicates<-
  Tracks_TRAF6 %>% 
  group_by(
    PROTEIN,
    IMAGE,
    COHORT
  ) %>% 
  summarize(
    LIFETIME = mean(LIFETIME),
    MAX_NORMALIZED_INTENSITY = mean(MAX_NORMALIZED_INTENSITY),
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = mean(MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  as.data.table()

#Plot max intensities of TRFA6
ggplot(
  data = Tracks_TRAF6,
  aes(
   x = COHORT,
   y = MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1,
   fill = COHORT
       )
)+
  geom_violin(
    alpha = 0.7,
  )+
  geom_boxplot(
    width = .2,
    outlier.shape = NA,
    fill = NA
  )+
  fill_palette(
    palette = color_violin
  )+
  geom_point(
    data = Replicates,
    position = position_jitter(height=0.2, width=0)
  )+
  scale_y_log10(
    breaks = scales::breaks_log(n = 10, base = 10)
  )+
  scale_x_discrete(
    labels = c("1x", "3x", "5x")
  )+
  labs(
    y = "max size of TRAF6 (log)"
  )+
  coord_cartesian(
    ylim = c(1.5, NA)
  )+
  theme_classic(base_size = 9)+
  theme(
    legend.position = "0",
    legend.title = element_blank(),
    axis.title.y = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",
                             size = 7)
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 4")

ggsave(
  "cl247_cl255_cl263_MAX-NORM-INT-TRAF6_violin.pdf",
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 50,
  width = 46
)
