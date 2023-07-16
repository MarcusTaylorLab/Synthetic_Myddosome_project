# Plot the maximal normalized intensity of all tracks as a proxy for puncta size --------

#start grouping and summarizing the data
Tracks_MyD88<-
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
    FRAMES_ADJUSTED == min(FRAMES_ADJUSTED)
  ) %>% 
  as.data.table()

#Add stats for the images taken to clean up violin
Replicates<-
  Tracks_MyD88 %>% 
  group_by(
    PROTEIN,
    IMAGE,
    COHORT
  ) %>% 
  summarize(
    LIFETIME = median(LIFETIME),
    MAX_NORMALIZED_INTENSITY = median(MAX_NORMALIZED_INTENSITY),
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = median(MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  ) %>% 
  as.data.table()

#Plot max intensities of MyD88
ggplot(
  data = Tracks_MyD88,
    aes(
      y = COHORT,
      x = MAX_NORMALIZED_INTENSITY,
      fill = COHORT)
)+
  geom_violin(
    alpha = 0.7,
    scale = "width"
  )+
  geom_boxplot(
    width = .1,
    outlier.shape = NA,
    fill = NA
  )+
  geom_point(
    data = Replicates,
    position = position_jitter(height=0.2, width=0)
  )+
  labs(
    y = "Cell Lines",
    x = "Max Normalized Intensity of Chimeric MyD88s"
  )+
  scale_x_continuous(
    trans = "log10",
    breaks = scales::breaks_log(n = 8, base = 10)
  )+
  fill_palette(
    palette = color_violin
  )+
  theme_classic(base_size = 22)+
  theme(
    legend.position = "0",
    legend.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black")
  )

setwd("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 3")

ggsave(
  "cl247_cl255_cl263_MAX-NORM-INT-MYD88_violin.pdf",
  scale = 3,
  units = "mm",
  height = 40,
  width = 60
)