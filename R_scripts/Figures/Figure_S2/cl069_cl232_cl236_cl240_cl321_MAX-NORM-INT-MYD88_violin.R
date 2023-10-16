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
    IMAGE,
    COHORT
  ) %>% 
  summarize(
    LIFETIME = median(LIFETIME),
    MAX_NORMALIZED_INTENSITY = mean(MAX_NORMALIZED_INTENSITY),
    MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1 = median(MAX_COMPLEMENTARY_NORMALIZED_INTENSITY_1)
  )

#Plot max intensities of MyD88
ggplot(
  data = Tracks_MyD88,
    aes(
      y = COHORT,
      x = MAX_NORMALIZED_INTENSITY,
      fill = COHORT,
      group = COHORT
      )
)+
  geom_violin(
    alpha = 0.7,
    scale = "width",
    linewidth = 0.4,
    color = "black"
    
  )+
  geom_boxplot(
    width = .2,
    outlier.shape = NA,
    fill = NA,
    linewidth = 0.4,
    color = "black"
  )+
  geom_point(
    data = Replicates,
    position = position_jitter(height=0.3, width=0),
    size = 0.75
  )+
  labs(
    y = "Cell lines",
    x = "Max size of chimeric oligomer (log-scale)"
  )+
  coord_cartesian(
    xlim = c(1,NA)
  )+
  scale_x_continuous(
    trans = "log10",
    breaks = scales::breaks_log(n=10,base=10)
  )+
  scale_y_discrete(
    labels = rev(c("WT", "cMyD88", bquote(cMyD88^TIR), bquote(cMyD88^bDLD), bquote(cMyD88^DHF91)))
  )+
  fill_palette(
    palette = color_violin
  )+
  theme_classic(base_size = 7)+
  theme(
    legend.position = "0",
    legend.title = element_blank(),
    axis.text.x = element_text(color = "black",
                               size = 6),
    axis.title.y = element_blank()
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure S2")

ggsave(
  "cl069_cl232_cl236_cl240_cl321__MAX-NORM-INT-MYD88_violin.pdf",
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 40,
  width = 90
)
