# Cell Summary Table for the Dwell Time Calculation -----------------------

#Table for the Dwell Time
Cell_Summary<-
  Table %>%
  filter(
    PROTEIN == "MyD88",
    NORMALIZED_INTENSITY >= 1
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  arrange(
    FRAME, 
    .by_group = TRUE #We order every group (ie every Track) by frames
  ) %>% 
  mutate(
    COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 2 #threshold at which recruitment is counted
    #Here I create a new column (with mutate) that will have the value 1 (for TRUE), 
    #if the condition above is satisfied or 0 (for FALSE) if the condition is not satisfied
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>%
  mutate(
    STREAK = cumsum(!COLOCALIZATION), 
    #Here i create a new column were I calculate the cumulative sum for all the rows in the column COLOCALIZATION,
    #But since TRUE==1 and FALSE==0 we have to add the ! before COLOCALIZATION
    #This results in all the rows following each other with a continuous "TRUE" get the same "STREAK value"
  ) %>% 
  group_by(
    COHORT,
    IMAGE,
    UNIVERSAL_TRACK_ID,
    STREAK, #group continuous frames above threshold together, since they have the same STREAK value
    COLOCALIZATION
  ) %>%
  summarise(
    DWELL_FRAMES = sum(COLOCALIZATION), #number of frames that complementary protein is above threshold in a continuous stretch
    DWELL_TIME = (DWELL_FRAMES-1)*4
  ) %>%
  filter(
    DWELL_FRAMES >= 1, #I want to look at transient recruitment events so even 1 frame is interesting for me, different cutoffs are prbly useful
  ) %>%
  as.data.table()

Cell_Summary$COHORT <- factor(
  Cell_Summary$COHORT, levels = c("MyD88-GFP-synTRAF6-BD-1x TRAF6",
                                  "MyD88-GFP-synTRAF6-BD-3x TRAF6",
                                  "MyD88-GFP-synTRAF6-BD-5x TRAF6"
  )
)

# Histogram --------------------------------------------------------
ggplot(
  data = Cell_Summary,
  aes(
    x = DWELL_TIME,
    y = ..density..*100,
    group = COHORT,
    fill = COHORT
  )
)+
  geom_histogram(
    binwidth = 4,
    #color ="black",
    center = 0,
    position = "dodge"
  )+
  scale_y_cut(
    breaks = 2.3,
    which = c(1,2),
    scales = c(1,3),
  )+
  scale_x_continuous(
    limits = c(NA,150)
  )+
  labs(
    x = "Lifetime of TRAF6 (s)",
    y = "% of recruitments"
  )+
  fill_palette(
    palette = color_violin
  )+
  theme_classic(base_size = 22)+
  theme(
    legend.position = "0",
    strip.text = element_blank(),
    axis.text = element_text(color = "black")
  )
