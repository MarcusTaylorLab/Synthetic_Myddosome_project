# Calculating the Recruitment time for stable TRAF& colocalization -----------------------

#Table for the Dwell Time
Recruitment_Time<-
  Table %>%
  filter(
    PROTEIN == "MyD88",
    FRAMES_ADJUSTED <= 100, #only take the frames without excessive bleaching
    FRAMES_SINCE_LANDING <= 200,
    STARTING_NORMALIZED_INTENSITY <= 6
    #filtering out noisy spots
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
    DWELL_TIME = lead(TIME_ADJUSTED, default = last(TIME_ADJUSTED)) -  TIME_ADJUSTED
  ) %>% 
  group_by(
    COHORT,
    IMAGE,
    UNIVERSAL_TRACK_ID,
    STREAK, #group continuous frames above threshold together, since they have the same STREAK value
    COLOCALIZATION
  ) %>%
  mutate(
    DWELL_FRAMES = sum(COLOCALIZATION) #number of frames that complementary protein is above threshold in a continuous stretch
  ) %>%
  filter(
    DWELL_FRAMES >= 3, #I define stable recruitment as three frames above threshold
  ) %>%
  summarise(
    RECRUITMENT_TIME = min(TIME_ADJUSTED)
  ) %>% 
  as.data.table()

Replicates_Time <- 
  Recruitment_Time %>% 
  group_by(
    COHORT,
    IMAGE
  ) %>% 
  summarise(
    RECRUITMENT_TIME = median(RECRUITMENT_TIME),
    EVENTS = n()
  ) %>% 
  as.data.table()

#plot the RECRUITMENT_TIME of different COHORTS
ggplot(
  data = Recruitment_Time,
  aes(
    y = COHORT,
    x = RECRUITMENT_TIME,
    fill = COHORT
    )
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
    data = Replicates_Time,
    position = position_jitter(height=0.2, width=0)
  )+
  labs(
    y = "Cell Lines",
    x = "Recruitment time of TRAF6 (s)"
  # )+
  # scale_x_continuous(
  #   trans = "log10",
  #   breaks = scales::breaks_log(n = 8, base = 10)
  )+
  fill_palette(
    palette = color_violin
  )+
  theme_classic(base_size = 22)+
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color = "black")
  )