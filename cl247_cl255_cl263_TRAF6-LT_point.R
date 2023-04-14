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
    #DWELL_TIME = shift(TIME_ADJUSTED, n = 1, type = "lead", fill = last(TIME_ADJUSTED)) -  TIME_ADJUSTED
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

# Cum Freq --------------------------------------------------------

#create labels for the x-axis
c <-as.character(seq(from = 0, to = 104, by = 16))
xlabels <- paste("<", c, sep = "")
xlabels <- c(xlabels)

ggplot(
  data = Cell_Summary,
  aes(
    x = DWELL_TIME,
    #y = 1-..y.., #unhash to get the inverse cumulative frequency
    group = COHORT,
    color = COHORT
  )
)+
  stat_ecdf(
    geom = "point",
    size = 2,
    pad = F #so no data points get added to the sides
  )+
  coord_cartesian(
    xlim = c(0,100),
    ylim = c(0,1)
  )+
  scale_x_continuous(
    #limits = c(0,100),
    breaks = seq(from = 0, to = 104, by = 16),
    labels = xlabels
  )+
  labs(
    x = "Lifetime of TRAF6 (s)",
    y = "Cumulative Density"
  )+
  scale_color_discrete(
    name = "Cell Line",
    type = color_violin,
    labels = c("1x", "3x", "5x")
  )+
  theme_bw(base_size = 22)+
  theme(
    legend.position = "0",
    strip.text = element_blank(),
    axis.text = element_text(color = "black")
  )

setwd("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 3")

ggsave(
  "cl247_cl255_cl263_TRAF6-LT_point.pdf",
  scale = 3,
  units = "mm",
  height = 40,
  width = 60
)
