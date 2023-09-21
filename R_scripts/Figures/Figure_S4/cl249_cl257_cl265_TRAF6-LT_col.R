# Cell Summary Table for the Dwell Time Calculation -----------------------

#Table for the Dwell Time
Cell_Summary<-
  Table %>%
  filter(
    PROTEIN == "MyD88",
    NORMALIZED_INTENSITY >= 4.5
  ) %>% 
  group_by(
    UNIVERSAL_TRACK_ID
  ) %>% 
  arrange(
    FRAME, 
    .by_group = TRUE #We order every group (ie every Track) by frames
  ) %>% 
  mutate(
    COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1.5 #threshold at which recruitment is counted
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
    CELL,
    UNIVERSAL_TRACK_ID,
    STREAK, #group continuous frames above threshold together, since they have the same STREAK value
    COLOCALIZATION
  ) %>%
  summarise(
    DWELL_FRAMES = sum(COLOCALIZATION), #number of frames that complementary protein is above threshold in a continuous stretch
    DWELL_TIME = (sum(COLOCALIZATION)-1)*4
  ) %>%
  filter(
    DWELL_FRAMES >= 1 #I want to look at transient recruitment events so even 1 frame is interesting for me, different cutoffs are prbly useful
  ) %>%
  mutate(
    CATEGORY_DWELL_TIME = 
      fcase(
        DWELL_TIME == 0, "<4 s",
        DWELL_TIME <= 8 & DWELL_TIME != 0, "4-8 s",
        DWELL_TIME >= 12 & DWELL_TIME <= 40, "12-40 s",
        DWELL_TIME > 40, ">40 s"
      )
  ) %>% 
  as.data.table()

Cell_Summary$CATEGORY_DWELL_TIME <- factor(
  Cell_Summary$CATEGORY_DWELL_TIME, levels = rev(c(
    "<4 s",
    "4-8 s",
    "12-40 s",
    ">40 s"
  )
  ))

LT_NO <- 
  Cell_Summary %>% 
  group_by(
    COHORT,
    CATEGORY_DWELL_TIME
  ) %>% 
  count(
    CATEGORY_DWELL_TIME, #Values to count
    COHORT,
    name = "N_CATEGORY_DWELL_TIME", #Name of the Column with the counts
    .drop = FALSE
  ) %>% 
  group_by(
    COHORT
  ) %>% 
  mutate(
    PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME),
    TOTAL_EVENTS = N_CATEGORY_DWELL_TIME
  ) %>% 
  as.data.table()

Mean_LT <-
  Cell_Summary %>% 
  group_by(
    COHORT
  ) %>% 
  summarise(
    LT_TRAF6 = mean(DWELL_TIME),
    SEM_LT_TRAF6 = sem(DWELL_TIME)
  )

Mean_Total <- 
  Cell_Summary %>% 
  group_by(
    COHORT,
    CATEGORY_DWELL_TIME
  ) %>% 
  count(
    CATEGORY_DWELL_TIME, #Values to count
    COHORT,
    name = "N_CATEGORY_DWELL_TIME", #Name of the Column with the counts
    .drop = FALSE
  ) %>% 
  group_by(
    COHORT
  ) %>% 
  mutate(
    PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME),
    TOTAL_EVENTS = N_CATEGORY_DWELL_TIME
  ) %>% 
  group_by(
    COHORT,
    CATEGORY_DWELL_TIME
  ) %>% 
  summarise(
    PCT_RECRUITMENT = mean(PCT_RECRUITMENT)
  ) %>% 
  as.data.table()

ggplot(
  data = Mean_Total,
  aes(
    y = COHORT,
    x = PCT_RECRUITMENT*100,
    color = COHORT,
    fill = CATEGORY_DWELL_TIME
  )
)+
  geom_col(
    width = 0.7,
    alpha = 0,
    size = 1
  )+
  geom_col(
    width = 0.7,
    size = .75,
    color = NA
  )+
  scale_fill_discrete(
    type = viridis(n=4, direction = -1)
  )+
  color_palette(
    palette = color_violin
  )+
  labs(
    x = "% of total recruitments",
  )+
  scale_y_discrete(
    labels = c(bquote(cMyD88^"1xA"), bquote(cMyD88^"3xA"), bquote(cMyD88^"5xA"))
  )+
  theme_classic(base_size = 9)+
  theme(
    legend.position = 0,
    plot.margin = margin(t = 0.25, unit = "native"),
    axis.text = element_text(color = "black",
                             size = 7),
    axis.title.y = element_blank()
  )+
  coord_cartesian(
    ylim = c(1,3),
    xlim = c(-6,100),
    clip = "off"
  )+
  annotate(
    "tile",
    y = 4.25,
    x = seq(10, 90, 20+(2/3)*10),
    width = 20,
    height = 0.7,
    fill = viridis(n=4, direction = 1)
  )+
  annotate(
    "text",
    y = 4.25,
    x = seq(10, 90, 20+(2/3)*10),
    label = c("<4", "4-8", "12-40", ">40"),
    size = 2,
    color = c("grey90", "grey90", "grey10", "grey10")
  )+
  annotate(
    "text",
    y = 4.25,
    x = -10,
    label = "TRAF6 LT (s)",
    size = 2,
    hjust = 1
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure S4")

ggsave(
  "cl249_cl257_cl265_TRAF6-LT_col.pdf",
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 50,
  width = 92
)
