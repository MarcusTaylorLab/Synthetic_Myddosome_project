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
        DWELL_TIME == 0, "0 s",
        DWELL_TIME <= 12 & DWELL_TIME != 0, "4-12 s",
        DWELL_TIME > 12, ">12 s"
      )
  ) %>% 
  as.data.table()

Cell_Summary$COHORT <- factor(
  Cell_Summary$COHORT, levels = c(
    "MyD88-GFP-synTRAF6-BD-1x TRAF6",
    "MyD88-GFP-synTRAF6-BD-3x TRAF6",
    "MyD88-GFP-synTRAF6-BD-5x TRAF6"
  )
)

Cell_Summary$CATEGORY_DWELL_TIME <- factor(
  Cell_Summary$CATEGORY_DWELL_TIME, levels = c(
    "0 s",
    "4-12 s",
    ">12 s"
  )
)

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
    PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME)
  ) %>% 
  group_by(
    COHORT,
    CATEGORY_DWELL_TIME
  ) %>% 
  summarise(
    PCT_RECRUITMENT = mean(PCT_RECRUITMENT)
  ) %>% 
  as.data.table()

color_fill <- c(
  "0 s" = "#f0f0f0",
  "4-12 s" = "#bdbdbd",
  ">12 s" = "#636363"
)

ggplot(
  data = Mean_Total,
  aes(
    x = COHORT,
    y = PCT_RECRUITMENT*100,
    color = COHORT,
    fill = CATEGORY_DWELL_TIME
  )
)+
  geom_col(
    width = 0.7,
    size = .75
  )+
  color_palette(
    palette = color_violin
  )+
  fill_palette(
    palette = color_fill
  )+
  labs(
    y = "% of total recruitments",
  )+
  theme_classic(base_size = 9)+
  theme(
    legend.position = "right",
    axis.title.x = element_blank(),
    axis.text = element_text(color = "black",
                             size = 7),
    axis.text.x = element_blank()
  ) +
  guides(
    color = "none",  # Remove legend for color aesthetic
    fill = guide_legend(title = "LT of TRAF6")  # Set legend title for fill aesthetic
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 4")

ggsave(
  "cl247_cl255_cl263_TRAF6-LT_col.pdf",
  scale = 1,
  units = "mm",
  height = 35,
  width = 60
)
