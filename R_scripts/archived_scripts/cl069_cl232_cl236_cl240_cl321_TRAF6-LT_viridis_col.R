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
    COLOCALIZATION = COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1 #threshold at which recruitment is counted
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

Cell_Summary$CATEGORY_DWELL_TIME <- factor(
  Cell_Summary$CATEGORY_DWELL_TIME, levels = c(
    "0 s",
    "4-12 s",
    ">12 s"
  )
)

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
    DWELL_TIME
  ) %>% 
  count(
    DWELL_TIME, #Values to count
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
    DWELL_TIME
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
    fill = DWELL_TIME
  )
)+
  geom_col(
    width = 0.7
  )+
  binned_scale(aesthetics = "fill",
               scale_name = "stepsn", 
               palette = function(x) viridis(5),
               breaks = c(0, 4, 8, 12, 40),
               right = FALSE,
               guide = "colorbar",
               labels = c(0, 4, 8, "12-40", "≥40")
  )+
  labs(
    x = "% of total recruitments",
  )+
  scale_y_discrete(
    labels = c(bquote(chMyD88^DHF91), bquote(chMyD88^bDLD), bquote(chMyD88^TIR), "chMyD88", "WT")
  )+
  theme_classic(base_size = 7)+
  guides(
    fill = guide_legend(title = "TRAF6 LT (s)",
                        title.position = "top",
                        title.hjust = 0
                        )  # Set legend title for fill aesthetic
  )+
  theme(
    legend.position = "top",
    legend.box.spacing = margin(0.5),
    legend.justification = "right",
    axis.text = element_text(color = "black",
                             size = 6),
    legend.text = element_text(color = "black",
                               size = 6,
                               hjust = -1
                              ),
    axis.title.y = element_blank(),
    legend.key.size = unit(2, "mm"),
    legend.title = element_text(size = 8)
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 2")

ggsave(
  "cl069_cl232_cl236_cl240_cl321_TRAF6-LT_col.pdf",
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 45,
  width = 45
)
