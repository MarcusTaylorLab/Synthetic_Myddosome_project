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
        DWELL_TIME == 0, "<4 s",
        DWELL_TIME <= 8 & DWELL_TIME != 0, "4-8 s",
        DWELL_TIME >= 12 & DWELL_TIME <= 40, "12-40 s",
        DWELL_TIME > 40, ">40 s"
      )
  ) %>% 
  as.data.table()

Cell_Summary$CATEGORY_DWELL_TIME <- factor(
  Cell_Summary$CATEGORY_DWELL_TIME, levels = c(
    "<4 s",
    "4-8 s",
    "12-40 s",
    ">40 s"
  )
)

Cell_Summary$COHORT <- factor(
  Cell_Summary$COHORT, levels = c(
    "MyD88 TRAF6",
    "MyD88-TRAF6-BD TRAF6",
    "MyD88-TIR-TRAF6-BD TRAF6",
    "BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6",
    "MyD88-DHF91-TRAF6-BD TRAF6"
  )
)

Mean_Cell <-
  Cell_Summary %>% 
  group_by(
    COHORT,
    IMAGE,
    CELL,
    CATEGORY_DWELL_TIME
  ) %>% 
  count(
    CATEGORY_DWELL_TIME, #Values to count
    COHORT,
    IMAGE,
    CELL,
    name = "N_CATEGORY_DWELL_TIME", #Name of the Column with the counts
    .drop = FALSE
  ) %>% 
  group_by(
    COHORT,
    IMAGE,
    CELL
  ) %>% 
  filter(
    sum(N_CATEGORY_DWELL_TIME) > 4 #more than 4 recruitment events per cell
  ) %>% 
  mutate(
    PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME),
    DATE = strsplit(IMAGE, " ")[[1]][1],
    SHAPE = fcase(
      DATE == "20230405", 21,
      DATE == "20230413", 22,
      DATE == "20220516", 23,
      DATE == "20220610", 24,
      DATE == "20220615", 25
    )
  ) %>%
  as.data.table()

Mean_Replicates <-
  Cell_Summary %>% 
  group_by(
    COHORT,
    IMAGE,
    CATEGORY_DWELL_TIME
  ) %>% 
  count(
    CATEGORY_DWELL_TIME, #Values to count
    COHORT,
    IMAGE,
    name = "N_CATEGORY_DWELL_TIME", #Name of the Column with the counts
    .drop = FALSE
  ) %>% 
  group_by(
    COHORT,
    IMAGE
  ) %>% 
  mutate(
    PCT_RECRUITMENT = N_CATEGORY_DWELL_TIME/sum(N_CATEGORY_DWELL_TIME),
    DATE = strsplit(IMAGE, " ")[[1]][1],
    SHAPE = fcase(
      DATE == "20230405", 21,
      DATE == "20230413", 22,
      DATE == "20220516", 23,
      DATE == "20220610", 24,
      DATE == "20220615", 25
    )
  ) %>%
  as.data.table()

Mean_Total <- 
  Mean_Replicates %>% 
  group_by(
    COHORT,
    CATEGORY_DWELL_TIME
  ) %>% 
  summarise(
    SEM_PCT_RECRUITMENT = sem(PCT_RECRUITMENT),
    PCT_RECRUITMENT = mean(PCT_RECRUITMENT)
  ) %>% 
  as.data.table()

# Violin plot --------------------------------------------------------
ggplot(
  data = Mean_Cell,
  aes(
    x = COHORT,
    y = PCT_RECRUITMENT*100)
)+
  geom_violin(
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT*100,
      fill = COHORT,
      group = COHORT
    ),
    alpha = 0.7,
    scale = "width",
    linewidth = 0.2,
    color = "black"
  )+
  geom_boxplot(
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT*100,
      fill = COHORT,
      group = COHORT
    ),
    width = .2,
    outlier.shape = NA,
    fill = NA,
    linewidth = 0.4,
    color = "black"
  )+
  geom_point(
    data = Mean_Replicates,
    aes(
      x = COHORT,
      y = PCT_RECRUITMENT*100,
      group = COHORT
    ),
    position = position_jitter(height=0.3, width=0),
    size = 0.75
  )+
  # stat_compare_means(
  #   data = Mean_Replicates,
  #   comparisons = my_comparison,
  #   hide.ns = TRUE,
  #   label = "p.signif",
  #   tip.length = 0
  # )+
  geom_pwc(
    data = Mean_Replicates,
    method = "wilcox.test",
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")),
    label = "p.signif",
    tip.length = 0.01,
    vjust = 0.5,
    hide.ns = "p"
  )+
  labs(
    y = "% recruitments per cell"
    )+
  scale_x_discrete(
    labels = c("WT", "cMyD88", bquote(cMyD88^TIR), bquote(cMyD88^bDLD), bquote(cMyD88^DHF91))
  )+
  scale_y_continuous(
    breaks = c(0,25,50,75,100)
  )+
  facet_rep_wrap(
    ~CATEGORY_DWELL_TIME,
    ncol = 4
  )+
  fill_palette(
    palette = color_violin
  )+
  theme_classic(base_size = 9
                )+
  theme(
    legend.position = "0",
    plot.title = element_text(hjust = -0.2, vjust = -7, size = 9),
    #strip.text = element_blank(),
    axis.text.x = element_text(color = "black",
                             size = 7,
                             angle = 45,
                             vjust = 0.5,),
    axis.title.x = element_blank(),
    strip.background = element_blank()
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure S2")

ggsave(
  "cl069_cl232_cl236_cl240_cl321_TRAF6-LT_all-recruitments_violin.pdf",
  plot = last_plot(),
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 80,
  width = 184
)
