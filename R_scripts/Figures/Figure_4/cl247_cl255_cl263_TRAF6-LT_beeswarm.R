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
      DATE == "20221129" & COHORT == "MyD88-GFP-synTRAF6-BD-1x TRAF6", 21,
      DATE == "20221207" & COHORT == "MyD88-GFP-synTRAF6-BD-1x TRAF6", 22,
      DATE == "20230317" & COHORT == "MyD88-GFP-synTRAF6-BD-1x TRAF6", 24,
      
      DATE == "20221123" & COHORT == "MyD88-GFP-synTRAF6-BD-3x TRAF6", 21,
      DATE == "20221207" & COHORT == "MyD88-GFP-synTRAF6-BD-3x TRAF6", 22,
      DATE == "20230309" & COHORT == "MyD88-GFP-synTRAF6-BD-3x TRAF6", 24,
      
      DATE == "20221123" & COHORT == "MyD88-GFP-synTRAF6-BD-5x TRAF6", 21,
      DATE == "20230309" & COHORT == "MyD88-GFP-synTRAF6-BD-5x TRAF6", 22,
      DATE == "20230317" & COHORT == "MyD88-GFP-synTRAF6-BD-5x TRAF6", 24
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
      DATE == "20221129" & COHORT == "MyD88-GFP-synTRAF6-BD-1x TRAF6", 21,
      DATE == "20221207" & COHORT == "MyD88-GFP-synTRAF6-BD-1x TRAF6", 22,
      DATE == "20230317" & COHORT == "MyD88-GFP-synTRAF6-BD-1x TRAF6", 24,
      
      DATE == "20221123" & COHORT == "MyD88-GFP-synTRAF6-BD-3x TRAF6", 21,
      DATE == "20221207" & COHORT == "MyD88-GFP-synTRAF6-BD-3x TRAF6", 22,
      DATE == "20230309" & COHORT == "MyD88-GFP-synTRAF6-BD-3x TRAF6", 24,
      
      DATE == "20221123" & COHORT == "MyD88-GFP-synTRAF6-BD-5x TRAF6", 21,
      DATE == "20230309" & COHORT == "MyD88-GFP-synTRAF6-BD-5x TRAF6", 22,
      DATE == "20230317" & COHORT == "MyD88-GFP-synTRAF6-BD-5x TRAF6", 24
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

# Just plot the long lived recruitments
Mean_Cell <- Mean_Cell %>% filter(CATEGORY_DWELL_TIME == ">12 s")
Mean_Replicates <- Mean_Replicates %>% filter(CATEGORY_DWELL_TIME == ">12 s")
Mean_Total <- Mean_Total  %>% filter(CATEGORY_DWELL_TIME == ">12 s")

 # Beeswarm plot --------------------------------------------------------
 ggplot(
     data = Mean_Cell,
     aes(
         x = COHORT,
         y = PCT_RECRUITMENT*100,
         fill = COHORT,
         color = COHORT
       )
   )+
     geom_beeswarm(
         aes(
             shape = Mean_Cell$SHAPE
             ),
         dodge.width = 1,
         cex = 2,
         priority = "density",
         size = 1.5,
         alpha = 1
       )+
     geom_point(
         data = Mean_Replicates,
         aes(
             x = COHORT,
             y = PCT_RECRUITMENT*100,
             group = COHORT,
             shape = Mean_Replicates$SHAPE
           ),
         size = 2,
         stroke = 1.5,
         color = "black",
         fill = NA,
         position = position_jitterdodge (
             dodge.width = 1,
             jitter.width = .7,
             seed = 7
           )
       )+
     geom_errorbar(
       data = Mean_Total,
       aes(
         ymin = (PCT_RECRUITMENT-SEM_PCT_RECRUITMENT)*100,
         ymax = (PCT_RECRUITMENT+SEM_PCT_RECRUITMENT)*100
       ),
       position = position_dodge(width = 1),
       width = 0.2,
       color = "black",
       size = 1
     )+
     geom_pwc(
         data = Mean_Replicates,
         method = "t.test",
         label = "p.signif",
         tip.length = 0,
         hide.ns = TRUE,
         label.size = 4,
         bracket.nudge.y = 0.08,
         vjust = 0.5
       )+
     scale_shape_identity(
       )+
     labs(
         y = "% long lived \n recruitments per cell",
       )+
     scale_y_continuous(
         limits = c(0, NA)
       )+
     fill_palette(
         palette = color_violin
       )+
     color_palette(
         palette = color_violin
       )+ 
     theme_classic(base_size = 9)+
     theme(
         legend.position = "0",
         #strip.text = element_blank(),
         axis.text = element_text(color = "black"),
         axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         strip.background = element_blank()
       )

 setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 4")

 ggsave(
     "cl247_cl255_cl263_TRAF6-LT_beeswarm.pdf",
     scale = 1,
     units = "mm",
     height = 35,
     width = 70
   )
 