# Recruitment of TRAF6 over Normalized Intensity --------------------------
Recruitment_List_data <- 
  Table %>%
  filter(
    NORMALIZED_INTENSITY >= 1, #filter out the noise
    FRAMES_ADJUSTED <=100, #only take the frames without excessive bleaching
    FRAMES_SINCE_LANDING <= 200,
    PROTEIN == "MyD88"
  ) %>% 
  mutate(
    RECRUITMENT = fcase(
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 2, "1",
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 <2, "0" #If TRAF6 is colocalized RECRUITMENT will be 1, otherwise 0
    ),
    NORMALIZED_INTENSITY = round(NORMALIZED_INTENSITY) #so we round to even integers
  ) %>% 
  group_by(
    COHORT,
    NORMALIZED_INTENSITY
  ) %>% 
  filter(
    n() >= 25 #so we only look at intensities where there are at least 25 events
  ) 

Recruitment_List_cohort <-
  Recruitment_List_data %>% 
  group_by(
    COHORT,
    NORMALIZED_INTENSITY
  ) %>%
  filter(
    n() >= 3 #so we only look at intensities where there are at least 3 events
  ) %>% 
  summarise(
    NORMALIZED_RECRUITMENT = (sum(RECRUITMENT == 1)/n())
  ) %>% 
  as.data.table()

Recruitment_List_Day <-
  Recruitment_List_data %>% 
  group_by(
    COHORT,
    IMAGE,
    NORMALIZED_INTENSITY
  ) %>% 
  filter(
    n() >= 3 #so we only look at intensities where there are at least 3 events
  ) %>% 
  summarise(
    NORMALIZED_RECRUITMENT = (sum(RECRUITMENT == 1)/n())
  ) %>%
  mutate(
    DATE = strsplit(IMAGE, " ")[[1]][1]
  ) %>% 
  as.data.table()

#rename COHORTS to a shorter character
cohort_names <- 
  c(
    "MyD88-TRAF6-BD TRAF6" = "MyD88-T6BM", 
    "MyD88 TRAF6" = "MyD88",
    "MyD88-TIR-TRAF6-BD TRAF6" = "TIR-T6BM",
    "MyD88-DHF91-TRAF6-BD TRAF6" = "DHF91-TIR-T6BM",
    "BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6" = "bDLD-TIR-T6BM"
)

#plot the percentage of TRAF6 recruitment over normalized Intensity
ggplot(
  data = Recruitment_List_Day,
  aes(
    x = NORMALIZED_INTENSITY,
    y = NORMALIZED_RECRUITMENT*100,
    group = DATE,
    color = DATE
  )
)+
  geom_line(
    linewidth = 1
  )+
  geom_line(
    data = Recruitment_List_cohort,
    aes(
      x = NORMALIZED_INTENSITY,
      y = NORMALIZED_RECRUITMENT*100,
      group = COHORT
    ),
    linewidth = 1,
    color = "black"
  )+
  scale_x_continuous(
    limits = c(0, 100),
    breaks = scales::breaks_width(20)
  )+
  scale_y_continuous(
    limits = c(0, 100)
  )+ 
  labs(
    x = "Size of chimeric oligomer",
    y = "% of oligomers \n colocalizing with TRAF6"
  )+
  facet_wrap(
    ~ COHORT,
   labeller = as_labeller(cohort_names)
  )+
  theme_classic(base_size = 22)+
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure S2")

ggsave(
  "cl247_cl255_cl263_NORM-INT_PCT-TRAF6-REC-DATE_path.pdf",
  scale = 3,
  units = "mm",
  height = 30,
  width = 60
)