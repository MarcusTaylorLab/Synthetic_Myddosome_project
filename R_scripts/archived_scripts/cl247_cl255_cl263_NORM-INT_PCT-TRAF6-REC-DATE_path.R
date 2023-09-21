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
    RECRUITMENT = case_when(
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 2 ~ "1",
      TRUE ~ "0" #If TRAF6 is colocalized RECRUITMENT will be 1, otherwise 0
    ),
    NORMALIZED_INTENSITY = 2*round(NORMALIZED_INTENSITY/2) #so we round to even integers
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
    "MyD88-GFP-synTRAF6-BD-1x TRAF6" = "MyD88-sT6BM-1x",
    "MyD88-GFP-synTRAF6-BD-3x TRAF6" = "MyD88-sT6BM-3x",
    "MyD88-GFP-synTRAF6-BD-5x TRAF6" = "MyD88-sT6BM-5x"
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
  geom_path(
    linewidth = 1
  )+
  geom_path(
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
    limits = c(0, 300),
    breaks = scales::breaks_width(50)
  )+
  scale_y_continuous(
    limits = c(0, 100)
  )+ 
  labs(
    x = "Normalized Intensity of chimeric MyD88",
    y = "% of Puncta colocalizing\n with TRAF6"
  )+
  facet_wrap(
    ~ COHORT,
    labeller = as_labeller(cohort_names)
  )+
  theme_classic(base_size = 22)+
  theme(
    legend.position = "bottom"
  )

setwd("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure S3")

ggsave(
  "cl247_cl255_cl263_NORM-INT_PCT-TRAF6-REC-DATE_path.pdf",
  scale = 3,
  units = "mm",
  height = 30,
  width = 60
)