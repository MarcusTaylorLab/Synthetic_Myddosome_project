# Recruitment of TRAF6 over Normalized Intensity --------------------------
Recruitment_List <- 
  Table %>%
  filter(
    NORMALIZED_INTENSITY >= 1, #filter out the noise
    FRAMES_ADJUSTED <=100, #only take the frames without excessive bleaching
    FRAMES_SINCE_LANDING <= 200,
    PROTEIN == "MyD88"
  ) %>% 
  mutate(
    RECRUITMENT = fcase(
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 >= 1.5, "1",
      COMPLEMENTARY_NORMALIZED_INTENSITY_1 < 1.5, "0" #If TRAF6 is colocalized RECRUITMENT will be 1, otherwise 0
    ),
    NORMALIZED_INTENSITY = 2*round(NORMALIZED_INTENSITY/2), #so we round to even integers
    COHORT = fcase(
      COHORT == "MyD88 TRAF6", "MyD88",
      COHORT == "MyD88-TRAF6-BD TRAF6", "cMyD88",
      COHORT == "MyD88-TIR-TRAF6-BD TRAF6", "cMyD88^{TIR}",
      COHORT == "BDLD_57H-MyD88-TIR-TRAF6-BD-GFP TRAF6", "cMyD88^{bDLD}",
      COHORT == "MyD88-DHF91-TRAF6-BD TRAF6", "cMyD88^{DHF91}"
    )
  ) %>% 
  group_by(
    COHORT,
    NORMALIZED_INTENSITY,
    IMAGE
  ) %>% 
  filter(
    n() >= 10 #so we only look at intensities where there are at least 10 events
  ) %>% 
  summarise(
    NORMALIZED_RECRUITMENT = (sum(RECRUITMENT == 1)/n()),
    SEM_NORMALIZED_RECRUITMENT = sem(RECRUITMENT),
    SD_NORMALIZED_RECRUITMENT = sd(RECRUITMENT)
  ) %>% 
  as.data.table()

Recruitment_List$COHORT <- factor(
  Recruitment_List$COHORT, levels = c(
    "MyD88",
    "cMyD88",
    "cMyD88^{TIR}",
    "cMyD88^{bDLD}",
    "cMyD88^{DHF91}"
  )
)

Mean_of_Means <-
  Recruitment_List %>% 
  group_by(
    COHORT,
    NORMALIZED_INTENSITY
  ) %>% 
  filter(
    n() >= 2 #so we only look at intensities where there are at least 2 replicates
  ) %>%
  summarise(
    counts = length(NORMALIZED_RECRUITMENT),
    SEM_NORMALIZED_RECRUITMENT = sem(NORMALIZED_RECRUITMENT),
    SD_NORMALIZED_RECRUITMENT = sd(NORMALIZED_RECRUITMENT),
    NORMALIZED_RECRUITMENT = mean(NORMALIZED_RECRUITMENT)
  )

color_pal <- 
  c(
    "MyD88" = "#44AA99",
    "cMyD88"= "#117733",
    "cMyD88^{TIR}" = "#332288",
    "cMyD88^{bDLD}" = "#AA4499",
    "cMyD88^{DHF91}" = "#882255"
  )

#plot the percentage of TRAF6 recruitment over normalized Intensity
ggplot(
)+
  geom_line(
    data = Mean_of_Means,
    aes(
      x = NORMALIZED_INTENSITY,
      y = NORMALIZED_RECRUITMENT*100,
      color = COHORT
    ),
    alpha = 1,
    linewidth = 1,
  )+
  geom_line(
    data = Recruitment_List,
    aes(
      x = NORMALIZED_INTENSITY,
      y = NORMALIZED_RECRUITMENT*100,
      group = IMAGE,
      color = COHORT
    ),
    alpha = 0.6,
    linewidth = 0.5
  )+
  color_palette(
    palette = color_pal
  )+
  fill_palette(
    palette = color_pal
  )+
  facet_rep_grid(
    .~ COHORT, scales = "free_x", space = "free_x", labeller = label_parsed
  )+
  ggh4x::facetted_pos_scales(
    x = list(
    COHORT == "MyD88" ~ scale_x_continuous(breaks = c(0, 10), expand = c(0.5,0)),
    COHORT == "cMyD88" ~ scale_x_continuous(breaks = c(0, 50, 100, 150)),
    COHORT == "cMyD88^{TIR}" ~ scale_x_continuous(breaks = c(0, 15, 30), expand = c(0.1,0)),
    COHORT == "cMyD88^{bDLD}" ~ scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300, 350), limits = c(0,350)),
    COHORT == "cMyD88^{DHF91}" ~ scale_x_continuous(breaks = c(0, 25, 50, 75, 100))
  )
  )+
  labs(
    x = "Size of chimeric oligomer",
    y = "% of TRAF6 \n pos. oligomers"
  )+
  theme_classic(base_size = 9)+
  theme(
    legend.position = "0",
    axis.text = element_text(color = "black",
                             size = 7),
    legend.title = element_blank(),
    strip.background = element_blank()
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure S4")

ggsave(
  "cl069_cl232_cl236_cl240_cl321_NORM-INT_PCT-TRAF6-REC-replicates_path.pdf",
  scale = 1,
  units = "mm",
  family = "Helvetica",
  height = 35,
  width = 184
)
