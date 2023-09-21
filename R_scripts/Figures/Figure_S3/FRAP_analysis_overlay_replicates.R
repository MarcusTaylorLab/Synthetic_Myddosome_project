library(pacman)

pacman::p_load(ggplot2, ggdark, data.table, dplyr, ggfx, viridis, ggridges, RColorBrewer, ggpubr, lemon)
filter <- dplyr::filter

#Standard error of the mean function
sem <- function(x) sd(x)/sqrt(length(x))

#load tables
TablePaths <- c(
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl232/20230124_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl232/20230202_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl232/20230425_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl236/20230124_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl236/20230202_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl236/20230425_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl240/20230124_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl240/20230202_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl240/20230425_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl321/20230713_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl321/20230731_Intensity.csv"
  
)

Table <- lapply(TablePaths, 
                fread)

Table <- rbindlist(Table,
                   fill = TRUE)

Tracks <-
  Table %>% 
  group_by(
    IMAGE
  ) %>% 
  mutate(
    MEAN_INT_BG_REMOVED = MEAN_INTENSITY-BACKGROUND_INTENSITY, 
    #I measured an ROI in the background and substract this value from every measured intensity for normalization
    RELATIVE_INTENSITY = (MEAN_INT_BG_REMOVED-MEAN_INT_BG_REMOVED[5])/(MEAN_INT_BG_REMOVED[1]-MEAN_INT_BG_REMOVED[5]),
    #I get the relative Intensity by subtracting the intensity value of row 5 
    #([5]; time of bleaching) from the respective intensity and dividing both by 
    #the difference of "prebleach" and "bleach". This results in prebleach = 1 and 
    #"bleach" = 0. 
    RATIO_REF = ((REF_INTENSITY[1]-BACKGROUND_INTENSITY[1])/(REF_INTENSITY-BACKGROUND_INTENSITY)),
    #I measured the intensity of an additional spot within the same cell that 
    #wasn't FRAPed to get the overall bleaching rate.
    RELATIVE_INTENSITY_ADJUSTED = RATIO_REF*RELATIVE_INTENSITY,
    #By multiplying the ratio of the reference intensity with the relative 
    #intensity you can correct for the overall bleaching.
    # TIME_ADJUSTED = (FRAMES-FRAMES[2])*4
    #whatever your frame rate was
    COHORT = fcase(
      COHORT == "cl232", "cMyD88",
      COHORT == "cl240", "cMyD88^{TIR}",
      COHORT == "cl321", "cMyD88^{bDLD}",
      COHORT == "cl236", "cMyD88^{DHF91}"
    )
  )

Tracks$COHORT <- factor(
  Tracks$COHORT, levels = c(
    "cMyD88",
    "cMyD88^{TIR}",
    "cMyD88^{bDLD}",
    "cMyD88^{DHF91}"
    )
)

#make a table with the individual cell values
Cells <- 
  Tracks %>%
  mutate(
    DATE = substring(IMAGE, 1, 8),
    TIME_ADJUSTED = (FRAMES-FRAMES[2])*4,
    TIME_SET_0 = TIME_ADJUSTED - TIME_ADJUSTED[5]
  ) %>% 
  filter(
    TIME_SET_0 <= 136
  )

#making a table with the means summarized by cohort and day
Means <- 
  Tracks %>%
  mutate(
    DATE = substring(IMAGE, 1, 8),
    TIME_ADJUSTED = (FRAMES-FRAMES[2])*4,
    TIME_SET_0 = TIME_ADJUSTED - TIME_ADJUSTED[5]
  ) %>% 
  filter(
    TIME_SET_0 <= 136
  ) %>% 
  group_by(
    DATE,
    COHORT,
    TIME_SET_0
  ) %>%
  summarize(
    RELATIVE_INTENSITY_ADJUSTED = mean(RELATIVE_INTENSITY_ADJUSTED),
    VARIANCE = var(RELATIVE_INTENSITY_ADJUSTED)
  )


#means 
MeanofMeans <- 
  Means %>% 
  group_by(
    COHORT,
    TIME_SET_0
  ) %>% 
  summarize(
    RELATIVE_INTENSITY_ADJUSTED = mean(RELATIVE_INTENSITY_ADJUSTED),
    SDDEV = sqrt(mean(VARIANCE)),
    SEM = sqrt(mean(VARIANCE))/sqrt(length(VARIANCE))
  )

color_pal <- 
  c(
    "cMyD88"= "#117733",
    "cMyD88^{TIR}" = "#332288",
    "cMyD88^{bDLD}" = "#AA4499",
    "cMyD88^{DHF91}" = "#882255"
  )

#Plot the data
ggplot(
  data = MeanofMeans,
  aes(
    y = RELATIVE_INTENSITY_ADJUSTED,
    x = TIME_SET_0
  )
)+
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = 0,
      ymin = -Inf,
      ymax = Inf
    ),
    alpha = 0.1,
    fill = "gray88"
  )+
  geom_vline(
    color = "gray48",
    linewidth = 0.9,
    linetype = 3,
    aes(
      xintercept = 0
    )
  )+
  geom_line(
    data = Cells,
    aes(
      y = RELATIVE_INTENSITY_ADJUSTED,
      x = TIME_SET_0,
      group = IMAGE
    ),
    color = "grey50"
  )+
  geom_line(
    data = MeanofMeans,
    aes(
      y = RELATIVE_INTENSITY_ADJUSTED,
      x = TIME_SET_0,
      color = COHORT
    ),
    linewidth = 1.5
  )+
  facet_wrap(
    ~COHORT,
    ncol = 4,
    labeller = label_parsed
  )+
  labs(
    x = "Time (s)",
    y = "Relative intensity (a.u.) ± s.e.m."
  )+
  annotate(
    geom = "text",
    label = "Bleach",
    x = 25,
    y = 1,
    color = "gray48"
  )+
  scale_y_continuous(
    breaks = c(0, 0.5, 1)
  )+
  scale_x_continuous(
    limits = c(NA,120)
  )+
  fill_palette(
    palette = color_pal
  )+
  color_palette(
    palette = color_pal
  )+ 
  theme_classic(
    base_size = 10
  )+
  theme(
    legend.position = "0",
    axis.text = element_text(colour = "black"),
    strip.background = element_blank()
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure S3")

ggsave(
  "cl232_cl236_cl240_cl231_FRAP-overlay_replicates_PATH.pdf",
  scale = 1,
  units = "mm",
  height = 50,
  width = 184
)

