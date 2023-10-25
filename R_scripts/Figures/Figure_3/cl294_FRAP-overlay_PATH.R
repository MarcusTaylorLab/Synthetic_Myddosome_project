library(pacman)

pacman::p_load(ggplot2, ggdark, data.table, dplyr, ggfx, viridis, ggridges, RColorBrewer, ggpubr, lemon)
filter <- dplyr::filter

#Standard error of the mean function
sem <- function(x) sd(x)/sqrt(length(x))

#load tables
TablePaths <- c(
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl294/20230613_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl294/20230928_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl294/20231005_Intensity.csv",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/FRAP_data_analysis/cl294/20231010_Intensity.csv"
)


Table <- lapply(TablePaths, 
                fread)

Table <- rbindlist(Table,
                   fill = TRUE)

Tracks <-
  Table %>% 
  group_by(
    IMAGE,
    CELL
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
    MEAN = mean(RELATIVE_INTENSITY_ADJUSTED),
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
    MEANOFM = mean(MEAN),
    SDDEV = sqrt(mean(VARIANCE)),
    SEM = sqrt(mean(VARIANCE))/sqrt(length(VARIANCE))
  )

#Plot the data
ggplot(
  data = MeanofMeans,
  aes(
    y = MEANOFM,
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
  geom_path(
    color = "black",
    linewidth = 0.75
  )+
  geom_ribbon(
    data = MeanofMeans,
    aes(
      x = TIME_SET_0,
      y = MEANOFM,
      ymin = MEANOFM - SEM,
      ymax = MEANOFM + SEM,
    ),
    fill = "black",
    alpha = 0.4
  )+
  labs(
    x = "Time (s)",
    y = "rel. int. (a.u.) ± s.e.m."
  )+
  annotate(
    geom = "text",
    label = "Bleach",
    x = 20,
    y = 1,
    color = "gray48"
  )+
  scale_y_continuous(
    breaks = c(0, 0.5, 1)
  )+
  scale_x_continuous(
    limits = c(NA,120)
  )+ 
  theme_classic(
    base_size = 10
  )+
  theme(
    legend.position = "0",
    axis.text = element_text(color = "black")
  )

setwd("/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 3")

ggsave(
  "cl294_FRAP-overlay_PATH.pdf",
  scale = 1,
  family = "Helvetica",
  units = "mm",
  height = 40,
  width = 60
)

