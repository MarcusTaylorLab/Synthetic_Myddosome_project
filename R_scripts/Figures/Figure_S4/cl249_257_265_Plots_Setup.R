library(pacman)

pacman::p_load(ggplot2, data.table, dplyr, ggfx, viridis, RColorBrewer, ggpubr, lemon, ggbreak, tidyfast, ggbeeswarm)
filter <- dplyr::filter

#Set paht of analyzed Images
TablePaths <- c(
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-1xA TRAF6/20230317 1.5nM_cl249_TRAF6_MyD88-GFP-synTRAF6-BD-1xA 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-3xA TRAF6/20230317 1.5nM_cl257_TRAF6_MyD88-GFP-synTRAF6-BD-3xA 001/Essential.csv.gz",
  "/Volumes/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-5xA TRAF6/20230317 1.5nM_cl265_TRAF6_MyD88-GFP-synTRAF6-BD-5xA 001/Essential.csv.gz"
)

Table <- lapply(TablePaths, 
                fread)

#fill=FALSE since all images are two color
Table <- rbindlist(Table, fill = FALSE)

#Create the same order for all the following plots
Table$COHORT <- factor(
  Table$COHORT, levels = c("MyD88-GFP-synTRAF6-BD-1xA TRAF6",
                           "MyD88-GFP-synTRAF6-BD-3xA TRAF6",
                           "MyD88-GFP-synTRAF6-BD-5xA TRAF6"
  )
)

#set a color palette for the plots
color_violin<-c(
  "MyD88-GFP-synTRAF6-BD-1xA TRAF6" = "#88CCEE", 
  "MyD88-GFP-synTRAF6-BD-3xA TRAF6" = "#999933",
  "MyD88-GFP-synTRAF6-BD-5xA TRAF6" = "#CC6677"
)

#Standard error of the mean function
sem <- function(x) sd(x)/sqrt(length(x))
