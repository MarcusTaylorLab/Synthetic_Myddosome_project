library(pacman)

pacman::p_load(ggplot2, data.table, dplyr, ggfx, viridis, RColorBrewer, ggpubr, lemon, ggbreak, tidyfast, ggbeeswarm)
filter <- dplyr::filter

#Set paht of analyzed Images
TablePaths <- c(
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-1x TRAF6/20221129 4nM_cl247_TRAF6_MyD88-GFP-synTRAF6-BD-1x 001/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-1x TRAF6/20221207 4nM_cl247_TRAF6_MyD88-GFP-synTRAF6-BD-1x 002/Essential.csv.gz",
  #"//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-1x TRAF6/20221207 4nM_cl247_TRAF6_MyD88-GFP-synTRAF6-BD-1x 006/Essential.csv.gz",
  #"//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-1x TRAF6/20230309 1.5nM_cl247_TRAF6_MyD88-GFP-synTRAF6-BD-1x 002/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-1x TRAF6/20230317 1.5nM_cl247_TRAF6_MyD88-GFP-synTRAF6-BD-1x 001/Essential.csv.gz",
  
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-3x TRAF6/20221123 4nM_cl255_TRAF6_MyD88-GFP-synTRAF6-BD-3x 001/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-3x TRAF6/20221207 4nM_cl255_TRAF6_MyD88-GFP-synTRAF6-BD-3x 001/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-3x TRAF6/20230309 1.5nM_cl255_TRAF6_MyD88-GFP-synTRAF6-BD-3x 001/Essential.csv.gz",
  
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-5x TRAF6/20221123 4nM_cl263_TRAF6_MyD88-GFP-synTRAF6-BD-5x 001/Essential.csv.gz",
  #"//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-5x TRAF6/20221207 4nM_cl263_TRAF6_MyD88-GFP-synTRAF6-BD-5x 001/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-5x TRAF6/20230309 1.5nM_cl263_TRAF6_MyD88-GFP-synTRAF6-BD-5x 003/Essential.csv.gz",
  "//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/data_analysis_good_images/MyD88-GFP-synTRAF6-BD-5x TRAF6/20230317 1.5nM_cl263_TRAF6_MyD88-GFP-synTRAF6-BD-5x 001/Essential.csv.gz"
)
  
Table <- lapply(TablePaths, 
                  fread)
  
#fill=FALSE since all images are two color
Table <- rbindlist(Table, fill = FALSE)
  
#Create the same order for all the following plots
Table$COHORT <- factor(
  Table$COHORT, levels = c("MyD88-GFP-synTRAF6-BD-5x TRAF6",
                           "MyD88-GFP-synTRAF6-BD-3x TRAF6",
                           "MyD88-GFP-synTRAF6-BD-1x TRAF6"
  )
)

#set a color palette for the plots
color_violin<-c(
  "MyD88-GFP-synTRAF6-BD-1x TRAF6" = "#88CCEE", 
  "MyD88-GFP-synTRAF6-BD-3x TRAF6" = "#999933",
  "MyD88-GFP-synTRAF6-BD-5x TRAF6" = "#CC6677"
)

#Standard error of the mean function
sem <- function(x) sd(x)/sqrt(length(x))

# cl247_cl255_cl263_MAX-NORM-INT-MYD88_violin ----------------------------------------

source("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 3/cl247_cl255_cl263_MAX-NORM-INT-MYD88_violin.R")

# cl247_cl255_cl263_MAX-NORM-INT-TRAF6_violin ----------------------------------------

source("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 3/cl247_cl255_cl263_MAX-NORM-INT-TRAF6_violin.R")

# cl247_cl255_cl263_NORM-INT_PCT-TRAF6-REC_path ----------------------------------------

source("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 3/cl247_cl255_cl263_NORM-INT_PCT-TRAF6-REC_path.R")

# cl247_cl255_cl263_TRAF6-LT_point ----------------------------------------

source("//data-tay/TAYLOR-LAB/Synthetic Myddosome Paper/Mock Figures/Figure 3/cl247_cl255_cl263_TRAF6-LT_beeswarm.R")
