# Structure of repository

This folder contains R_scripts which are used to generate plots from the output of the image analysis pipeline. 
The repository is made up of three main sections. 

- ```Figures```: These are the up-to-date scripts which are used to generate all R-plots for the different figures
- ```Scripts_ELISA```: These scripts are used to generate tables with the correct format for further R analysis from the raw ELISA data
- ```archived_scripts```: These scripts were previously used to generate plots for figures, but are currently not used anymore
- ```tempate_scripts```: Some plots are being used multiple times with different cell lines for different figures. These scripts can be used as templates for future analysis.

## Detailed information on structure of ```Figures``` folder

This folder contains a folder each for a main and a supplementary figure. The script ```color_palette.R``` contains the HEX codes of the color scheme used across the whole paper.

The nomenclautre of naming the files is as follows: 
**cell-lines-used-in-plot**_**x-axis**_**y-axis**_**type-of-plot**.R

and exemplified by these two examples:
```cl247_cl255_cl263_MAX-NORM-INT-TRAF6_violin.R```
```cl069_cl232_cl236_cl240_cl321_NORM-INT_PCT-TRAF6-REC_path.R```

For each Figure there is a ```_setup.R``` (i.e. ```cl249_257_265_Plots_Setup.R``` script in the corresponding folder. This script loads all the required packages for the plots, defines functions needed for the analysis and reads the ```Essential.csv.gz``` files for each image of the respective cell line. The ```Essential.csv.gz``` is the most important output of the image analysis pipeline as it combines all the measured values in one table.
By using ```_setup.R``` the repetition of code is greatly reduced, as well as the time consuming reading of the ```Essential.csv.gz``` files.
