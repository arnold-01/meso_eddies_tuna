The model is a GLM-GAMs model developed in R using the libraries (library(MASS), library(mgcv), library(ggplot2), library(RODBC), library(plotrix)).
jl for the species albacore, bigyes tuna, swordfish, and yellowfin tuna. The model is developed separately for the four species (glm_alb_pred, glm_bet_pred, glm_swo_pred, and glm_yft_pred).
Install the latest version of the R programming language on your system (https://cran.r-project.org/doc/manuals/r-patched/R-admin.html) and virtual studio (jupyter, python) (https://code.visualstudio.com/download).
The simulations and execution code are available in the files fusion_data_GLM_GAMs.R, samplaing_area.ipynb, prediction_description_variable.ipynb, and temporal_spatial_distribution.ipynb.
The anticyclonic_2001_2022.csv and cyclonic_parametres_2001_2022.csv files contain annual, monthly, and daily mesoscale-eddy parameter data in the BCLME, southeast Atlantic Ocean.
The tuna01_22zaf.csv and tuna01_22nam.csv files contain catch/effort data for four pelagic species (albacore, bigeye tuna, swordfish, and yellowfin tuna) from January-December 2001 to 2022, respectively, in the South African zone and the Namibian zone.
The tuna_data_namzaf0122.csv files are a merger of (tuna01_22zaf.csv and tuna01_22nam.csv) to create a single study area. 
The tuna_data_clean.csv file is the fusion of all prediction data (mesoscale-eddies parameters) and standardization of the four pelagic species, according to factors (year, month, longitude, latitude).
The complete documentation is available on request .MIT License, Copyright (c) 2025 Arnold EKOUE
