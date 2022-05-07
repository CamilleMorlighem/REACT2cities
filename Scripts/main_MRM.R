# ------------------------------------------------------------------------------------------------------------------
# Modelling of the intra-urban malaria risk using socioeconomic and environmental factors with a random forest (RF) model
# in two sub-Saharan cities: Kampala (Uganda) and Dar es Salaam (Tanzania)
#  
# More specifically, we used as 
# - Dependent variable : a malaria prevalence metric: the Plasmodium falciparum Parasite Rate standardized over the two-to-ten age range (PfPR2-10)
# - Predictors: 
#       - Normalized Difference Vegetation Index (NDVI) and Normalized Difference Wetness Index (NDWI)
#       - Elevation (Shuttle Radar Topography Mission (SRTM))
#       - variables from  different geospatial datasets: 
#             - pseudo climatic variables (COSMO) 
#             - local climate zone (LCZ), each variable is the proportion of coverage of one specific LCZs 
#             - land cover and land use (LULC), developed at very high-resolution, each variable is the proportion of a LC/LU
#             
# 
# This script performs all the following steps to model and predict PfPR2-10: 
#   1) select malaria surveys 
#   2) extract covariates 
#   3) RF modelling 
#   4) RF predictions
#--------------------------------------------------------------------------------------------------------------------------


rm(list=ls())

# -------------------------------
#  Loads libraries 
# -------------------------------
library(XLConnect)
library(sf)
library(utils)      #Malaria data 
library(foreach)
library(raster)
# -------------------------------
library(spex)       #Prediction grid 
# -------------------------------
library(dplyr)
library(exactextractr)
library(parallel)
library(doParallel) #Covariates extraction 
library(iterators)
library(sp)
library(units)
library(plyr)
# -------------------------------
library(stats)
library(corrplot)   #Descriptive analysis
# -------------------------------
library(ggplot2) 
library(mlr)        #RF modelling 
library(parallelMap) 
library(pdp)
# -------------------------------
library(ggpubr)     #BoxPlot results 
# -------------------------------
library(stars)      #RF predictions


# -------------------------------
#  Parameters user can change 
# -------------------------------

# Define parallelization parameters 
parallel::detectCores()
n_cores <- parallel::detectCores() - 3
my.cluster <- parallel::makeCluster(n_cores, type = "PSOCK")

# Register cluster to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

# Define cities 
ssa_cities = list(c("Kampala", "Kamp", 32636), c("Dar es salaam", "DES", 32737)) 

# Define criteria for malaria surveys 
max_age = 16 ; survey_date = c(2005,2016) ; select_DHS = F

# Define aim(s) of covariates extraction: 1- "training" ; 2- "prediction" (both can be used)
cov_extract_aims = c( "prediction")

# Define mode:   1- All covariates without variable selection (ALL) ; 2- Implement a Variable selection : Recursive features selection (RFE)
cov_selection_mode = "ALL"

# Define random forest (RF) parameters for the spatial cross validation in mlr package
iters = 2 ; folds = 5 ; reps = 10 # value for the paper 

# Define var group :  1 - All covariates from 3 geospatial datasets (all3Geo) - LU/LC, LCZ, COSMO ; 2 - Comparison of different geospatial datasets (GeoSpDt) ; 
var_group_list = c("all3Geo", "GeoSpDt") ; variables_group = var_group_list[1]  

# Define pred mode for all3Geo:   1- All covariates without variable selection (ALL) ; 2- Recursive features selection model (RFE) with iteration i to be used ("RFE-[i]")
pred_mode = "ALL"


# -------------------------------
#  Loads paths
# -------------------------------

# Find main directory of the project 
dirname(rstudioapi::getSourceEditorContext()$path)
Dir = dirname(rstudioapi::getSourceEditorContext()$path)
Dir = sub(pattern='/[^/]*$', "/", x=Dir)

# Define subdirectories 
dir_malaria_data = paste0(Dir, "Results/Malaria_data/")
dir_pred_grid = paste(Dir, c("Data/Prediction_grid/"), sep="")
dir_geo_variables = paste(Dir, c("Data/Variables/"), sep="")
dir_covariates = paste(Dir, c("Results/Covariates/"), sep="")
dir_rf_models = paste(Dir, c("Results/RF_modelling/"), sep="")
dir_rf_pred =  paste0(Dir, c("Results/Prediction/"))
dir_admin_areas = paste0(Dir, c("Data/Admin/"))

# Define path to malaria DB 
malaria_db = paste0(Dir, "Data/Malaria_data/Malaria_data_v1_SDJ_20200520.xlsx") 

# Load functions 
invisible(sapply(list.functions<-list.files(file.path(paste0(Dir,"/Scripts/Functions/")), full.names=TRUE), source))


# -------------------------------
#  Fixed parameters 
# -------------------------------

# Fix the crs of the project 
crs_project <- CRS("+init=epsg:4326")

# Fix cities extent : LU map is the smaller extent of all variables --> defines the extent of city 
ssa_extents = foreach(i=1:length(ssa_cities), .combine = "c") %do% {
  paste0(dir_geo_variables, "LULC/", ssa_cities[i][[1]][1], "/LU/", ssa_cities[i][[1]][2], "_LU.shp")
}


# -------------------------------
# Select malaria data 
# -------------------------------

select.malaria.data(cities = ssa_cities, crs_REACT = crs_project, malaria_data_db = malaria_db, DirOut = dir_malaria_data,
                    Max_up_age = max_age, survey_period = survey_date, extents = ssa_extents, use_DHS = select_DHS)


# -------------------------------
# Covariates extraction 
# -------------------------------

foreach(aim_i = 1:length(cov_extract_aims)) %do% {

  covariates.extraction(cities=ssa_cities, crs_REACT = crs_project, aim = cov_extract_aims[aim_i], n.cores = n_cores, Dir_input_malaria_data = dir_malaria_data,
                        Dir_input_var = dir_geo_variables, Dir_out_cov = dir_covariates, Dir_input_pts = dir_pred_grid)

}

# -------------------------------
# RF modelling 
# -------------------------------

rf.modelling(cities = ssa_cities, my_cores = n_cores, mode = cov_selection_mode, reps_rf = reps, folds_rf = folds, iters_rf = iters,
             var_group = variables_group, input_cov = dir_covariates, input_malaria = dir_malaria_data, output_file_own = dir_rf_models)


# -------------------------------
# RF predictions 
# -------------------------------

rf.prediction(cities = ssa_cities, crs_REACT = crs_project, input_cov = dir_covariates, input_malaria = dir_malaria_data, output_rf = dir_rf_models, output_pred = dir_rf_pred, input_admin = dir_admin_areas, Dir_input_pts = dir_pred_grid, var_group = "all3Geo", mode = pred_mode)
  

# -------------------------------
# Opt: Descriptive analysis
# -------------------------------

check.correlation(cities = ssa_cities, crs_REACT = crs_project, input_cov = dir_covariates, input_malaria = dir_malaria_data)

