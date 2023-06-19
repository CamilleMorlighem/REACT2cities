# ------------------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# More specifically, we used as
# - Dependent variable : the Plasmodium falciparum Parasite Rate standardized over the two-to-ten age range (PfPR2-10)
# - Predictors:
#       - pseudo climatic variables (COSMO)
#       - local climate zone (LCZ), each variable is the proportion of coverage of one specific LCZs
#       - land cover and land use (LULC), developed at very high-resolution, each variable is the proportion of a LC/LU
#       - climatic variables
#       - ancillarly variables (NDVI, NDWI, elevation)
# - Model : random forest modelling built with
#       - Spatial cross validation
#       - Recursive Feature Elimination
#       - Spatial optimisation methods (optional)
#
# This script performs a simulation study to test spatial optimisation methods on simulated displaced data. It is implemented
# with the following steps :
#   1) Divide non-DHS data into training and test sets and simulate displacement on training set using DHS procedure
#   2) Generate duplicates for spatial optimisation methods
#   3) Extract covariates for original data, displaced data, and duplicates of displaced data
#   4) RF modelling using spatial optimisation methods
#   5) Compare predictive performance of methods using test set
#--------------------------------------------------------------------------------------------------------------------------

# -------------------------------
#  Loads libraries
# -------------------------------
library(XLConnect)
library(dplyr)
library(sf)
library(raster)
library(rgdal)
library(stringr)
library(parallel)
library(tidyr)
library(tidyverse)
library(pracma)
library(foreach)
library(exactextractr)
library(plyr)
library(foreach)
library(iterators)

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
ssa_cities = list(c("Dar es salaam", "DES", 32737))

# Define criteria for malaria surveys
max_age = 16 ; survey_date = c(2005,2016) ; select_DHS = T

# Define random forest (RF) parameters for the spatial cross validation in mlr package
iters = 2 ; folds = 5 ; reps = 10 # value for the paper

# -------------------------------
#  Loads paths
# -------------------------------

# Find main directory of the project
dirname(rstudioapi::getSourceEditorContext()$path)
Dir = dirname(rstudioapi::getSourceEditorContext()$path)
Dir = sub(pattern='/[^/]*$', "/", x=Dir)

# Define subdirectories
dir_malaria_data = paste0(Dir, "Results/Malaria_data/")
dir_geo_variables = paste(Dir, c("Data/Variables/"), sep="")
dir_covariates = paste(Dir, c("Results/Covariates/"), sep="")
dir_rf_models = paste(Dir, c("Results/RF_modelling/"), sep="")

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

ssa_admin2 = foreach(i=1:length(ssa_cities), .combine = "c") %do% {
    paste0(dir_admin_areas, ssa_cities[i][[1]][1], "/gadm_2.shp")
}


# -------------------------------
# 1) Simulate displacement
# -------------------------------
simulate.displacement(cities = ssa_cities, crs_REACT = crs_project, malaria_data_db = malaria_db, DirOut = dir_malaria_data, Max_up_age = max_age, survey_period = survey_date, extents = ssa_extents, admin_files = ssa_admin2)

# -------------------------------
# 2) Generate duplicates
# -------------------------------
generate.duplicates(cities = ssa_cities, crs_REACT = crs_project, Dir_input_malaria_data = dir_malaria_data, admin_files = ssa_admin2, simulation = T)

# -------------------------------
# 3) Extract covariates
# -------------------------------

for (method in c("m1", "m0", "TM")){

    covariates.extraction(cities=ssa_cities, crs_REACT = crs_project, aim = "training", n.cores = n_cores, Dir_input_malaria_data = dir_malaria_data,
                          Dir_input_var = dir_geo_variables, Dir_out_cov = dir_covariates, method = method, simulation  = T)

}

# -------------------------------
# 4) RF modelling
# -------------------------------

for (method in c("TM", "m0", "m1", "m2")){

    rf.modelling(cities = ssa_cities, n.cores = n_cores, reps_rf = reps, folds_rf = folds, iters_rf = iters, save_result = "yes",
                 Dir_input_cov = dir_covariates, Dir_input_malaria_data = dir_malaria_data, Dir_output_rf = dir_rf_models, method = method, simulation = T)

}

# -------------------------------
# 5) Compare performance
# -------------------------------
validate.simulation(cities = ssa_cities, Dir_input_cov = dir_covariates, Dir_input_malaria_data = dir_malaria_data, Dir_output_rf = dir_rf_models)

