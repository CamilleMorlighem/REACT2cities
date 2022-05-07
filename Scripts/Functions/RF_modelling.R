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
# This function includes the steps to assess different RF models with different combination of variables 
# as well as different variable selection (RFE or not) 
# using a RF model built in spatial cross validation (with mlr and range R packages)
#
# It create for each city: 
# - "RFCovImp... .csv" : a table (csv file) with the RF results (Goodness-of-fit indices, RF parameters, etc.), 
#    each line containes the result for one simulation of the spatial cross-validation  
# - a table (csv file) with the importance for each covariates for each simulation coming from the spatial cross-validation 
# - a barplot (pdf file) with the average of covariates importance (across all simulations)
# - partial dependence plots (pdf files) for each covariate included in the model (each line grey results form a simulation)
# - rds files with the partial dependence values
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities                 List with one vector c(city name, city short name, city epsg of projected crs) per city 
#' @param my_cores               Number of cores used for parallelelization 
#' @param mode                   Either "ALL" (all covariates are used) or "RFE" (implements a variable selection). Default is NULL 
#' @param reps_rf                Number of repetitions in the spatial cross-validation (default 10)
#' @param folds_rf               Number of folds used in the spatial cross-validation (default 5)
#' @param iters_rf               Number of sub-folds used for hyperparameters tuning inside the spatial cross-validation (default 2)
#' @param var_group              Either "all3Geo" (uses all LULC, LCZ and COSMO covariates), "BEST" (uses best selected covariates) or "GeoSpDt" (compares groups of covariates)
#' @param best_var               If var_group is "BEST", list with one vector/city containing the covariates in the city best RFE model. Default is NULL
#' @param save_result            "yes" or "no" (default "yes")  
#' @param input_cov              Path to input covariate tables 
#' @param input_malaria          Path to input selected malaria data points 
#' @param output_file_own        Path to output RF models

#' @return                       Nothing 
#' @export                       # Per city: "RFCovImp... .csv" : a table with RF results
#                                             A csv with covariate importance for each simulation coming from the spatial cross-validation 
#                                             A barplot (pdf) with the average covariates importance (across all simulations)
#                                             Partial dependence plots (pdf) for each covariate included in the model 
#                                             rds files with the partial dependence values
#--------------------------------------------------------------------------------------------------------------------------


rf.modelling = function(cities, my_cores, mode = NULL, reps_rf = 10, folds_rf = 5, iters_rf = 2, var_group, save_result = 'yes', input_cov, input_malaria, output_file_own, best_var = NULL){
    
  # Define file path  
  file_path = paste0(output_file_own,var_group,"/")  
  dir.create(file_path, showWarnings = F, recursive = T)
  
  # ------------------------------------------------------------------------------------------------------------------
  #                 A. Modelling step 
  # ------------------------------------------------------------------------------------------------------------------
  
  # Choose to implement:  - a spatial cross-validation with a Hyperparameter tuning within the sp cv (SpCVHpt)
  #                       - tuning within the folds of the Spatial CV 2 (HYt_insideCV), by dividing the each fold by iters_rf number                        
  # ------------------------------------------------------------------------------------------------------------------
  cv_rf = "SpCVHpt"
  HYt_mode = "HYt_insideCV"
  
  
  # ------------------------------------------------------------------------------------------------------------------
  # A.a  Modelling part to compare different geospatial datasets "GeoSpDt"
  # ------------------------------------------------------------------------------------------------------------------
  
  # For loop on the different models to compare (models included in GeoSpDt_list)
  # ------------------------------------------------------------------------------------------------------------------
  if(var_group == "GeoSpDt"){
    # Define the for loop nb to start with 
    if(mode == "ALL"){  start_GeoSpDt_list = 1 
    } else if(mode == "RFE"){  start_GeoSpDt_list = 2 }
    
    # Define models to compare (var group list)
    # --------------------------------------------------------------------------------------------------------------------
    GeoSpDt_list = list(
      # base model to compare with the other models
      list(var_group_test= "Base", var_list=c("PfPR2_10", "AVG_NDVI_100", "AVG_NDWI_100", "SRTM_30")),
      
      # LCZdp - LCZ distance & propotion
      list(var_group_test= "LCZdp", var_list=c("PfPR2_10", "AVG_NDVI_100", "AVG_NDWI_100", "SRTM_30",
                                               "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland", #"LCZ_mineral",
                                               "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                               
                                               "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland", #"Dist_mineral",
                                               "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands")),
      
      list(var_group_test = "LCZdp_COSMO", var_list=c("PfPR2_10", "AVG_NDVI_100", "AVG_NDWI_100", "SRTM_30",
                                                      "AVG_QV2M", "AVG_RH2M", "AVG_T2M", "AVG_TS", "AVG_PREC",
                                                      "MAX_PREC", "MIN_PREC",
                                                      "AVG_TSI", "AVG_TSI_RH",
                                                      "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                      "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                      "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland",
                                                      "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands")),
      
      list(var_group_test = "LCZdp_LULC", var_list=c("PfPR2_10", "AVG_NDVI_100", "AVG_NDWI_100", "SRTM_30",
                                                     
                                                     "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                     "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                     "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland",
                                                     "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands",
                                                     
                                                     "LC_bare_ground", "LC_buildings","LC_low_veg"  ,"LC_trees" , "LC_water" , "LU_ACS" ,
                                                     "LU_wetlands","LU_planned","LU_informal")),
      
      list(var_group_test = "LCZdp_LULC_COSMO", var_list=c("PfPR2_10", "AVG_NDVI_100", "AVG_NDWI_100", "SRTM_30",
                                                           
                                                           "LC_bare_ground", "LC_buildings","LC_low_veg"  ,"LC_trees" , "LC_water" , "LU_ACS" ,
                                                           "LU_wetlands","LU_planned","LU_informal",
                                                           
                                                           "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                           "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                           "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland",  
                                                           "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands",
                                                           
                                                           "AVG_QV2M", "AVG_RH2M", "AVG_T2M", "AVG_TS",
                                                           "AVG_PREC", "MAX_PREC", "MIN_PREC",  "AVG_TSI", "AVG_TSI_RH"))
    )
    
    # For loop on the different models to compare
    # --------------------------------------------------------------------------------------------------------------------
    for (var_group_test_i in start_GeoSpDt_list:length(GeoSpDt_list)){     # var_group_test_i = 1
      
      # Define model type and list of variables
      # -------------------------------------------------------------------------------------------------------------------
      GeoSpDt_modtype = GeoSpDt_list[[var_group_test_i]]$var_group_test ; GeoSpDt_modtype
      var_group_test_list = GeoSpDt_list[[var_group_test_i]]$var_list 
      
      # For loop on cities
      # -------------------------------------------------------------------------------------------------------------------
      for (city_nb in 1:length(cities)){
        
        # Define city & input file path & output file name
        # ------------------------------------------------------------------------------------------------------------------
        city = cities[city_nb][[1]][1]
        city_short = cities[city_nb][[1]][2]
        
        # Define (and create) output file path according to the model type
        # -------------------------------------------------------------------------------------------------------------------
        output_file = paste0(file_path,city, '/', mode,"/",GeoSpDt_modtype,"/") ; output_file
        
        # Check for file existence and create it if needed
        if(!dir.exists(file.path(output_file))){ dir.create(file.path(output_file), recursive = TRUE)   }
        
        output_file_name = paste0(cv_rf,"_",var_group,"_",mode,"_",GeoSpDt_modtype,"_",city) ; output_file_name 
  
        # Load Malaria and covariates data
        # ------------------------------------------------------------------------------------------------------------------
        # Load shapefile of points with malaria prevalence, as 'sf' object
        malaria_sf = st_read(paste0(input_malaria, city, "/malaria_data_",city_short,".shp"))
        malaria_sf <- malaria_sf[,c("ID","PfPR2_10")]
        malaria_sf$ID = as.numeric(as.character(malaria_sf$ID))
        
        # store the coordinates of malaria data
        coords = as.data.frame(cbind(st_coordinates(malaria_sf), ID = malaria_sf$ID))
        
        # Load covariates csv into a dataframe
        cov_dataset = read.csv2(paste0(input_cov, city, "/", city_short, "_full_cov.csv"), sep =";")
        cov_dataset = dplyr::left_join(st_drop_geometry(malaria_sf), cov_dataset, by=c("ID"= "ID"))
        
        # ------------------------------------------------------------------------------------------------------------------
        #       Variable selection
        # ------------------------------------------------------------------------------------------------------------------
        # Select variable according the variable list of the model to be assessed
        colnames(cov_dataset)
        variables= dplyr::select(cov_dataset,PfPR2_10, var_group_test_list)
        var_list = colnames(variables) ; var_list
        
        # Verify for NA atfer variable selection - and update the data
        data_updated = na.omit(cov_dataset[,c("ID",var_list)])
        coords_updated = coords[which(coords$ID %in%  data_updated$ID),c("X","Y")]
        variables = data_updated[,var_list]
        
        # define predictor list - variable list without depend variable
        predictor_dataset = dplyr::select(variables,-PfPR2_10)
        predictor_list = colnames(predictor_dataset) ; predictor_list
        
        # ------------------------------------------------------------------------------------------------------------------
        #       RF in Spatial Cross-Validation + hyperparameters (HP) tuning with mlr R package & Ranger package
        #       Covariate importance & partial dependence plot computation
        # ------------------------------------------------------------------------------------------------------------------
        
        # A - No variable seletion - use of all covariates
        # ------------------------------------------------------------------------------------------------------------------
        if(mode == "ALL"){
  
          # self-created 'RF.CovImp.pdp()' function :
          # - Run a random forest (RF) model in spatial cross validation with hyperparameter tuning
          # - store RF results, model characteristics and goodness-of-fit indices in a csv file
          # - Extract covariates importance and store these values store in a csv file
          # - Create a barplot figure with the importance of each covariate - in pdf
          # - Create a dependence plot for each covariate of the model - in pdf
  
          RF_CovImp_pdp = RF.CovImp.pdp(variables = variables, coordinates = coords_updated, dependent_var = "PfPR2_10",
                                          var_group = var_group, city = city, mode = mode,
                                          folds_rf = folds_rf, reps_rf = reps_rf, iters_rf = iters_rf,
                                          nb_cores = my_cores,
                                          output_file = output_file, output_file_name = output_file_name, save_result = save_result, HYt_mode = HYt_mode)
  
  
  
  
        # B - Implement a manual recursive feature selection
        # ------------------------------------------------------------------------------------------------------------------
        } else if(mode == "RFE"){
          for (test_var in 1:(length(predictor_list)-2)){ # test_var = 1
  
            if( length(var_list) > 3 ){
              # Define the RFE step
              nb_test = paste0("RFE-",test_var)
  
              # Define output file name
              output_file_name = paste0(cv_rf,"_",var_group,"_",nb_test,"_",city) ; output_file_name
  
              # RF model, Covariate Imporance & partial dependence plot computation
              # ---------------------------------------------------------------------------------------------------------------
  
              RF_CovImp_pdp = RF.CovImp.pdp(variables = variables, coordinates = coords_updated, dependent_var = "PfPR2_10",
                                              var_group = var_group, city = city, mode = mode, nb_test = nb_test,
                                              folds_rf = folds_rf, reps_rf = reps_rf, iters_rf = iters_rf,
                                              nb_cores = my_cores,
                                              output_file = output_file, output_file_name = output_file_name, save_result = save_result,
                                              HYt_mode = HYt_mode)
  
  
              # Delete the less important covariates
              # ---------------------------------------------------------------------------------------------------------------
              CovImp_summary = RF_CovImp_pdp ; Cov_minImp = CovImp_summary$var[which(CovImp_summary$mean == min(CovImp_summary$mean))]
              variables = dplyr::select(variables, -Cov_minImp)
              var_list = colnames(variables) ; var_list
            }
          }
        }
      }
    }
  
  
  # -----------------------------------------------------------------------------------------------------------------------
  # A.b  Model with all variables from 3 geospatial datasets (all3Geo) 
  # -----------------------------------------------------------------------------------------------------------------------
  
  } else if(var_group == "all3Geo"){
    
    # for loop on cities
    for (city_nb in 1:length(cities)){
      
      # Define city
      # ------------------------------------------------------------------------------------------------------------------
      city = cities[city_nb][[1]][1]
      city_short = cities[city_nb][[1]][2]
      
      # Define output file path according to the model type
      # ------------------------------------------------------------------------------------------------------------------
      output_file = paste0(output_file_own,"/",var_group,"/",city, "/", mode,"/") ; output_file
      
      # Check for file existance and create it if needed
      if(!dir.exists(file.path(output_file))){
        dir.create(file.path(output_file), recursive = TRUE)   }
      output_file_name = paste0(var_group,"_",mode,"_",city) ; output_file_name 
      
      # Load Malaria and covariates data
      # ------------------------------------------------------------------------------------------------------------------
      # Load shapefile of points with malaria prevalence, as 'sf' object
      malaria_sf = st_read(paste0(input_malaria, city, "/malaria_data_",city_short,".shp"))
      malaria_sf <- malaria_sf[,c("ID","PfPR2_10")]
      malaria_sf$ID = as.numeric(as.character(malaria_sf$ID))
      
      # store the coordinates of malaria data
      coords = as.data.frame(cbind(st_coordinates(malaria_sf), ID = malaria_sf$ID))
      
      # Load covariates csv into a dataframe
      cov_dataset = read.csv2(paste0(input_cov, city, "/", city_short, "_full_cov.csv"), sep =";")
      cov_dataset = dplyr::left_join(st_drop_geometry(malaria_sf), cov_dataset, by=c("ID"= "ID"))
      
      # ------------------------------------------------------------------------------------------------------------------
      #       Variable selection
      # ------------------------------------------------------------------------------------------------------------------
    
      var_list=c("PfPR2_10", "AVG_NDVI_100", "AVG_NDWI_100", "SRTM_30",
                 
                 "LC_bare_ground", "LC_buildings","LC_low_veg"  ,"LC_trees" , "LC_water" , "LU_ACS" ,
                 "LU_wetlands","LU_planned","LU_informal",
                 
                 "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                 "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                 "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland",  
                 "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands",
                 
                 "AVG_QV2M", "AVG_RH2M", "AVG_T2M", "AVG_TS",
                 "AVG_PREC", "MAX_PREC", "MIN_PREC",  "AVG_TSI", "AVG_TSI_RH")
        
      
      variables= dplyr::select(cov_dataset,var_list)
      
      # Verify for NA value after covariate selection & update the data
      data_updated = na.omit(cov_dataset[,c("ID",var_list)])
      coords_updated = coords[which(coords$ID %in%  data_updated$ID),c("X","Y")]
      variables = data_updated[,var_list]
      
      # define predictor list and dataset - variable without dependent variable
      predictor_dataset = dplyr::select(variables,-PfPR2_10)
      predictor_list = colnames(predictor_dataset)
      
      # ------------------------------------------------------------------------------------------------------------------
      #       Spatial Cross-Validation + hyperparameters (HP) tuning with mlr R package & Ranger package
      #       Covariate importance & partial dependence plot computation
      # ------------------------------------------------------------------------------------------------------------------
      
      # A - No variable seletion - use of all covariates
      # ------------------------------------------------------------------------------------------------------------------
      if(mode == "ALL"){
        
        # RF model, Covariate Imporance & partial dependence plot computation
        # ---------------------------------------------------------------------------------------------------------------
  
        RF_CovImp_pdp = RF.CovImp.pdp(variables = variables, coordinates = coords_updated, dependent_var = "PfPR2_10",
                                        var_group = var_group, city = city, mode = mode,
                                        folds_rf = folds_rf, reps_rf = reps_rf, iters_rf = iters_rf,
                                        nb_cores = my_cores,
                                        output_file = output_file, output_file_name = output_file_name, save_result = save_result, HYt_mode = HYt_mode)    
        
      
      # B - Implement a manual recursive features selection (RFE)
      # ------------------------------------------------------------------------------------------------------------------
      } else if(mode == "RFE"){
        for (test_var in 1:(length(predictor_list)-1)){
          
          # Define the step
          nb_test = paste0("RFE-",test_var)
          
          # Define output file name
          output_file_name = paste0(var_group,"_",nb_test,"_",city) ; output_file_name
          
          # RF model, Covariate Imporance & partial dependence plot computation
          # ---------------------------------------------------------------------------------------------------------------
  
          RF_CovImp_pdp = RF.CovImp.pdp(variables = variables, coordinates = coords_updated, dependent_var = "PfPR2_10",
                                          var_group = var_group, city = city, mode = mode,
                                          folds_rf = folds_rf, reps_rf = reps_rf, iters_rf = iters_rf,
                                          nb_cores = my_cores, nb_test = nb_test, 
                                          output_file = output_file, output_file_name = output_file_name, save_result = save_result, HYt_mode = HYt_mode)    
          
          # Delete the less important covariate
          # ------------------------------------------------------------------------------------------------------------------
          CovImp_summary = RF_CovImp_pdp ; Cov_minImp = CovImp_summary$var[which(CovImp_summary$mean == min(CovImp_summary$mean))]
          variables = dplyr::select(variables, -Cov_minImp)
          var_list = colnames(variables) ; var_list
        }
      }
    }
  
  
  
  # -----------------------------------------------------------------------------------------------------------------------
  # A.c  Best models selected by RFE
  # -----------------------------------------------------------------------------------------------------------------------
  
  } else if (var_group == "BEST"){
    
    mode ="NULL"
    
    for (city_nb in 1:length(cities)){
      
      # Define city & input file path & output file name
      # ------------------------------------------------------------------------------------------------------------------
      city = cities[city_nb][[1]][1]
      city_short = cities[city_nb][[1]][2]
      
      # define output file and name
      # ------------------------------------------------------------------------------------------------------------------
      output_file_name = paste0(cv_rf,"_",var_group,"_",city) ; output_file_name 
      output_file = paste0(file_path,city,"/") ; output_file
      # Check for file existence and create it if needed
      if(!dir.exists(file.path(output_file))){ dir.create(file.path(output_file), recursive = TRUE)   }
      
      # Load Malaria and covariates data
      # ------------------------------------------------------------------------------------------------------------------
      # Load shapefile of points with malaria prevalence, as 'sf' object
      malaria_sf = st_read(paste0(input_malaria, city, "/malaria_data_",city_short,".shp"))
      malaria_sf <- malaria_sf[,c("ID", "PfPR2_10")] ; malaria_sf$ID = as.numeric(as.character(malaria_sf$ID))
      
      # store the coordinates (lon and lat) of malaria data
      coords = as.data.frame(cbind(st_coordinates(malaria_sf), ID = malaria_sf$ID))
      
      # Load covariates csv into a dataframe
      cov_dataset = read.csv2(paste0(input_cov, city, "/", city_short, "_full_cov.csv"), sep =";")
      cov_dataset = dplyr::left_join(st_drop_geometry(malaria_sf), cov_dataset, by=c("ID"= "ID"))
      
      # ------------------------------------------------------------------------------------------------------------------
      #       Variable selection
      # ------------------------------------------------------------------------------------------------------------------
      # Select variable according the variable list of the model to be assessed
      
      var_group_test_list = best_var[city_nb][[1]]
      
      variables= dplyr::select(cov_dataset,PfPR2_10, all_of(var_group_test_list))
      var_list = colnames(variables) ; var_list
      
      # Verify for NA atfer variable selection - and update the data
      data_updated = na.omit(cov_dataset[,c("ID",var_list)])
      coords_updated = coords[which(coords$ID %in%  data_updated$ID),c("X","Y")]
      variables = data_updated[,var_list]
      
      # define predictor list - variable list without depend variable
      predictor_dataset = dplyr::select(variables,-PfPR2_10)
      predictor_list = colnames(predictor_dataset) ; predictor_list
      
      
      RF_CovImp_pdp = RF.CovImp.pdp(variables = variables, coordinates = coords_updated, dependent_var = "PfPR2_10",
                                      var_group = var_group, city = city, mode = mode,
                                      folds_rf = folds_rf, reps_rf = reps_rf, iters_rf = iters_rf,
                                      nb_cores = my_cores, nb_test = NULL, 
                                      output_file = output_file, output_file_name = output_file_name, save_result = save_result, HYt_mode = HYt_mode, 
                                      save_RF_object = "yes")  
       
    }
  }
  
  
  
  
  # -----------------------------------------------------------------------------------------------------------------------
  #                  B. Plot results - evalution RFE 
  # -----------------------------------------------------------------------------------------------------------------------
  
  list_GOF = c("R.squared", "RMSE", "MAE")
  
  for (city_nb in 1:length(cities)){
    
    # Define city & input file path & output file name
    # ------------------------------------------------------------------------------------------------------------------
    city = cities[city_nb][[1]][1]
    city_short = cities[city_nb][[1]][2]
  
    if(var_group == "GeoSpDt"){
      
      if(mode == "ALL"){
        
        # Define output file and name
        output_file_name = paste0(cv_rf,"_",var_group,"_",mode,"_",city) ; output_file_name
        output_file = paste0(file_path,city, "/", mode,"/summary/") ; output_file
      
        # Check for file existence and create it if needed
        if(!dir.exists(file.path(output_file))){ dir.create(file.path(output_file), recursive = TRUE)   }
        
        # Create table
        # --------------------------------------------------------------------------------------------------------------------
      
        Result_tot = 0 ; Result_tot = Result_tot[-1]
        
        # For loop on the different models to compare
        # --------------------------------------------------------------------------------------------------------------------
        for (var_group_test_i in start_GeoSpDt_list:length(GeoSpDt_list)){     # var_group_test_i = 1
          
          # Define model type and list of variables
          # -------------------------------------------------------------------------------------------------------------------
          GeoSpDt_modtype = GeoSpDt_list[[var_group_test_i]]$var_group_test
          
          # Define (and create) output file path according to the model type
          # -------------------------------------------------------------------------------------------------------------------
          input_file = paste0(file_path,city, '/', mode,"/",GeoSpDt_modtype,"/") ; input_file
          
          input_file_name = list.files(input_file, pattern = "RFResults") ; input_file_name
          input_file_name = input_file_name[grep(input_file_name, pattern = paste0(city))] ; input_file_name
          input_file_name = input_file_name[grep(input_file_name, pattern = ".csv")]
          
          Result_temp = read.csv2(paste0(input_file, input_file_name))
  
          Result_temp$mode = mode
          Result_temp$GeoSpDt_modtype = GeoSpDt_modtype
  
          Result_temp$R2_mean = mean(Result_temp$R.squared)
          Result_temp$R2_sd = sd(Result_temp$R.squared)
          Result_temp$RMSE_mean = mean(Result_temp$RMSE)
          Result_temp$MSE_mean = mean(Result_temp$MSE)
          Result_temp$MAE_mean = mean(Result_temp$MAE)
          
          Result_temp$R2_median = median(Result_temp$R.squared)
          Result_temp$RMSE_median = median(Result_temp$RMSE)
          Result_temp$MSE_median = median(Result_temp$MSE)
          Result_temp$MAE_median = median(Result_temp$MAE)
          
          Result_tot = rbind(Result_tot,Result_temp)
          
        }
        
        # Summary table 
        # ------------------------------------------------------------------------------
        
        Result_tot$rfe_modtype = as.factor(Result_tot$GeoSpDt_modtype)
        rfe_modtype_list = levels(Result_tot$rfe_modtype)
       
        for (h in rfe_modtype_list){
          
          Result_tot_temp = Result_tot[which(Result_tot$rfe_modtype == h),]
          
          col_nbs = which(colnames(Result_tot_temp) %in% list_GOF)
          RF_result_summary = apply(Result_tot_temp[,col_nbs],2,summary)
          
          write.csv2(RF_result_summary, paste0(output_file,"Summary_performance_",var_group,"_",h,"_",city,".csv"))
          
        }
  
      
        # ------------------------------------------------------------------------------------------------------------------------------------
        #       Plot
        # ------------------------------------------------------------------------------------------------------------------------------------
  
        RF.plot.result(RFResult = Result_tot, factor = "GeoSpDt_modtype", stat_plot_list = c("Mean","Median"),
                        output_file = output_file, output_file_name = output_file_name)
        
      
      
      } else if(mode == "RFE"){
        
        # Create table
        Result_tot = 0 ; Result_tot = Result_tot[-1]
        
        # For loop on the different models to compare
        # --------------------------------------------------------------------------------------------------------------------
        for (var_group_test_i in start_GeoSpDt_list:length(GeoSpDt_list)){     # var_group_test_i = 1
          
          # Define model type and list of variables
          # -------------------------------------------------------------------------------------------------------------------
          GeoSpDt_modtype = GeoSpDt_list[[var_group_test_i]]$var_group_test
          
          # Define output file and name
          output_file_name = paste0(cv_rf,"_",var_group,"_",GeoSpDt_modtype,"_",mode,"_",city) ; output_file_name
          # Define output file
          output_file = paste0(file_path ,city, '/', mode,"/summary/",GeoSpDt_modtype,"/") ; output_file
        
          # Check for file existence and create it if needed
          if(!dir.exists(file.path(output_file))){ dir.create(file.path(output_file), recursive = TRUE)   }
          
          # Define (and create) output file path according to the model type
          # -------------------------------------------------------------------------------------------------------------------
          input_file = paste0(file_path,city, '/', mode,"/",GeoSpDt_modtype,"/") ; input_file
          input_file_name_list = list.files(input_file, pattern = "RFResults")
          input_file_name_list = input_file_name_list[grep(input_file_name_list, pattern = paste0(city))]
          input_file_name_list = input_file_name_list[grep(input_file_name_list, pattern = ".csv")]
          
          # Create table
          Result_tot = 0 ; Result_tot = Result_tot[-1]
          
          for (j in 1:length(input_file_name_list)){ # j = 1
            
            input_file_name = input_file_name_list[j]
            
            Result_temp = read.csv2(paste0(input_file, input_file_name))
            Result_temp$mode = mode
            Result_temp$GeoSpDt_modtype = GeoSpDt_modtype
  
            Result_temp$R2_mean = mean(Result_temp$R.squared)
            Result_temp$R2_sd = sd(Result_temp$R.squared)
            Result_temp$RMSE_mean = mean(Result_temp$RMSE)
            Result_temp$MSE_mean = mean(Result_temp$MSE)
            Result_temp$MAE_mean = mean(Result_temp$MAE)
            
            Result_temp$R2_median = median(Result_temp$R.squared)
            Result_temp$RMSE_median = median(Result_temp$RMSE)
            Result_temp$MSE_median = median(Result_temp$MSE)
            Result_temp$MAE_median = median(Result_temp$MAE)
            
            Result_tot = rbind(Result_tot,Result_temp)
          }
          
          if(is.null(Result_tot$nb_test)){
            Result_tot$nb_test =  Result_tot$test_nb
          }
          
          RF.plot.result(RFResult = Result_tot, factor = "nb_test", stat_plot_list = c("Mean","Median", "byFactor"),
                         output_file = output_file, output_file_name = output_file_name)
          
        
          # ---- Summary table ----
          Result_tot$nb_test = as.factor(Result_tot$nb_test)
          rfe_modtype_list = levels(Result_tot$nb_test)
          
          for (h in rfe_modtype_list){   #  h = "RFE-1"
            
            Result_tot_temp = Result_tot[which(Result_tot$nb_test == h),]
            
            col_nbs = which(colnames(Result_tot_temp) %in% list_GOF )
            RF_result_summary = apply(Result_tot_temp[,col_nbs],2,summary)
            RF_result_summary = as.data.frame(RF_result_summary)
            write.csv2(RF_result_summary, paste0(output_file,"Summary_performance_",output_file_name,"_",h,".csv"))
            
          }
        }
      }
    }
  }
} 
  
