# ------------------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# This function creates plots to compare the evaluation metrics of the different models

# It does for each city:
#
# 1 - Boxplots to compare different model types
# 2 - Fill the data within the summary table of RFE (with mean and median for the best fre model)
# 3 - Dependence plots
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities                 List with one vector c(city name, city short name, city epsg of projected crs) per city
#' @param input_cov              Path to input covariates
#' @param input_malaria          Path to input selected malaria data points
#' @param output_file_own        Path to input RF models
#' @param rfe_mode_plot          Either "ALL" (boxplots with models based on all covariates), "RFE" (boxplots based on RFE cov) or "ALL_RFE" (boxplots based on both)
#' @param LCZ_mode               Either "LCZonly" (uses only proportions Lcz), "LCZdponly" (uses only proportions and distances LCZ altogether, default) or "allLCZ" (uses proportions only and proportions and distances combined)

#' @return                       Nothing
#' @export                       # Per city: Boxplots (pdf) to compare different RF models
#'                                           Dependency plots separately for the two cities and combined
#--------------------------------------------------------------------------------------------------------------------------


plot.final.results.DP.BP = function(cities, input_cov, input_malaria, output_file_own, rfe_mode_plot, LCZ_mode = "LCZdponly") {

  cv_rf = "SpCVHpt" ; HYt_mode = "HYt_insideCV" ; var_group = "GeoSpDt"

  file_path = paste0(output_file_own,"/",var_group, '/')
  dir.create(file_path, showWarnings = F, recursive = T)


  # Define models to compare (var group list) GeoSpDt_list
  # ------------------------------------------------------------------------------------------------------------------

  if(LCZ_mode == "allLCZ"){

    GeoSpDt_list = list(
      # base model to compare with the other models
      list(var_group_test= "Base", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM")),

      list(var_group_test= "LCZ", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                             "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland", "LCZ_mineral",
                                             "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands")),

      list(var_group_test = "LCZ_COSMO", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                                    "COSMO_AVG_QV2M", "COSMO_AVG_RH2M", "COSMO_AVG_T2M", "COSMO_DRYAVG_TS", "COSMO_AVG_pp_avg",
                                                    "COSMO_MAX_pp_max", "COSMO_MIN_pp_min",
                                                    "TSI", "TSI_RH",
                                                    "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                    "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands")),

      list(var_group_test = "LCZ_LULC", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                                   "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                   "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                   "LC_bare_ground", "LC_building","LC_low_veg"  ,"LC_tall_veg" , "LC_water" , "LU_ACS" ,
                                                   "LU_Wetlands","LU_Planned","LU_Informal")),

      list(var_group_test = "LCZ_LULC_COSMO", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                                         "LC_bare_ground", "LC_building","LC_low_veg"  ,"LC_tall_veg" , "LC_water" , "LU_ACS" ,
                                                         "LU_Wetlands","LU_Planned","LU_Informal",
                                                         "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                         "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                         "COSMO_AVG_QV2M", "COSMO_AVG_RH2M", "COSMO_AVG_T2M", "COSMO_DRYAVG_TS",
                                                         "COSMO_AVG_pp_avg", "COSMO_MAX_pp_max", "COSMO_MIN_pp_min",  "TSI", "TSI_RH")),

      list(var_group_test= "LCZdp", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                               "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland", "LCZ_mineral",
                                               "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",

                                               "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland", "Dist_mineral",
                                               "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands")),

      list(var_group_test = "LCZdp_COSMO", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                                      "COSMO_AVG_QV2M", "COSMO_AVG_RH2M", "COSMO_AVG_T2M", "COSMO_DRYAVG_TS", "COSMO_AVG_pp_avg",
                                                      "COSMO_MAX_pp_max", "COSMO_MIN_pp_min",
                                                      "TSI", "TSI_RH",
                                                      "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                      "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                      "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland", "Dist_mineral",
                                                      "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands")),

      list(var_group_test = "LCZdp_LULC", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",

                                                     "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                     "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                     "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland", "Dist_mineral",
                                                     "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands",

                                                     "LC_bare_ground", "LC_building","LC_low_veg"  ,"LC_tall_veg" , "LC_water" , "LU_ACS" ,
                                                     "LU_Wetlands","LU_Planned","LU_Informal")),

      list(var_group_test = "LCZdp_LULC_COSMO", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",

                                                           "LC_bare_ground", "LC_building","LC_low_veg"  ,"LC_tall_veg" , "LC_water" , "LU_ACS" ,
                                                           "LU_Wetlands","LU_Planned","LU_Informal",

                                                           "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                           "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                           "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland", "Dist_mineral",
                                                           "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands",

                                                           "COSMO_AVG_QV2M", "COSMO_AVG_RH2M", "COSMO_AVG_T2M", "COSMO_DRYAVG_TS",
                                                           "COSMO_AVG_pp_avg", "COSMO_MAX_pp_max", "COSMO_MIN_pp_min",  "TSI", "TSI_RH"))
    )



  } else if(LCZ_mode == "LCZonly"){

    GeoSpDt_list = list(
      # base model to compare with the other models
      list(var_group_test= "Base", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM")),

      list(var_group_test= "LCZ", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                             "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland", "LCZ_mineral",
                                             "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands")),

      list(var_group_test = "LCZ_COSMO", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                                    "COSMO_AVG_QV2M", "COSMO_AVG_RH2M", "COSMO_AVG_T2M", "COSMO_DRYAVG_TS", "COSMO_AVG_pp_avg",
                                                    "COSMO_MAX_pp_max", "COSMO_MIN_pp_min",
                                                    "TSI", "TSI_RH",
                                                    "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                    "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands")),

      list(var_group_test = "LCZ_LULC", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                                   "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                   "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                   "LC_bare_ground", "LC_building","LC_low_veg"  ,"LC_tall_veg" , "LC_water" , "LU_ACS" ,
                                                   "LU_Wetlands","LU_Planned","LU_Informal")),

      list(var_group_test = "LCZ_LULC_COSMO", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                                         "LC_bare_ground", "LC_building","LC_low_veg"  ,"LC_tall_veg" , "LC_water" , "LU_ACS" ,
                                                         "LU_Wetlands","LU_Planned","LU_Informal",
                                                         "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                         "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                         "COSMO_AVG_QV2M", "COSMO_AVG_RH2M", "COSMO_AVG_T2M", "COSMO_DRYAVG_TS",
                                                         "COSMO_AVG_pp_avg", "COSMO_MAX_pp_max", "COSMO_MIN_pp_min",  "TSI", "TSI_RH")) )




  } else if(LCZ_mode == "LCZdponly"){

    GeoSpDt_list = list(
      # base model to compare with the other models
      list(var_group_test= "Base", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM")),

      list(var_group_test= "LCZdp", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                               "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland", "LCZ_mineral",
                                               "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",

                                               "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland", "Dist_mineral",
                                               "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands")),

      list(var_group_test = "LCZdp_COSMO", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",
                                                      "COSMO_AVG_QV2M", "COSMO_AVG_RH2M", "COSMO_AVG_T2M", "COSMO_DRYAVG_TS", "COSMO_AVG_pp_avg",
                                                      "COSMO_MAX_pp_max", "COSMO_MIN_pp_min",
                                                      "TSI", "TSI_RH",
                                                      "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                      "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                      "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland", "Dist_mineral",
                                                      "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands")),

      list(var_group_test = "LCZdp_LULC", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",

                                                     "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                     "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                     "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland", "Dist_mineral",
                                                     "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands",

                                                     "LC_bare_ground", "LC_building","LC_low_veg"  ,"LC_tall_veg" , "LC_water" , "LU_ACS" ,
                                                     "LU_Wetlands","LU_Planned","LU_Informal")),

      list(var_group_test = "LCZdp_LULC_COSMO", var_list=c("PfPR2_10", "NDVI", "NDWI", "SRTM",

                                                           "LC_bare_ground", "LC_building","LC_low_veg"  ,"LC_tall_veg" , "LC_water" , "LU_ACS" ,
                                                           "LU_Wetlands","LU_Planned","LU_Informal",

                                                           "LCZ_compact", "LCZ_indu", "LCZ_informal", "LCZ_lowland",
                                                           "LCZ_mineral", "LCZ_open", "LCZ_sparse",  "LCZ_trees", "LCZ_water", "LCZ_wetlands",
                                                           "Dist_compact", "Dist_indu", "Dist_informal", "Dist_lowland", "Dist_mineral",
                                                           "Dist_open", "Dist_sparse", "Dist_trees", "Dist_water", "Dist_wetlands",

                                                           "COSMO_AVG_QV2M", "COSMO_AVG_RH2M", "COSMO_AVG_T2M", "COSMO_DRYAVG_TS",
                                                           "COSMO_AVG_pp_avg", "COSMO_MAX_pp_max", "COSMO_MIN_pp_min",  "TSI", "TSI_RH"))
    )
  }



  # -----------------------------------------------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------------------------------------------
  #     A - Get variable datasets to get the mean
  # -----------------------------------------------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------------------------------------------

  # For loop on cities
  # -------------------------------------------------------------------------------------------------------------------
  for (city_nb in 1:length(cities)){

    # Define city & input file path & output file name
    # ------------------------------------------------------------------------------------------------------------------
    city = cities[city_nb][[1]][1]
    city_short = cities[city_nb][[1]][2]

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

    # assign cov_dataset to different variable name for each city
    assign(x = paste0("cov_dataset_",city_short), value = cov_dataset)

    # Get the mean of PfPR
    # -----------------------------------------------------------------------------------------------------------------------
    if(city_short == "DES") { mean_Pf_DES = mean(cov_dataset_DES$PfPR2_10)
    }  else if(city_short == "Kamp") { mean_Pf_KAM = mean(cov_dataset_Kamp$PfPR2_10) }

  }



  # -----------------------------------------------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------------------------------------------
  #     B - Plot results - comparison of all models - ALL and best RFE model
  # -----------------------------------------------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------------------------------------------

  # -------------------------------------------------------------------------------------------------------------------
  # Get the results data
  # -------------------------------------------------------------------------------------------------------------------

  # Define Goodness-of-fit list
  list_GOF = c("R.squared", "RMSE", "MAE")

  for (city_nb in 1:length(cities)){

    city = cities[city_nb][[1]][1]
    city_short = cities[city_nb][[1]][2]

    if(rfe_mode_plot == "RFE" | rfe_mode_plot == "ALL_RFE"){

      # Load summary table
      summary_result_rfe = read.csv2(file = paste0(file_path, city, "/RFE/summary/",var_group,"_RFE_summary.csv"))

    }

    # Create table
    # --------------------------------------------------------------------------------------------------------------------
    Result_tot = 0 ; Result_tot = Result_tot[-1]

    # For loop on the different models to compare
    # --------------------------------------------------------------------------------------------------------------------
    for (var_group_test_i in 1:length(GeoSpDt_list)){     # var_group_test_i = 1

      # Define model type and list of variables
      # -------------------------------------------------------------------------------------------------------------------
      GeoSpDt_modtype = GeoSpDt_list[[var_group_test_i]]$var_group_test

      # Define (and create) output file path according to the model type
      # -------------------------------------------------------------------------------------------------------------------

      # RFE
      # -------------------------------------------------------------------------------------------------------------------
      if(rfe_mode_plot == "RFE" | rfe_mode_plot == "ALL_RFE"){
        if(GeoSpDt_modtype != "Base"){

          input_file_RFE = paste0(file_path,city, "/RFE/",GeoSpDt_modtype,"/") ; input_file_RFE

          # pick the best model
          step1 = summary_result_rfe[which(summary_result_rfe$model == GeoSpDt_modtype),]
          step2 = step1[which(step1$..city == city | step1$..city == city_short),]

          input_file_name_RFE = list.files(input_file_RFE, pattern = "RFResults")
          input_file_name_RFE = input_file_name_RFE[grep(input_file_name_RFE, pattern = paste0(city))]
          input_file_name_RFE = input_file_name_RFE[grep(input_file_name_RFE, pattern = ".csv")]
          input_file_name_RFE = input_file_name_RFE[grep(input_file_name_RFE, pattern = paste0(step2$RFE))]

          Result_temp_RFE = read.csv2(paste0(input_file_RFE, input_file_name_RFE))


          # store into the database
          # ---------------------------------------------------------------------------------------------------------------------------------------------------------
          Result_temp_RFE$plot_name = paste0(GeoSpDt_modtype," (RFE)")
          Result_temp_RFE$mode = "RFE"

          Result_temp_RFE$R2_mean = mean(Result_temp_RFE$R.squared)
          Result_temp_RFE$R2_sd = sd(Result_temp_RFE$R.squared)
          Result_temp_RFE$RMSE_mean = mean(Result_temp_RFE$RMSE)
          Result_temp_RFE$MSE_mean = mean(Result_temp_RFE$MSE)
          Result_temp_RFE$MAE_mean = mean(Result_temp_RFE$MAE)

          Result_temp_RFE$R2_median = median(Result_temp_RFE$R.squared)
          Result_temp_RFE$RMSE_median = median(Result_temp_RFE$RMSE)
          Result_temp_RFE$MSE_median = median(Result_temp_RFE$MSE)
          Result_temp_RFE$MAE_median = median(Result_temp_RFE$MAE)

          if(city == "DES"){  meanPfPR = mean_Pf_DES }
          if(city == "Kampala"){  meanPfPR = mean_Pf_KAM }

          Result_temp_RFE$NMAE = 1-(Result_temp_RFE$MAE/meanPfPR)
          Result_temp_RFE$NMAE_median = median(Result_temp_RFE$NMAE)
          Result_temp_RFE$NRMSE = 1-((Result_temp_RFE$RMSE)^2/meanPfPR)^(1/2)
          Result_temp_RFE$NRMSE_median = median(Result_temp_RFE$NRMSE)

          # because of an error I made in the Hy tuning outside SpCV ! To be deleted if data are chanegd
          if(is.null(Result_temp_RFE$nb_test)){
            Result_temp_RFE$nb_test =  Result_temp_RFE$test_nb
            Result_temp_RFE = Result_temp_RFE[,-which(colnames(Result_temp_RFE) == "test_nb")]
          }

        }
      }


      # ALL
      # -------------------------------------------------------------------------------------------------------------------
      input_file_ALL = paste0(output_file_own,"/",var_group,'/', city, "/ALL/",GeoSpDt_modtype,"/") ; input_file_ALL

      # Select the right file
      input_file_name_ALL = list.files(input_file_ALL, pattern = "RFResults")
      input_file_name_ALL = input_file_name_ALL[grep(input_file_name_ALL, pattern = paste0(city))]
      input_file_name_ALL = input_file_name_ALL[grep(input_file_name_ALL, pattern = ".csv")]

      Result_temp_ALL = read.csv2(paste0(input_file_ALL, input_file_name_ALL))

      Result_temp_ALL$plot_name = paste0(GeoSpDt_modtype," (ALL)")
      Result_temp_ALL$nb_test = 'ALL'
      Result_temp_ALL$mode = "ALL"

      Result_temp_ALL$R2_mean = mean(Result_temp_ALL$R.squared)
      Result_temp_ALL$R2_sd = sd(Result_temp_ALL$R.squared)
      Result_temp_ALL$RMSE_mean = mean(Result_temp_ALL$RMSE)
      Result_temp_ALL$MSE_mean = mean(Result_temp_ALL$MSE)
      Result_temp_ALL$MAE_mean = mean(Result_temp_ALL$MAE)

      Result_temp_ALL$R2_median = median(Result_temp_ALL$R.squared)
      Result_temp_ALL$RMSE_median = median(Result_temp_ALL$RMSE)
      Result_temp_ALL$MSE_median = median(Result_temp_ALL$MSE)
      Result_temp_ALL$MAE_median = median(Result_temp_ALL$MAE)

      if(city_short == "DES"){  meanPfPR = mean_Pf_DES }
      if(city_short == "Kamp"){  meanPfPR = mean_Pf_KAM }

      Result_temp_ALL$NMAE = 1-(Result_temp_ALL$MAE/meanPfPR)
      Result_temp_ALL$NMAE_median = median(Result_temp_ALL$NMAE)
      Result_temp_ALL$NRMSE = 1-((Result_temp_ALL$RMSE)^2/meanPfPR)^(1/2)
      Result_temp_ALL$NRMSE_median = median(Result_temp_ALL$NRMSE)

      # select the columns
      list_column = c("R.squared", "RMSE", "MSE", "MAE", "plot_name","MAE_mean", "mode"   ,
                      "R2_mean"     ,     "R2_sd"      ,      "RMSE_mean"     ,   "MSE_mean"      ,   "MAE_mean"     ,
                      "R2_median" ,       "RMSE_median"    ,  "MSE_median" ,
                      "MAE_median"    ,   "NMAE"        ,     "NMAE_median"    , "NRMSE"         ,   "NRMSE_median" )


      col_nbs_ALL = which(colnames(Result_temp_ALL) %in% list_column )
      Result_temp_ALL = Result_temp_ALL[,col_nbs_ALL]
      if(GeoSpDt_modtype != "Base" & rfe_mode_plot != "ALL" ){
        col_nbs_RFE = which(colnames(Result_temp_RFE) %in% list_column )
        Result_temp_RFE = Result_temp_RFE[,col_nbs_RFE]

      }

      # ---------------------------
      # If RFE + ALL
      # ----------------------------
      if(rfe_mode_plot == "ALL_RFE"){

        # bind two dt
        if(GeoSpDt_modtype != "Base"){  Result_temp = rbind(Result_temp_ALL, Result_temp_RFE)      }
        if(GeoSpDt_modtype == "Base"){  Result_temp = Result_temp_ALL  }

        # bind with all
        Result_tot = rbind(Result_tot,Result_temp) ; rm(Result_temp) #; rm(Result_temp_ALL)
        if(GeoSpDt_modtype != "Base"){rm(Result_temp_RFE)}

      }

      # ---------------------------
      # If only RFE
      # ----------------------------
      if(rfe_mode_plot == "RFE"){
        if(GeoSpDt_modtype == "Base"){ Result_temp = Result_temp_ALL  }
        if(GeoSpDt_modtype != "Base"){  Result_temp = Result_temp_RFE  }

        # bind with all
        Result_tot = rbind(Result_tot,Result_temp) ; rm(Result_temp) #; rm(Result_temp_ALL)
      }

      # ---------------------------
      # If only ALL
      # ---------------------------
      if(rfe_mode_plot == "ALL"){

        Result_temp = Result_temp_ALL

        # bind with all
        Result_tot = rbind(Result_tot,Result_temp) ; rm(Result_temp) ; rm(Result_temp_ALL)
      }
    }
    assign(x = paste0("Result_tot",city_short), value = Result_tot)

  }



  # -------------------------------------------------------------------------------------------------------------------
  # Boxplot
  # -------------------------------------------------------------------------------------------------------------------

  for (city_nb in 1:length(cities)){

    # Define city & input file path & output file name
    # ------------------------------------------------------------------------------------------------------------------
    city = cities[city_nb][[1]][1]
    city_short = cities[city_nb][[1]][2]

    if(rfe_mode_plot == "RFE" | rfe_mode_plot == "ALL_RFE"){
      summary_result_rfe = read.csv2(file = paste0(file_path, "/", city, "/RFE/summary/",var_group,"_RFE_summary.csv"))
    }

    # Define output file and name
    # ----------------------------------------------------
    output_file_name = paste0(cv_rf,"_",var_group,"_",rfe_mode_plot,"_",LCZ_mode,"_",city) ; output_file_name
    output_file = paste0(output_file_own,"/",var_group, '/', city, "/RFE/summary/Boxplot/") ; output_file

    # Check for file existence and create it if needed
    if(!dir.exists(file.path(output_file))){ dir.create(file.path(output_file), recursive = TRUE)   }

    # ------------------------------------------------------------------------------------------------------------------------------------
    #       Plot
    # ------------------------------------------------------------------------------------------------------------------------------------

    if(city_short == "DES"){Result_tot = Result_totDES}
    if(city == "Kampala"){Result_tot = Result_totKamp}

    # Plot for  c("NMAE", "NRMSE")
    # ----------------------------------------
    rangeNRMSE = range(Result_totDES$NRMSE, Result_totKamp$NRMSE)
    rangeNMAE = range(Result_totDES$NMAE, Result_totKamp$NMAE)
    rangeR2 = range(Result_totDES$R.squared, Result_totKamp$R.squared)


    RF.plot.comp.result(RFResult = Result_tot, orientation = "horizontal",
                        stat_plot_list = c("Median"), GoF_list = c("NMAE", "NRMSE", "R.squared"),
                        rangeNMAE = rangeNMAE, rangeNRMSE = rangeNRMSE,rangeR2=rangeR2,
                        output_file = output_file, output_file_name = output_file_name, city = city, plot_mode = rfe_mode_plot)


  }




  # ------------------------------------------------------------------------------------------------------------------
  # PAIRWISE WILCOXSON TEST
  # ------------------------------------------------------------------------------------------------------------------

  Table_DES = Result_totDES
  Table_Kampala = Result_totKamp

  pairwise.wilcox.test(Table_Kampala$R.squared, Table_Kampala$plot_name,
                       p.adjust.method = "BH")
  pairwise.wilcox.test(Table_Kampala$RMSE, Table_Kampala$plot_name,
                       p.adjust.method = "BH")
  pairwise.wilcox.test(Table_Kampala$MAE, Table_Kampala$plot_name,
                       p.adjust.method = "BH")


  pairwise.wilcox.test(Table_DES$R.squared, Table_DES$plot_name,
                       p.adjust.method = "BH")
  pairwise.wilcox.test(Table_DES$RMSE, Table_DES$plot_name,
                       p.adjust.method = "BH")
  pairwise.wilcox.test(Table_DES$MAE, Table_DES$plot_name,
                       p.adjust.method = "BH")



  # ------------------------------------------------------------------------------------------------------------------
  # Fill the table with the list of covariates
  # ------------------------------------------------------------------------------------------------------------------


  # For loop to fill the data
  for (city_nb in 1:length(cities)){

    # Define city & input file path & output file name
    # ------------------------------------------------------------------------------------------------------------------
    city = cities[city_nb][[1]][1]
    city_short = cities[city_nb][[1]][2]

    summary_result_rfe_filled = read.csv2(file = paste0(file_path,'/', city, "/RFE/summary/",var_group,"_RFE_summary.csv"))
    summary_result_rfe_filled$RFE_cov = NA

    # create column to store mean and medium values of GoF
    for (f in 1:length(list_GOF)){ #f = 1
      summary_result_rfe_filled$temp_mean = NA
      summary_result_rfe_filled$temp_md = NA

      colnames(summary_result_rfe_filled)[which(colnames(summary_result_rfe_filled) == 'temp_mean')] = paste0("Mean_",list_GOF[f])
      colnames(summary_result_rfe_filled)[which(colnames(summary_result_rfe_filled) == 'temp_md')] = paste0("Median_",list_GOF[f])
    }

    # Define output file and name
    output_file_name = paste0(cv_rf,"_",var_group,"_ALL-RFE_",city) ; output_file_name
    output_file = paste0(file_path,city, "/RFE/summary/Boxplot/") ; output_file


    # For loop on the different models to compare
    # --------------------------------------------------------------------------------------------------------------------
    for (var_group_test_i in 2:length(GeoSpDt_list)){     # var_group_test_i = 2

      # Define model type and list of variables
      # -------------------------------------------------------------------------------------------------------------------
      GeoSpDt_modtype = GeoSpDt_list[[var_group_test_i]]$var_group_test

      # Define (and create) output file path according to the model type
      # -------------------------------------------------------------------------------------------------------------------

      # RFE
      # -------------------------------------------------------------------------------------------------------------------

      # fill cov

      input_file_RFE = paste0(file_path,city, "/RFE/",GeoSpDt_modtype,"/") ; input_file_RFE
      step1 = summary_result_rfe_filled[which(summary_result_rfe_filled$model == GeoSpDt_modtype),]
      step2 = step1[which(step1$..city == city | step1$..city == city_short),]

      input_file_name_RFE = list.files(input_file_RFE, pattern = "RFResults")
      input_file_name_RFE = input_file_name_RFE[grep(input_file_name_RFE, pattern = paste0(city))]
      input_file_name_RFE = input_file_name_RFE[grep(input_file_name_RFE, pattern = ".csv")]
      input_file_name_RFE = input_file_name_RFE[grep(input_file_name_RFE, pattern = paste0(step2$RFE))]

      Result_temp_RFE = read.csv2(paste0(input_file_RFE, input_file_name_RFE))

      summary_result_rfe_filled$RFE_cov[which(summary_result_rfe_filled$model == GeoSpDt_modtype & summary_result_rfe_filled$city == city)] = paste(Result_temp_RFE$covariates[1])

      # fill performance

      inputfile_perf_RFE = paste0(file_path, city, "/RFE/summary/",GeoSpDt_modtype,"/")
      list_files = list.files(inputfile_perf_RFE)
      list_files = list_files[grep(list_files, pattern = ".csv")]
      list_files = list_files[grep(list_files, pattern = paste0(step2$RFE))]
      # get the result
      perf_RFE = read.csv2(paste0(inputfile_perf_RFE, list_files))
      perf_RFE$MSE[which(perf_RFE$X == 'Mean')]

      summary_result_rfe_filled$Mean_R.squared[which(summary_result_rfe_filled$model == GeoSpDt_modtype & summary_result_rfe_filled$city == city)] = perf_RFE$R.squared[which(perf_RFE$X == 'Mean')]
      summary_result_rfe_filled$Mean_RMSE[which(summary_result_rfe_filled$model == GeoSpDt_modtype & summary_result_rfe_filled$city == city)] = perf_RFE$RMSE[which(perf_RFE$X == 'Mean')]
      summary_result_rfe_filled$Mean_MAE[which(summary_result_rfe_filled$model == GeoSpDt_modtype & summary_result_rfe_filled$city == city)] = perf_RFE$MAE[which(perf_RFE$X == 'Mean')]
      summary_result_rfe_filled$Median_R.squared[which(summary_result_rfe_filled$model == GeoSpDt_modtype & summary_result_rfe_filled$city == city)] = perf_RFE$R.squared[which(perf_RFE$X == 'Median')]
      summary_result_rfe_filled$Median_RMSE[which(summary_result_rfe_filled$model == GeoSpDt_modtype & summary_result_rfe_filled$city == city)] = perf_RFE$RMSE[which(perf_RFE$X == 'Median')]
      summary_result_rfe_filled$Median_MAE[which(summary_result_rfe_filled$model == GeoSpDt_modtype & summary_result_rfe_filled$city == city)] = perf_RFE$MAE[which(perf_RFE$X == 'Median')]

    }
  }

  write.csv2(summary_result_rfe_filled, file = paste0(file_path,city, "/RFE/summary/",var_group,"_RFE_summary_covFilled.csv"))


  # ---------------------------------------------------------------------------------------
  # Dependence plot for the two countries
  # ---------------------------------------------------------------------------------------

  #best set of cov
  GeoSpDt_modtype = "LCZdp_LULC_COSMO"

  summary_result_rfe = data.frame()

  for (city_nb in 1:length(cities)){

    city = cities[city_nb][[1]][1]
    city_short = cities[city_nb][[1]][2]
    sum_result_rfe = read.csv2(file = paste0(file_path,city, "/RFE/summary/",var_group,"_RFE_summary.csv"))

    summary_result_rfe = rbind(summary_result_rfe, sum_result_rfe)

  }


  step1 = summary_result_rfe[which(summary_result_rfe$model == GeoSpDt_modtype),]

  model_DES = as.character(step1[which(step1$..city == "DES"),]$X10var_Model)
  model_KAM = as.character(step1[which(step1$..city == "Kampala"),]$X10var_Model)

  # Define output file and name
  output_file_name = paste0(cv_rf,"_",var_group,"_",GeoSpDt_modtype,"_RFE-best") ; output_file_name
  output_file = paste0(file_path, "DP_common/") ; output_file

  # Check for file existence and create it if needed
  if(!dir.exists(file.path(output_file))){ dir.create(file.path(output_file), recursive = TRUE)   }

  DepPlot_KAM = readRDS(paste0(file_path, "Kampala", "/RFE/",GeoSpDt_modtype,"/", "DepPlotValues_",cv_rf,"_",var_group,"_",model_KAM[1],"_Kampala.rds"))
  DepPlot_DES = readRDS(paste0(file_path, "Dar es salaam", "/RFE/",GeoSpDt_modtype,"/", "DepPlotValues_",cv_rf,"_",var_group,"_",model_DES[1],"_Dar es salaam.rds"))
  names(DepPlot_KAM)[which(names(DepPlot_KAM) == "COSMO_2014.06.09_AVG_TSI")] = "AVG_TSI" ; names(DepPlot_KAM)[which(names(DepPlot_KAM) == "COSMO_2014.06.09_AVG_AVG_RH2M")] = "AVG_RH2M"
  names(DepPlot_KAM)[which(names(DepPlot_KAM) == "COSMO_2014.06.09_AVG_QV2M")] = "AVG_QV2M" ; names(DepPlot_KAM)[which(names(DepPlot_KAM) == "COSMO_2014.06.09_AVG_AVG_T2M")] = "AVG_T2M"
  names(DepPlot_KAM)[which(names(DepPlot_KAM) == "LU_Informal")] = "LU_informal"
  names(DepPlot_DES)[which(names(DepPlot_DES) == "TSI")] = "AVG_TSI" ; names(DepPlot_DES)[which(names(DepPlot_DES) == "LU_Informal")] = "LU_informal"

  predictor_list_common = intersect(names(DepPlot_DES), names(DepPlot_KAM)) ; predictor_list_common

  for (var_nb in 1:length(predictor_list_common)){
    var_name = predictor_list_common[var_nb]
    var_DpPlot_KAMdataset = DepPlot_KAM[[var_name]]
    var_DpPlot_DESdataset = DepPlot_DES[[var_name]]

    nb_model = length(var_DpPlot_KAMdataset)-1

    # Plot
    pdf(paste0(output_file,"DdtPlot_",output_file_name,"_",var_name,"_bothCities.pdf"))

    par(mar=c(5,5,3,2))

    range_yValues = range(var_DpPlot_KAMdataset[,2:nb_model], var_DpPlot_DESdataset[,2:nb_model])
    range_xValues = range(var_DpPlot_KAMdataset[,1], var_DpPlot_DESdataset[,1])

    plot(y = "", x = "", xlab = "", ylab="", xlim=range_xValues, ylim=range_yValues,
         main = paste(var_name), cex.axis = 2, cex.main = 2)

    # create rug
    dec_KAM = quantile(na.omit(cov_dataset_Kamp[,var_name]), prob = seq(0, 1, length = 11), type = 5)
    dec_KAM = dec_KAM[-c(1,length(dec_KAM))]
    rug(x=dec_KAM, col=rgb(0.7,0.5,0,0, alpha = 0.6), lwd="1.5", side = 3)

    dec_DES = quantile(na.omit(cov_dataset_DES[,var_name]), prob = seq(0, 1, length = 11), type = 5)
    dec_DES = dec_DES[-c(1,length(dec_DES))]
    rug(x=dec_DES, col=rgb(0,0.5,0.7,0, alpha = 0.6), lwd="1.5")

    #p = 1
    for (p in 1:nb_model){
      lines(x=var_DpPlot_KAMdataset[,1], y=var_DpPlot_KAMdataset[,(p+1)], col=rgb(0.7,0.5,0,0,alpha = 0.2), lwd = 1.2)
      lines(x=var_DpPlot_DESdataset[,1], y=var_DpPlot_DESdataset[,(p+1)], col=rgb(0,0.5,0.7,0,alpha = 0.2), lwd = 1.2)
    }

    summary_dpt_KAM = apply(var_DpPlot_KAMdataset[,2:nb_model],1,function(x) { cbind(mean(x), sd(x)) } )
    summary_dpt_DES = apply(var_DpPlot_DESdataset[,2:nb_model],1,function(x) { cbind(mean(x), sd(x)) } )

    lines(x = var_DpPlot_KAMdataset[,1], y = summary_dpt_KAM[1,], col=rgb(0.7,0.5,0), lwd = 2)
    lines(x = var_DpPlot_KAMdataset[,1], y = (summary_dpt_KAM[1,]+summary_dpt_KAM[2,]) , lty = 2, col=rgb(0.7,0.5,0), lwd = 2)
    lines(x = var_DpPlot_KAMdataset[,1], y = (summary_dpt_KAM[1,]-summary_dpt_KAM[2,]) , lty = 2, col=rgb(0.7,0.5,0), lwd = 2)

    lines(x = var_DpPlot_DESdataset[,1], y = summary_dpt_DES[1,], col=rgb(0,0.5,0.7), lwd = 2)
    lines(x = var_DpPlot_DESdataset[,1], y = (summary_dpt_DES[1,]+summary_dpt_DES[2,]) , lty = 2, col=rgb(0,0.5,0.7), lwd = 2)
    lines(x = var_DpPlot_DESdataset[,1], y = (summary_dpt_DES[1,]-summary_dpt_DES[2,]) , lty = 2, col=rgb(0,0.5,0.7), lwd = 2)

    dev.off()

  }

  pdf(paste0(output_file,"ALL_DdtPlot_",output_file_name,"_RF.pdf"),
      width = 7)

  par(mar=c(5,4,3,2), mfrow = c(2,2))

  for (var_nb in 1:length(predictor_list_common)){

    var_name = predictor_list_common[var_nb]

    var_DpPlot_KAMdataset = DepPlot_KAM[[var_name]]
    var_DpPlot_DESdataset = DepPlot_DES[[var_name]]

    range_yValues = range(var_DpPlot_KAMdataset[,2:nb_model], var_DpPlot_DESdataset[,2:nb_model])
    range_xValues = range(var_DpPlot_KAMdataset[,1], var_DpPlot_DESdataset[,1])

    plot(y = "", x = "", xlab = "", ylab="", xlim=range_xValues, ylim=range_yValues,
         main = paste(var_name))

    # create rug
    dec_KAM = quantile(na.omit(cov_dataset_Kamp[,var_name]), prob = seq(0, 1, length = 11), type = 5)
    dec_KAM = dec_KAM[-c(1,length(dec_KAM))]
    rug(x=dec_KAM, col=rgb(0.7,0.5,0,0, alpha = 0.6), lwd="1", side = 3)

    dec_DES = quantile(na.omit(cov_dataset_DES[,var_name]), prob = seq(0, 1, length = 11), type = 5)
    dec_DES = dec_DES[-c(1,length(dec_DES))]
    rug(x=dec_DES, col=rgb(0,0.5,0.7,0, alpha = 0.6), lwd="1")

    #p = 1
    for (p in 1:nb_model){
      lines(x=var_DpPlot_KAMdataset[,1], y=var_DpPlot_KAMdataset[,(p+1)], col=rgb(0.7,0.5,0,0,alpha = 0.2))
      lines(x=var_DpPlot_DESdataset[,1], y=var_DpPlot_DESdataset[,(p+1)], col=rgb(0,0.5,0.7,0,alpha = 0.2))
    }

    summary_dpt_KAM = apply(var_DpPlot_KAMdataset[,2:nb_model],1,function(x) { cbind(mean(x), sd(x)) } )
    summary_dpt_DES = apply(var_DpPlot_DESdataset[,2:nb_model],1,function(x) { cbind(mean(x), sd(x)) } )

    lines(x = var_DpPlot_KAMdataset[,1], y = summary_dpt_KAM[1,], col=rgb(0.7,0.5,0))
    lines(x = var_DpPlot_KAMdataset[,1], y = (summary_dpt_KAM[1,]+summary_dpt_KAM[2,]) , lty = 2, col=rgb(0.7,0.5,0))
    lines(x = var_DpPlot_KAMdataset[,1], y = (summary_dpt_KAM[1,]-summary_dpt_KAM[2,]) , lty = 2, col=rgb(0.7,0.5,0))
    lines(x = var_DpPlot_DESdataset[,1], y = summary_dpt_DES[1,], col=rgb(0,0.5,0.7))
    lines(x = var_DpPlot_DESdataset[,1], y = (summary_dpt_DES[1,]+summary_dpt_DES[2,]) , lty = 2, col=rgb(0,0.5,0.7))
    lines(x = var_DpPlot_DESdataset[,1], y = (summary_dpt_DES[1,]-summary_dpt_DES[2,]) , lty = 2, col=rgb(0,0.5,0.7))

  }
  dev.off()
  par(mar=c(5,4,3,2), mfrow = c(1,1))



  # ---------------------------------------------------------------------------------------
  # Dependence plot for the two countries SEPARETELy
  # ---------------------------------------------------------------------------------------

  for (city_nb in 1:length(cities)){

    city = cities[city_nb][[1]][1]
    city_short = cities[city_nb][[1]][2]

    output_file = paste0(file_path, city, "/RFE/summary/DP/") ; output_file
    if(!dir.exists(file.path(output_file))){ dir.create(file.path(output_file), recursive = TRUE)   }

    if (city == "Kampala"){
      DepPlot = readRDS(paste0(file_path, city, "/RFE/",GeoSpDt_modtype,"/", "DepPlotValues_",cv_rf,"_",var_group,"_",model_KAM[1],"_Kampala.rds"))
      names(DepPlot)[which(names(DepPlot) == "COSMO_2014.06.09_AVG_TSI")] = "AVG_TSI" ; names(DepPlot)[which(names(DepPlot) == "COSMO_2014.06.09_AVG_RH2M")] = "AVG_RH2M"
      names(DepPlot)[which(names(DepPlot) == "COSMO_2014.06.09_AVG_QV2M")] = "AVG_QV2M" ; names(DepPlot)[which(names(DepPlot) == "COSMO_2014.06.09_AVG_T2M")] = "AVG_T2M" ; names(DepPlot)[which(names(DepPlot) == "COSMO_2014.06.09_AVG_TS")] = "AVG_TS"
      names(DepPlot)[which(names(DepPlot) == "LU_Informal")] = "LU_informal" ;   names(DepPlot)[which(names(DepPlot) == "LU_Planned")] = "LU_planned"

    } else {
      DepPlot = readRDS(paste0(file_path, "Dar es salaam", "/RFE/",GeoSpDt_modtype,"/", "DepPlotValues_",cv_rf,"_",var_group,"_",model_DES[1],"_Dar es salaam.rds"))
      names(DepPlot)[which(names(DepPlot) == "TSI")] = "AVG_TSI" ; names(DepPlot)[which(names(DepPlot) == "LU_Informal")] = "LU_informal"
      names(DepPlot)[which(names(DepPlot) == "SRTM")] = "SRTM_30" ; names(DepPlot)[which(names(DepPlot) == "LC_tall_veg")] = "LC_trees"
    }

    predictor_list_common = names(DepPlot)

    for (var_nb in 1:length(predictor_list_common)){

      var_name = predictor_list_common[var_nb]
      var_DpPlot = DepPlot[[var_name]]
      nb_model = length(var_DpPlot)-1

      # Plot
      pdf(paste0(output_file,"DdtPlot_",output_file_name,"_",var_name,"_",city,".pdf"))
      par(mar=c(5,5,3,2))

      range_yValues = range(var_DpPlot[,2:nb_model])
      range_xValues = range(var_DpPlot[,1])

      plot(y = "", x = "", xlab = "", ylab="", xlim=range_xValues, ylim=range_yValues,
           main = paste(var_name), cex.axis = 2, cex.main = 2)

      # create rug
      if(city == 'Kampala'){
        dec_KAM = quantile(na.omit(cov_dataset_Kamp[,var_name]), prob = seq(0, 1, length = 11), type = 5)
        dec_KAM = dec_KAM[-c(1,length(dec_KAM))]
        rug(x=dec_KAM, col=rgb(0.7,0.5,0,0, alpha = 0.6), lwd="1.5", side = 3)
      }
      if(city_short == 'DES'){
        dec_DES = quantile(na.omit(cov_dataset_DES[,var_name]), prob = seq(0, 1, length = 11), type = 5)
        dec_DES = dec_DES[-c(1,length(dec_DES))]
        rug(x=dec_DES, col=rgb(0,0.5,0.7,0, alpha = 0.6), lwd="1.5")
      }
      #p = 1
      for (p in 1:nb_model){

        lines(x=var_DpPlot[,1], y=var_DpPlot[,(p+1)], col=rgb(0.7,0.5,0,0,alpha = 0.2), lwd = 1.2)

      }

      summary_dpt = apply(var_DpPlot[,2:nb_model],1,function(x) { cbind(mean(x), sd(x)) } )

      if(city == 'Kampala'){
        lines(x = var_DpPlot[,1], y = summary_dpt[1,], col=rgb(0.7,0.5,0), lwd = 2)
        lines(x = var_DpPlot[,1], y = (summary_dpt[1,]+summary_dpt[2,]) , lty = 2, col=rgb(0.7,0.5,0), lwd = 2)
        lines(x = var_DpPlot[,1], y = (summary_dpt[1,]-summary_dpt[2,]) , lty = 2, col=rgb(0.7,0.5,0), lwd = 2)
      }

      if(city_short == 'DES'){
        lines(x = var_DpPlot[,1], y = summary_dpt[1,], col=rgb(0,0.5,0.7), lwd = 2)
        lines(x = var_DpPlot[,1], y = (summary_dpt[1,]+summary_dpt[2,]) , lty = 2, col=rgb(0,0.5,0.7), lwd = 2)
        lines(x = var_DpPlot[,1], y = (summary_dpt[1,]-summary_dpt[2,]) , lty = 2, col=rgb(0,0.5,0.7), lwd = 2)
      }

      dev.off()

    }
  }
}
