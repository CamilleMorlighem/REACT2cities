# ------------------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# This code uses a random forest regressor to model malaria risk, with the following :
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
# It create for each city:
# - "RFCovImp... .csv" : a table (csv file) with the RF results (Goodness-of-fit indices, RF parameters, etc.),
#    each line containes the result for one simulation of the spatial cross-validation
# - a table (csv file) with the importance for each covariates for each simulation coming from the spatial cross-validation
# - a barplot (pdf file) with the average of covariates importance (across all simulations)
# - partial dependence plots (pdf files) for each covariate included in the model (each line grey results form a simulation)
# - rds files with the partial dependence values
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities                 List with one vector c(city name, city short name, city epsg of projected crs) per city
#' @param n.cores                Number of cores used for parallelelization
#' @param reps_rf                Number of repetitions in the spatial cross-validation (default 10)
#' @param folds_rf               Number of folds used in the spatial cross-validation (default 5)
#' @param iters_rf               Number of sub-folds used for hyperparameters tuning inside the spatial cross-validation (default 2)
#' @param save_result            "yes" or "no" (default "yes")
#' @param Dir_input_cov          Path to input covariate tables
#' @param Dir_input_malaria_data Path to input selected malaria data points
#' @param Dir_output_rf          Path to output RF models
#' @param method                 Method to implement, either "m0", "m1", "m2" or "TM" (useful for simulation study) (Default "m0")
#' @param simulation             Boolean indicates whether simulation study is conducted (Default F)

#' @return                       Nothing
#' @export                       # Per city: "RFCovImp... .csv" : a table with RF results
#                                             A csv with covariate importance for each simulation coming from the spatial cross-validation
#                                             A barplot (pdf) with the average covariates importance (across all simulations)
#                                             Partial dependence plots (pdf) for each covariate included in the model
#                                             rds files with the partial dependence values
#--------------------------------------------------------------------------------------------------------------------------


rf.modelling = function(cities, n.cores, reps_rf = 5, folds_rf = 5, iters_rf = 2, save_result = 'yes', Dir_input_cov, Dir_input_malaria_data, Dir_output_rf, method = "m0", simulation = F){

    for (city_vect in cities){

        city = city_vect[1]
        city_short = city_vect[2]
        myepsg = city_vect[3]

        subdir = '/'
        save_suffix = ''

        if (simulation) {

            subdir = 'simul_displ/'

            if (method == "m0") { #use displaced data
                save_suffix = "_displ"
            }
        }

        if (method == "m1" | method == "m2"){

            dupl_subdir = 'duplicates/'
            save_suffix = ''
            joinby = c("Orig_ID" = "ID")

        } else {

            dupl_subdir = ''
            joinby = c("ID" = "ID")

        }

        # 1. Load malaria and covariates data
        #-----------------------------------------------

        # Load shapefile of points with malaria prevalence, as 'sf' object
        malaria_data = paste0(Dir_input_malaria_data, city, '/', subdir, dupl_subdir, "malaria_data_", city_short, save_suffix, '.csv') %>%
            read.csv2()  %>%
            mutate(ID = as.numeric(ID),
                   PfPR2_10 = as.numeric(PfPR2_10)) %>%
            dplyr::select(-geometry)

        if (simulation) {
            #Load training set
            training_set = paste0(Dir_input_malaria_data, city, '/', subdir, "malaria_data_", city_short, '_training_set.csv') %>%
                read.csv2() %>%
                dplyr::select(ID)

            malaria_data = malaria_data %>%
                right_join(training_set, by = joinby)

        }

        malaria_sf = paste0(Dir_input_malaria_data, city, '/', subdir, dupl_subdir, "malaria_data_", city_short, save_suffix, '.shp') %>%
            read_sf() %>%
            mutate(ID = as.numeric(ID)) %>%
            dplyr::select(ID) %>%
            right_join(malaria_data) %>%
            dplyr::select(ID, PfPR2_10) ; nrow(malaria_sf)

        # store the coordinates of malaria data
        coordinates_malaria = malaria_sf %>%
            st_coordinates() %>%
            cbind(ID = malaria_sf$ID) %>%
            as.data.frame()

        # Load covariates csv into a dataframe
        cov_df = paste0(Dir_input_cov, city, '/', subdir, dupl_subdir, city_short, "_full_cov", save_suffix, ".csv") %>%
            read.csv2(sep =";") %>%
            dplyr::select(-matches("dist")) %>%
            right_join(st_drop_geometry(malaria_sf), by=c("ID"= "ID")) %>%
            na.omit() ;

        # After Verifying for NA value, update the data
        coords = coordinates_malaria[which(coordinates_malaria$ID %in% cov_df$ID), c("X","Y")]
        cov_dataset = dplyr::select(cov_df, -ID)

        # define predictor list and dataset - variable without dependent variable
        predictor_dataset = dplyr::select(cov_dataset,-PfPR2_10)
        predictor_list = colnames(predictor_dataset)

        # 2. RF modelling
        #-----------------------------------------------

        if (method == "m0" | method == "m1" | method == "TM") {

            Dir_output_rf_rfe = paste0(Dir_output_rf, city, "/", subdir, method, "/RFE/")
            dir.create(Dir_output_rf_rfe, showWarnings = F, recursive = T)

            metrics_df = data.frame()

            for (nb_test in 1:(length(predictor_list)-1)){

                CovImp_summary = RF.CovImp.pdp(variables = cov_dataset, coordinates = coords,
                                               city = city, mode = "RFE", nb_test = nb_test,
                                               city_short = city_short, Dir_output_rf = Dir_output_rf_rfe,
                                               folds_rf = folds_rf, reps_rf = reps_rf, iters_rf = iters_rf,
                                               nb_cores = n.cores,
                                               save_RF_object = save_result,  save_result = save_result)

                # Delete the less important covariate
                Cov_minImp = CovImp_summary$var[which(CovImp_summary$mean == min(CovImp_summary$mean))]
                cov_dataset = dplyr::select(cov_dataset, -Cov_minImp)
                var_list = colnames(cov_dataset) ; var_list

                # Get mean R2
                df = paste0(Dir_output_rf_rfe, 'RFE_', nb_test, '/', "RF_Results_", city_short, ".csv") %>%
                    read.csv2(sep = ";") %>%
                    get_eval_metrics(nb_test)

                metrics_df = rbind(metrics_df, df)
            }

            write.csv2(metrics_df, file = paste0(Dir_output_rf_rfe,"/RFE_metrics_", city_short,".csv"))

        } else if (method == "m2"){

            # add orig ID attribute
            cov_df = malaria_data %>%
                dplyr::select(ID, Orig_ID) %>%
                inner_join(cov_df) ; nrow(cov_df)

            Dir_output_rf_all_base = paste0(Dir_output_rf, city, "/", subdir, method, "/", "ALL", save_suffix, "/")
            dir.create(Dir_output_rf_all_base, recursive = T, showWarnings = F)

            metrics_df = data.frame()

            # repeats over 1000 iterations
            for (nb_test in 1:1000){

                Dir_output_rf_all = paste0(Dir_output_rf_all_base, "All_", nb_test, "/")
                dir.create(Dir_output_rf_all, showWarnings = F, recursive = T)

                # group df by original ID (group all 9 duplicated points) and select randomly one point per group
                cov_df_one_per_group <- cov_df %>%
                    group_by(Orig_ID) %>%
                    sample_n(1) %>%
                    data.frame()

                cov_dataset = cov_df_one_per_group %>%
                    dplyr::select(-ID, -Orig_ID)

                # extract coordinates
                coordinates = inner_join(cov_df_one_per_group, malaria_sf)
                coords = coordinates$geometry %>%
                    st_coordinates() %>%
                    data.frame()

                #save points selection to csv
                write.csv2(coordinates, file = paste0(Dir_output_rf_all,"/malaria_data_pts_selected_", city_short,".csv"))

                RF_CovImp_pdp = RF.CovImp.pdp(variables = cov_dataset, coordinates = coords,
                                              city = city, city_short = city_short, mode = "ALL", nb_test = nb_test,
                                              folds_rf = folds_rf, reps_rf = reps_rf, iters_rf = iters_rf,
                                              nb_cores = n.cores, Dir_output_rf = Dir_output_rf_all,
                                              save_RF_object = save_result,  save_result = save_result)

                # Get mean R2
                df = paste0(Dir_output_rf_all, "RF_Results_", city_short, ".csv") %>%
                    read.csv2(sep = ";") %>%
                    get_eval_metrics(nb_test)

                metrics_df = rbind(metrics_df, df)
            }

            write.csv2(metrics_df, file = paste0(Dir_output_rf_all_base,"/method_2_metrics_", city_short,".csv"))

            # find best iteration
            results_iteration = paste0(Dir_output_rf_all_base,"/method_2_metrics_", city_short,".csv") %>%
                read.csv2()

            best_iter = results_iteration %>%
                slice_min(mean_RMSE) %>%
                dplyr::select(nb_test)

            # extract displaced data of this iteration
            cov_df = paste0(Dir_output_rf_all_base, "/All_", best_iter, "/malaria_data_pts_selected_", city_short, ".csv") %>%
                read.csv2(sep = ";") %>%
                dplyr::select(-c(X, geometry, Orig_ID)) %>%
                mutate(ID = as.numeric(ID)) %>%
                na.omit()

            # After Verifying for NA value, update the data
            coords = coordinates_malaria[which(coordinates_malaria$ID %in% cov_df$ID), c("X","Y")]
            cov_dataset = dplyr::select(cov_df, -ID)

            # define predictor list and dataset - variable without dependent variable
            predictor_dataset = dplyr::select(cov_dataset,-PfPR2_10)
            predictor_list = colnames(predictor_dataset)

            Dir_output_rf_rfe = paste0(Dir_output_rf, city, "/", subdir, method, "/RFE", save_suffix, "/")
            dir.create(Dir_output_rf_rfe, showWarnings = F, recursive = T)

            metrics_df = data.frame()

            # performs RFE
            for (nb_test in 1:(length(predictor_list)-1)){

                CovImp_summary = RF.CovImp.pdp(variables = cov_dataset, coordinates = coords,
                                               city = city, city_short = city_short, mode = "RFE", nb_test = nb_test,
                                               folds_rf = folds_rf, reps_rf = reps_rf, iters_rf = iters_rf,
                                               nb_cores = n.cores, Dir_output_rf = Dir_output_rf_rfe,
                                               save_RF_object = save_result,  save_result = save_result)

                # Delete the less important covariate
                Cov_minImp = CovImp_summary$var[which(CovImp_summary$mean == min(CovImp_summary$mean))]
                cov_dataset = dplyr::select(cov_dataset, -Cov_minImp)
                var_list = colnames(cov_dataset) ; var_list

                # Get mean R2
                df = paste0(Dir_output_rf_rfe, 'RFE_', nb_test, '/', "RF_Results_", city_short, ".csv") %>%
                    read.csv2(sep = ";") %>%
                    get_eval_metrics(nb_test)

                metrics_df = rbind(metrics_df, df)
            }

            write.csv2(metrics_df, file = paste0(Dir_output_rf_rfe,"/RFE_metrics_", city_short,".csv"))
        }
    }



}
