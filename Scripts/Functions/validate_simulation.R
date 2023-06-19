# ------------------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# This code validates the simulation study: predictive performance is computed on the testing non-displaced sets
#
# It create for each city:
# - a table (csv file) with the predictive performance of all methods (m0, m1, m2 TM)
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities                 List with one vector c(city name, city short name, city epsg of projected crs) per city
#' @param Dir_input_cov          Path to input covariate tables
#' @param Dir_input_malaria_data Path to input selected malaria data points
#' @param Dir_output_rf          Path to RF models

#' @return                       Nothing
#' @export                       # Per city: a table with results of the simulation study
#--------------------------------------------------------------------------------------------------------------------------


validate.simulation = function(cities, Dir_input_cov, Dir_input_malaria_data, Dir_output_rf){

    subdir = "simul_displ/"

    for (city_vect in cities){

        city = city_vect[1] ; city
        city_short = city_vect[2]
        city_i = match(list(city_vect), cities)

        all_metrics = data.frame()

        for (method in c( "m0", "TM", "m1", "m2")){

            # Read RFE results
            Dir_final_results = paste0(Dir_output_rf, city,  "/", subdir)
            rfe_results = read.csv2(paste0(Dir_final_results, method, "/RFE/", "RFE_metrics_", city_short, ".csv"))

            best_iter = rfe_results %>%
                slice_min(mean_RMSE)

            # Get best rfe iteration
            best_iter_df = read.csv2(paste0(Dir_final_results, method, "/RFE/", "RFE_", best_iter$nb_test, "/RF_Results_", city_short, ".csv"))

            # Get best covariates
            best_var = best_iter_df$covariates[1] %>%
                paste0(",ID") %>%
                strsplit(",")

            best_var = best_var[[1]]


            # Load data
            # ---------------------------------------------------------------
            # upload malaria data points
            test_data = paste0(Dir_input_malaria_data, city, "/", subdir, '/malaria_data_', city_short, '_test_set.csv') %>%
                read.csv2() %>%
                #dplyr::filter(UP_AGE <= 16) %>%
                mutate(ID = as.numeric(ID)) %>%
                dplyr::select(ID, PfPR2_10)

            newdata = paste0(Dir_input_cov, city, '/', subdir,  city_short, "_full_cov.csv") %>%
                read.csv2(sep =";") %>%
                mutate(ID = as.numeric(ID)) %>%
                dplyr::select(best_var) %>%
                dplyr::right_join(test_data, by=c("ID"= "ID")) %>%
                na.omit()

            # Coordinates of new data
            coords = paste0(Dir_input_malaria_data, city, '/malaria_data_', city_short, '.shp') %>%
                st_read() %>%
                mutate(ID = as.numeric(ID)) %>%
                right_join(newdata) %>%
                dplyr::select(ID, geometry)

            # upload rf results
            RF_CovImp_pdp = paste0(Dir_final_results, method, "/RFE/", "RFE_", best_iter$nb_test, "/RF_output_", city_short, ".rds") %>%
                readRDS()

            # Prediction
            # ---------------------------------------------------------------
            pred_rf_tot = lapply(RF_CovImp_pdp$models, function(x){predict(x, newdata=newdata)}$data)
            pred_rf_tot_df = do.call(cbind, pred_rf_tot) %>%
                as.data.frame()

            # Avg by row
            pred_avg_df = rowMeans(pred_rf_tot_df)
            pred_rf_cat = cbind("cat" = newdata$ID, "response" = pred_avg_df, "obs" = newdata$PfPR2_10) %>%
                as.data.frame() %>%
                mutate(residuals = abs(obs - response))

            # Metrics
            # ---------------------------------------------------------------
            metrics_df = data.frame(m = method)
            metrics_df = metrics_df %>%
                mutate(mse = mean(pred_rf_cat$residuals^2)) %>%
                mutate(rmse = sqrt(mse)) %>%
                mutate(mae = mean(pred_rf_cat$residuals)) %>%
                mutate(r2 = 1 - mse/var(pred_rf_cat$obs))

            all_metrics = rbind(metrics_df, all_metrics)

        }

        write.csv2(all_metrics, paste0(Dir_final_results, "validation_test.csv"))

    }
}


