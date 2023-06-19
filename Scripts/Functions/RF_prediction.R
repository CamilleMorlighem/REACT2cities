# ------------------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# This code uses a random forest regressor to predict malaria risk, with the following :
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
# This function predicts pfpr2-10 using the best set of covariates for each city and plots the predictions
#
# It create for each city:
# - a predictive map stored as a raster
# - a png plot of the predictive map
# - 6 different plots of the predictions (some as the raw predictive maps, others where the predictions are aggregated per admin area)
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities                 List with one vector c(city name, city short name, city epsg of projected crs) per city
#' @param crs_REACT              Geographic coordinate reference system
#' @param Dir_input_cov          Path to input covariate tables
#' @param Dir_input_malaria_data Path to input selected malaria data points
#' @param Dir_output_rf          Path to input RF models
#' @param Dir_output_pred        Path to predicted outputs
#' @param Dir_input_pts          Path to output RF models
#' @param Dir_input_admin        Path to adm boundaries shp
#' @param method                 Method to implement, either "m0", "m1", "m2" or "TM" (useful for simulation study) (Default "m0")
#' @param simulation             Boolean indicates whether simulation study is conducted (Default F)

#' @return                       Nothing
#' @export                       # Per city: A predictive map stored as a raster
#                                            A png plot of the predictive map
#                                            6 different plots of the predictions (some as the raw predictive maps, others where the predictions are aggregated per admin area)
#--------------------------------------------------------------------------------------------------------------------------



rf.prediction = function(cities, crs_REACT, Dir_input_cov, Dir_input_malaria_data, Dir_output_rf, Dir_output_pred, Dir_input_pts, Dir_input_admin, method = "m0", simulation = F){

  subdir = '/'
  save_suffix = ''

  if (simulation) {

    subdir = 'simul_displ/'

    if (method != "TM") { #use displaced data
      save_suffix = "_displ"
    }
  }

  if (method == "m1" | method == "m2"){

    dupl_subdir = 'duplicates/'
    save_suffix = ''

  } else {

    dupl_subdir = ''

  }

  for (city_vect in cities){

    city = city_vect[1]
    city_short = city_vect[2]
    city_i = match(list(city_vect), cities)

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
    malaria_data_sf = paste0(Dir_input_malaria_data, city, "/", subdir, '/malaria_data_', city_short, save_suffix, '.shp') %>%
      st_read() %>%
      mutate(ID = as.numeric(ID))

    malaria_data_sf = paste0(Dir_input_cov, city, '/', subdir, "/", city_short, "_full_cov.csv") %>%
      read.csv2(sep =";") %>%
      dplyr::select(-matches("dist")) %>%
      mutate(ID = as.numeric(ID)) %>%
      dplyr::right_join(malaria_data_sf, by=c("ID"= "ID")) %>%
      na.omit() %>%
      dplyr::select(geometry)

    # upload rf results
    RF_CovImp_pdp = paste0(Dir_final_results, method, "/RFE/", "RFE_", best_iter$nb_test, "/RF_output_", city_short, ".rds") %>%
      readRDS()

    # New data
    newdata = paste0(Dir_input_cov, city,"/prediction/",city_short,"_full_cov.csv") %>%
      read.csv2() %>%
      na.omit() %>%
      dplyr::select(best_var) %>%
      mutate(ID = as.numeric(ID))

    # Coordinates of new data
    coords = paste0(Dir_input_pts, city, "/", city_short, '_1km_pts/', city_short, "_1km_pts.shp") %>%
      st_read() %>%
      mutate(ID = as.numeric(idx))

    # Prediction
    # ---------------------------------------------------------------
    pred_rf_tot = lapply(RF_CovImp_pdp$models, function(x){predict(x, newdata=newdata)}$data)
    pred_rf_tot_df = do.call(cbind, pred_rf_tot) %>%
      as.data.frame()

    # Avg by row
    pred_avg_df = rowMeans(pred_rf_tot_df)
    pred_rf_cat = cbind("cat" = newdata$ID, "response" = pred_avg_df) %>%
      as.data.frame()

    # Data to save
    Data_output = cbind(newdata$ID, pred_avg_df, pred_rf_tot_df)
    colnames(Data_output) = c("cat","avg_response", paste0("response_model",1:25))

    rect_buffers = paste0(Dir_input_pts, city, "/", city_short, '_1km_vect/', city_short, "_1km_vect.shp") %>%
      st_read() %>%
      st_transform(crs =crs_REACT)

    # Load prediction grid
    grid_pred = paste0(Dir_input_pts,  city, "/",  city_short, "_1km_grid.tif") %>%
      raster() %>%
      projectRaster(crs = crs_REACT)

    grid_pred = rasterize(rect_buffers, grid_pred)

    new_raster = grid_pred
    values(new_raster) = NA

    for (i in newdata$ID){ # i =1
      df = pred_rf_cat %>% filter(cat == i)
      if(nrow(df)!=0){
        idx <- Which(grid_pred==i, cells=TRUE)
        raster::values(new_raster)[idx] <- pred_rf_cat[pred_rf_cat$cat == i,]$response
      }
    }

    new_raster = trim(new_raster)
    dir.create(paste0(Dir_output_pred, city, "/"), recursive = T, showWarnings = F)
    writeRaster(new_raster, paste0(Dir_output_pred, city, "/", city, "_", method, ".tif"), overwrite =T )

    # Plot 1 - plot raster 1km
    # ---------------------------------------------------------------------------$
    # transform raster into dataframe to plot in ggplot
    raster_pred_df <- new_raster %>%
      xyFromCell(1:ncell(new_raster)) %>%
      data.frame() %>%
      cbind(values(new_raster))

    colnames(raster_pred_df) = c("x","y","PfPr_Pred")

    if (city_short == "DES"){
      lp = c(0.8, 0.8)

    } else if (city_short == "Ouaga"){
      lp = c(0.2, 0.2)

    } else if (city_short == "Dakar"){
      lp = c( 0.15, 0.8)

    } else {
      lp = c(0.8, 0.15)
    }

    theme_set(theme_bw())

    th <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.text.x = element_text(size = rel(1.6), color = "grey25"),
                axis.text.y = element_text(size = rel(1.6), color = "grey25"),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                legend.text = element_text(size = rel(1), color = "grey25"),
                legend.title = element_text(size = rel(1.1), color = "grey25", margin = ggplot2::margin(t = 0, r = 0, b = 5, l = 0)),
                legend.key.height = unit(0.2, 'cm'),
                legend.key.width = unit(0.6, 'cm'),
                legend.position = lp)

    max_pfpr = round(max(raster_pred_df$PfPr_Pred, na.rm =T), 2)
    min_pfpr = round(min(raster_pred_df$PfPr_Pred, na.rm =T), 2)

    p = ggplot() + th +
      geom_tile(data = raster_pred_df , aes(x = x, y = y, fill = PfPr_Pred)) +
      coord_map(projection = "mercator")+  #coord_equal() +
      scale_y_continuous(expand = expansion(mult = c(0, 0), add = c(0, 0))) +
      scale_x_continuous(expand = expansion(mult = c(0, 0), add = c(0, 0))) +
      scale_fill_gradientn(name = expression(paste(italic("Pf"), Pr[2-10])), na.value = "white",
                           limits = c(min_pfpr, max_pfpr),
                           breaks = c(min_pfpr, max_pfpr),
                           colours = rev(terrain.colors(10))) +
      ggsn::scalebar(transform = TRUE, x.min = min(raster_pred_df$x), x.max = max(raster_pred_df$x),
               y.min =  min(raster_pred_df$y), y.max =  max(raster_pred_df$y),
               dist_unit = "km", dist = 2.5,
               location = "bottomright", height =0.01, border.size = 0, st.size = 3.5,st.dist=0.02)

    ggsave(paste0(Dir_output_pred, '/', city, "/",  method, "_1km_predictions_",city,".png"),
           height=8, width=12,
           units="in", dpi=350, plot= p, bg = "transparent")

    # Plot 2 - plot raster 1km  + malaria points distribution
    # ---------------------------------------------------------------------------
    # project malaria data
    malaria_sf_pr = malaria_data_sf %>%
      st_as_sf() %>%
      st_transform(crs(new_raster))

    p = ggplot() +
      th +
      geom_raster(data = raster_pred_df , aes(x = x, y = y, fill = PfPr_Pred)) +
      coord_quickmap() +
      coord_equal() +
      geom_sf(data = malaria_sf_pr) +
      coord_sf(datum = st_crs(4326)) +
      scale_fill_gradientn(name = expression(paste(italic("Pf"), Pr[2-10])), na.value = "white",
                           limits = c(0, as.integer(max_pfpr)),
                           breaks = c(0, as.integer(max_pfpr)),
                           colours = rev(terrain.colors(10)))

    ggsave(paste0(Dir_output_pred, city, '/', method, "_1km_predictions_with_pts_",city,".png"),
           height=8, width=12,
           units="in", dpi=350, plot= p, bg = "transparent")

    if (city_short %in% c("Kamp", "DES")){

      for (level in list("4", "5")){

        # Load admin
        # --------------------------------------------------
        admin_bd = paste0(Dir_input_admin, "/", city, "/", city_short, "_adm", level, ".shp") %>%
          st_read() %>%
          st_make_valid() %>%
          st_transform(crs(new_raster))

        # Aggregate value at admin level
        # --------------------------------------------------
        # raster to polygons
        raster_pred_shp = new_raster %>%
          rasterToPolygons() %>%
          st_as_sf()

        admin_bd_intersects = raster_pred_shp %>% st_intersects(admin_bd)
        admin_subs_sf = admin_bd[unlist(admin_bd_intersects),]

        # aggregate
        raster_pred_shp_agg =  raster_pred_shp %>%
          aggregate(admin_subs_sf,
                    FUN = mean,
                    join = function(x, y) st_is_within_distance(x, y, dist = 0.3))

        names(raster_pred_shp_agg) = c("PfPr_pred", "geometry")

        # Plot 3 - aggregated mean - admin
        # --------------------------------------------------
        p = ggplot() + th +
          geom_sf(data = raster_pred_shp_agg, aes(fill = PfPr_pred), size=.1) +
          scale_fill_gradientn(name = expression(paste("Mean  ",PfPr[2-10])), na.value = "white",
                               limits = c(min_pfpr, max_pfpr), colours = rev(terrain.colors(10)))

        ggsave(paste0(Dir_output_pred, city, '/', method, "_admin_", level, "_pred_",city,".png"),
               height=8, width=12,
               units="in", dpi=350, plot= p, bg = "transparent")

        # Plot 4 - Raster 1km + admin
        # --------------------------------------------------
        p = ggplot() + th +
          geom_raster(data = raster_pred_df , aes(x = x, y = y, fill = PfPr_Pred)) +
          coord_quickmap()+
          scale_fill_gradientn(name = expression(paste("Predicted ",PfPr[2-10])), na.value = "white",
                               limits = c(1,15), colours = rev(terrain.colors(10)))+
          geom_sf(data = admin_subs_sf, fill=NA, size =.1)

        ggsave(paste0(Dir_output_pred, city, '/', method, "_1_km_predictions__admin_", level, "_",city,".png"),
               height=8, width=12,
               units="in", dpi=350, plot= p, bg = "transparent")

      }

    }
  }

}
