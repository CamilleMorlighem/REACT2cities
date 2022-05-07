# ------------------------------------------------------------------------------------------------------------------
#  Modelling of the intra-urban malaria risk using socioeconomic and environmental factors with a random forest (RF) model
#  in two sub-Saharan cities: Kampala (Uganda) and Dar es Salaam (Tanzania)
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
#             - Environmental change variables
#             - Socioeconomic variables 
# 
# This function predicts pfpr2-10 using the best RFE models for each city and plots the predictions 
#
# It create for each city: 
# - a predictive map stored as a raster 
# - a png plot of the predictive map 
# - a csv file with the predicted values of the raster cells 
# - 6 different plots of the predictions (some as the raw predictive maps, others where the predictions are aggregated per admin area)
# It also creates a plot showing the raster cell predictions and aggregated predictions/adm area for the 2 cities 
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities                 List with one vector c(city name, city short name, city epsg of projected crs) per city 
#' @param crs_REACT              Geographic coordinate reference system  
#' @param var_group              Type of RF model to use for predictions. Either "all3Geo" (uses all LULC, LCZ and COSMO covariates) or "BEST" (uses best selected covariates) 
#' @param input_cov              Path to input covariate tables 
#' @param input_malaria          Path to input selected malaria data points
#' @param output_rf              Path to input RF models
#' @param output_pred            Path to predicted outputs 
#' @param Dir_input_pts          Path to output RF models
#' @param input_admin            Path to output RF models
#' @param mode                   If var_group is "all3geo", character indicating whether to use the RFE model (and the number of the iteration i) ("RFE-[i]") or all covariates ("ALL"). Default is NULL.

#' @return                       Nothing 
#' @export                       # Per city: A predictive map stored as a raster 
#                                            A png plot of the predictive map 
#                                            A csv file with the predicted values of the raster cells 
#                                            6 different plots of the predictions (some as the raw predictive maps, others where the predictions are aggregated per admin area)
#                                A plot with the raster cell predictions and aggregated predictions of the 2 cities 
#--------------------------------------------------------------------------------------------------------------------------


rf.prediction = function(cities, crs_REACT, input_cov, input_malaria, output_rf, output_pred, Dir_input_pts, input_admin, var_group, mode = NULL){
  
  
  cv_rf = "SpCVHpt" ; HYt_mode = "HYt_insideCV" ; 
  
  for (city_nb in 1:length(cities)){
    
    # Define city & input file path & output file name
    # ------------------------------------------------------------------------------------------------------------------
    city = cities[city_nb][[1]][1]
    city_short = cities[city_nb][[1]][2]
    
    # define output file and name
    # ------------------------------------------------------------------------------------------------------------------
    output_file_name = paste0(cv_rf,"_",var_group,"_",city) ; output_file_name 
    output_file = paste0(output_pred,city, '/') ; output_file
    
    # Check for file existence and create it if needed
    if(!dir.exists(file.path(output_file))){ dir.create(file.path(output_file), recursive = TRUE)   }
    
    # Load Malaria and covariates data
    # ------------------------------------------------------------------------------------------------------------------
    malaria_sf = st_read(paste0(input_malaria, city, "/malaria_data_",city_short,".shp"))
    malaria_sf <- malaria_sf[,c("ID","PfPR2_10")] ; malaria_sf$ID = as.numeric(as.character(malaria_sf$ID))
    
    # store the coordinates (lon and lat) of malaria data
    coords = as.data.frame(cbind(st_coordinates(malaria_sf), ID = malaria_sf$ID))
    
    # Load covariates csv into a dataframe
    cov_dataset = read.csv2(paste0(input_cov, city, "/", city_short, "_full_cov.csv"), sep =";")
    var_group_test_list = colnames(cov_dataset)
    cov_dataset = dplyr::left_join(st_drop_geometry(malaria_sf), cov_dataset, by=c("ID"= "ID"))
    
    # ------------------------------------------------------------------------------------------------------------------
    #       Variable selection
    # ------------------------------------------------------------------------------------------------------------------
    # Select variable according the variable list of the model to be assessed
    
    # if (var_group == 'BEST' | (var_group == "all3Geo" & mode != "ALL")){
    #   var_group_test_list = best_var[city_nb][[1]]
    #   var_group_test_list = c(var_group_test_list, "ID")
    #   
    # } 
    
    if (var_group=="all3Geo"){

      if (mode == "ALL"){
        rf_res = read.csv2(paste0(output_rf, var_group, "/", city, "/", mode, "/RFResults_all3Geo_", mode,"_", city, ".csv"))

      } else {
        rf_res = read.csv2(paste0(output_rf, var_group, "/", city, "/", "RFE", "/RFResults_all3Geo_", mode,"_", city, ".csv"))

      }

    } else {
      rf_res = read.csv2(paste0(output_rf,"BEST/", city, "/Results_",output_file_name,".csv"))
      #RF_CovImp_pdp = readRDS(paste0(output_rf,"BEST/", city, "/RFoutput_",output_file_name,"_celia.rds"))

    }
    
    var_group_test_list = rf_res$covariates[[1]]
    var_group_test_list = strsplit(var_group_test_list, ",")[[1]]
    var_group_test_list = c(var_group_test_list, "ID")
    
    variables= dplyr::select(cov_dataset,PfPR2_10, all_of(var_group_test_list))
    var_list = colnames(variables) ; var_list
    
    # Verify for NA atfer variable selection - and update the data
    data_updated = na.omit(cov_dataset[,c("ID",var_list)])
    coords_updated = coords[which(coords$ID %in%  data_updated$ID),c("X","Y")]
    variables = data_updated[,var_list] ; colnames(variables)
    
    # define predictor list - variable list without depend variable
    predictor_dataset = dplyr::select(variables,-PfPR2_10)
    predictor_list = colnames(predictor_dataset) ; predictor_list
    
    # ------------------------------------------------------------------------------------------------------------------
    #       Prediction
    # ------------------------------------------------------------------------------------------------------------------
    # Load data 
    # ---------------------------------------------------------------
    # upload rf results 
    if (var_group=="all3Geo"){
      
      if (mode == "ALL"){
        RF_CovImp_pdp = readRDS(paste0(output_rf, var_group, "/", city, "/", mode, "/RFoutput_all3Geo_", mode,"_", city, ".rds"))
        
      } else {
        RF_CovImp_pdp = readRDS(paste0(output_rf, var_group, "/", city, "/", "RFE", "/RFoutput_all3Geo_", mode,"_", city, ".rds"))
        
      }
      
    } else {
      RF_CovImp_pdp = readRDS(paste0(output_rf,"BEST/", city, "/RFoutput_",output_file_name,".rds"))
      #RF_CovImp_pdp = readRDS(paste0(output_rf,"BEST/", city, "/RFoutput_",output_file_name,"_celia.rds"))
      
    }
    
    # New data 
    newdata = read.csv2(paste0(input_cov,city,"/prediction/",city_short,"_full_cov.csv"), sep=';', dec=",") ; newdata
    #newdata = read.csv2(paste0(input_cov,city,"/prediction/",city_short,"_full_cov_celia.csv"), sep=';', dec=",") ; newdata$ID = newdata$cat ; newdata = select(newdata, -c("cat"))
    newdata = na.omit(newdata) ; nrow(newdata)
    newdata = newdata[, var_group_test_list]
    
    coords = st_read(paste(Dir_input_pts, city, "/", city_short, '_1km_pts/', city_short, "_1km_pts.shp", sep=""))
    coords$idx = as.numeric(coords$idx) ; newdata$ID = as.numeric(newdata$ID)
    
    #newdata_coords = st_as_sf(left_join(newdata, coords, by = c("ID" = "idx")))
    #newdata_with_id = newdata
    
    # Load grid 
    grid_pred = raster(paste0(Dir_input_pts,city,"/GRID_1KM_", city_short))
    
    # Prediction 
    # ---------------------------------------------------------------
    pred_rf_tot = lapply(RF_CovImp_pdp$models, function(x){predict(x, newdata=newdata)}$data)
    pred_rf_tot_df = as.data.frame(do.call(cbind,pred_rf_tot))
    
    # Avg by row
    pred_avg_df = rowMeans(pred_rf_tot_df)
    pred_rf_cat = as.data.frame(cbind("cat" = newdata$ID, "response" = pred_avg_df))
    
    # Data to save 
    Data_output = cbind(newdata$ID, pred_avg_df, pred_rf_tot_df)
    colnames(Data_output) = c("cat","avg_response",paste0("response_model",1:50))
    write.csv2(Data_output, file = paste0(output_file,"Data_prediction_",output_file_name,".csv"))
    
    rect_buffers = st_read(paste(Dir_input_pts, city, "/", city_short, '_1km_vect/', city_short, "_1km_vect.shp", sep=""))
    rect_buffers = st_transform(rect_buffers, crs =crs_REACT)
    grid_pred = projectRaster(grid_pred, crs = crs_REACT)
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
    
    writeRaster(new_raster,paste0(output_pred, city, '/map_pred_avg_', city_short,".tif"), overwrite=T) 
    
    # plot 
    png(paste0(output_file,"Map_pred_AVG_",output_file_name,".png"), 
        width = 880, height = 880, units = "px", pointsize = 12)
    par(mar=c(c(5, 4, 4, 6)))
    plot(new_raster, cex.axis = 2, 
         legend.width=1.3,
         axis.args=list(at=seq(from = 0, by = 20/10, to = 25), cex.axis=2),
         legend.args=list(text='PfPR2-10', side = 3, line=1.2, cex=1.8))
    dev.off()
    par(mar=c(c(5, 4, 4, 2)))
    
    # To run on each model 
    # ---------------------------------------------------------
    # for(h in 1:length(RF_CovImp_pdp$models)){   # h=1
    #   
    #   # predict 
    #   pred_rf = predict(RF_CovImp_pdp$models[[h]], newdata = newdata[,2:ncol(newdata)])
    #   pred_rf_cat = cbind("cat" = newdata$ID,pred_rf$data)
    #   
    #   new_raster = grid_pred
    #   values(new_raster) = NA
    #   for (i in newdata$ID){ # i =1
    #     #ID = match(i,pred_rf_cat$cat)
    #     df = pred_rf_cat %>% filter(cat == i)
    #     if(nrow(df)!=0){
    #       idx <- Which(grid_pred==i, cells=TRUE) 
    #       raster::values(new_raster)[idx] <- pred_rf_cat[pred_rf_cat$cat == i,]$response
    #       #grid_pred[i] = pred_rf_cat$response[i]
    #     }
    #   }
    #   
    #   new_raster = trim(new_raster)
    #   plot(new_raster)
    #   
    #   # write raster
    #   writeRaster(new_raster,paste0(output_file,"Raster_pred_",h,"_",output_file_name,".tif"))
    #   
    #   # Plot map 
    #   png(paste0(output_file,"Map_pred_",h,"_",output_file_name,".png"), 
    #       width = 880, height = 880, units = "px", pointsize = 12)
    #   par(mar=c(c(5, 4, 4, 6)))
    #   plot(new_raster, cex.axis = 2, 
    #        legend.width=1.3,
    #        axis.args=list(at=seq(from = 0, by = 20/10, to = 25), cex.axis=2),
    #        legend.args=list(text='PfPR2-10', side = 3, line=1.2, cex=1.8))
    #   dev.off()
    #   par(mar=c(c(5, 4, 4, 2)))
    #   
    # }
    
    # -------------------------------------------------
    # Plot prediction 
    # -------------------------------------------------
    
    # Load raster with prediction - 1km raster - average value over 50 models  
    # --------------------------------------------------
    # load as raster 
    raster_pred = raster(paste0(output_pred, city, '/map_pred_avg_', city_short,".tif"))
    # load with Stars package
    raster_pred_stars = read_stars(paste0(output_pred, city, '/map_pred_avg_', city_short,".tif"))
    
    # transform raster into dataframe to plot in ggplot
    xy <- data.frame(xyFromCell(new_raster, 1:ncell(new_raster)))
    raster_pred_df = cbind(xy, values(new_raster))
    colnames(raster_pred_df) = c("x","y","PfPr_Pred")
    
    # Load admin 
    # --------------------------------------------------
    if(city_short == 'DES'){
      admin_4 = st_make_valid(st_read(paste0(input_admin, city, "/TZA_adm4_2002.shp")))
      admin_4_crop = st_make_valid(st_read(paste0(input_admin, city, "/DES_adm4_2002.shp")))
      admin_5 = st_make_valid(st_read(paste0(input_admin, city, "/DES_adm5.shp")))
      admin_5 = st_transform(admin_5, crs_REACT)
      
      
    } else if(city == 'Kampala'){
      admin_4 = st_make_valid(st_read(paste0(input_admin, city, "/UGA_adm4_2002.shp")))
      admin_4_crop = st_make_valid(st_read(paste0(input_admin, city, "/Kamp_adm4_2002.shp")))
      admin_5 = st_make_valid(st_read(paste0(input_admin, city, "/Kamp_adm5.shp")))
      admin_5 = st_transform(admin_5, crs_REACT)
    }
    
    admin_4_pr = st_make_valid(st_transform(admin_4, crs_REACT)) 
    admin_4_crop_pr = st_transform(admin_4_crop, crs_REACT)
    
    # Aggregate value at admin level 5
    # --------------------------------------------------
    # raster to shp 
    raster_pred_shp = rasterToPolygons(raster_pred)
    
    # aggregate 
    raster_pred_shp_agg = aggregate(st_as_sf(raster_pred_shp), admin_5$geometry, FUN = mean, 
                                    join = function(x, y) st_is_within_distance(x, y, dist = 0.3))
    names(raster_pred_shp_agg) = c("PfPr_pred", "geometry")
    
    
    # Aggregate value at admin level 4
    # --------------------------------------------------
    # raster to shp 
    raster_pred_shp = rasterToPolygons(raster_pred)
    raster_pred_sf = st_as_sf(raster_pred_shp)
    
    admin_intersect_sf = st_intersects(raster_pred_sf, admin_4_pr)
    # subsetting
    admin_subs_sf = admin_4_pr[unlist(admin_intersect_sf),]
    
    # aggregate 
    raster_pred_shp_agg4 = aggregate(st_as_sf(raster_pred_shp), admin_subs_sf$geometry, FUN = mean, 
                                     join = function(x, y) st_is_within_distance(x, y, dist = 0.3))
    names(raster_pred_shp_agg4) = c("PfPr_pred", "geometry")
    
    
    # --------------------------------------------------
    # Plots  
    # --------------------------------------------------
    max_pfpr = as.integer(floor(max(raster_pred_df$PfPr_Pred, na.rm =T)))
    min_pfpr = as.integer(floor(min(raster_pred_df$PfPr_Pred, na.rm =T)))
    max_pfpr = round(max(raster_pred_df$PfPr_Pred, na.rm =T))
    min_pfpr = round(min(raster_pred_df$PfPr_Pred, na.rm =T))
    
    # Define output file 
    output_file_plot = paste0(output_file, "/final_plots/")
    
    # Check for file existence and create it if needed
    if(!dir.exists(file.path(output_file_plot))){ dir.create(file.path(output_file_plot), recursive = TRUE)   }
    
    # define theme used in ggplot
    th <- 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.text.x = element_text(size=16),
            axis.text.y = element_text(size=16), 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            legend.text = element_text(size=15),
            legend.title = element_text(size=17) ) 
    
    #max_pfpr = 15
    
    # Plot 1 - aggregated mean - admin 5 
    # --------------------------------------------------
    
    p1 = ggplot() + th + 
      geom_sf(data = raster_pred_shp_agg, aes(fill = PfPr_pred), size=.1) + 
      scale_fill_gradientn(name = expression(paste("Mean ", italic("Pf"), Pr[2-10])), na.value = "white", 
                           limits = c(0, max_pfpr),
                           breaks = c(0, max_pfpr), colours = rev(terrain.colors(10)))
    
    
    
    ggsave(paste0("Plot_MeanPfPr_Admin5_",city,".png"),
           height=8, width=12,
           units="in", dpi=350, plot= p1, bg = "transparent",
           path = paste(output_file_plot))
    
    # Plot 2 - Raster 1km + admin 5 
    # --------------------------------------------------
    p2 = ggplot() + th + 
      geom_raster(data = raster_pred_df , aes(x = x, y = y, fill = PfPr_Pred)) +
      coord_quickmap()+
      scale_fill_gradientn(name = expression(paste("Mean ", italic("Pf"), Pr[2-10])), na.value = "white", limits = c(0, max_pfpr),
                           breaks = c(0, max_pfpr), colours = rev(terrain.colors(10))) +
      geom_sf(data = admin_5,fill=NA, size =.1)
    
    ggsave(paste0("Plot_PredPfPr_raster_Admin5_",city,".png"),
           height=8, width=12,
           units="in", dpi=350, plot= p2, bg = "transparent",
           path = paste(output_file_plot))
    
    # Plot 3 - Plot raster 1km 
    # ---------------------------
    
    p3 = ggplot() + th +
      geom_raster(data = raster_pred_df , aes(x = x, y = y, fill = PfPr_Pred)) +
      coord_quickmap()+
      scale_fill_gradientn(name = expression(paste("Mean ", italic("Pf"), Pr[2-10])), na.value = "white", 
                           limits = c(0, max_pfpr),
                           breaks = c(0, max_pfpr), colours = rev(terrain.colors(10)))
    
    ggsave(paste0("Plot_PredPfPr_raster_",city,".png"),
           height=8, width=12,
           units="in", dpi=350, plot= p3, bg = "transparent",
           path = paste(output_file_plot))
    
    
    # Plot 4 - plot raster 1km  + malaria points distribution 
    # ---------------------------------------------------------------------------
    
    # project malaria data 
    malaria_sf_pr = st_transform(malaria_sf, crs(raster_pred))
    
    p4 = ggplot() + th +
      geom_raster(data = raster_pred_df , aes(x = x, y = y, fill = PfPr_Pred)) +
      coord_quickmap()+
      geom_sf(data = malaria_sf_pr) +
      scale_fill_gradientn(name = expression(paste("Mean ", italic("Pf"), Pr[2-10])), na.value = "white", 
                           limits = c(0, max_pfpr),
                           breaks = c(0, max_pfpr), colours = rev(terrain.colors(10)))
    
    ggsave(paste0("Plot_PredPfPr_raster_malariaPt_",city,".png"),
           height=8, width=12,
           units="in", dpi=350, plot= p4, bg = "transparent",
           path = paste(output_file_plot))
    
    
    # Plot 5 - aggregated mean - admin 4
    # ---------------------------------------------------------------------------
    
    p5 = ggplot() + th + 
      geom_sf(data = raster_pred_shp_agg4, aes(fill = PfPr_pred), size=.1) + 
      scale_fill_gradientn(name = expression(paste("Mean ", italic("Pf"), Pr[2-10])), na.value = "white", 
                           limits = c(0, max_pfpr),
                           breaks = c(0, max_pfpr), colours = rev(terrain.colors(10)))
    
    
    ggsave(paste0("Plot_MeanPfPr_Admin4_",city,".png"),
           height=8, width=12,
           units="in", dpi=350, plot= p5, bg = "transparent",
           path = paste(output_file_plot))
    
    
    # Plot 6 - Raster 1km + admin 4 
    # --------------------------------------------------
    p6 = ggplot() + th + 
      geom_raster(data = raster_pred_df , aes(x = x, y = y, fill = PfPr_Pred)) +
      coord_quickmap()+
      scale_fill_gradientn(name = expression(paste("Mean ", italic("Pf"), Pr[2-10])), na.value = "white", 
                           limits = c(0, max_pfpr),
                           breaks = c(0, max_pfpr), colours = rev(terrain.colors(10))) +
      geom_sf(data = admin_subs_sf,fill=NA, size =.1)
    
    ggsave(paste0("Plot_PredPfPr_raster_Admin4_",city,".png"),
           height=8, width=12,
           units="in", dpi=350, plot= p6, bg = "transparent",
           path = paste(output_file_plot))
    
    
    
    if (city_short =="DES"){
      
      pp1 = p4 ; pp2 = p1
      
    } else {
      
      pp3 = p4 ; pp4 = p1
      
    }
    
    
    
  }  
  
  
  # gg1 = ggarrange(pp1, pp2, 
  #                 labels = c("a)", "b)"),
  #                 align = c("v"),
  #                 ncol = 2, nrow = 1, 
  #                 common.legend = T, legend = F)
  # 
  # gg2 = ggarrange(pp3, pp4, 
  #                 labels = c("c)", "d)"),
  #                 align = c("v"),
  #                 ncol = 2, nrow = 1, 
  #                 common.legend = T, legend = F)
  # 
  # fig = ggpubr::ggarrange(gg1, gg2,
  #                         #labels = c("a)", "b)", "c)", "d)"),
  #                         align = c("h"),
  #                         ncol = 1, nrow = 2, 
  #                         common.legend = T, legend = "right") 
  
  fig = ggpubr::ggarrange(pp1, pp2, pp3, pp4, 
                          labels = c("a)", "b)", "c)", "d)"),
                          #align = c("h"),
                          ncol = 2, nrow = 2, 
                          common.legend = T, legend = "right") 
  
  
  ggsave(paste0("Pred_map",".tiff"),
         height=18, width=20,
         units="in", dpi=300, plot= fig, bg = "transparent",
         path = paste(output_pred)) 
  
  
}
  
 