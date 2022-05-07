#-------------------------------------------------------------------------------------------------------------------------
# Multi-satellite environmental and socio-economic predictors of vector-borne diseases in African cities: 
# malaria as an example
#
# This function extracts covariates using the coordinates of the malaria surveys for both the TRAINING and PREDICTION models
# It creates, for each city, a csv file for each group of covariates : 
# - The cosmo variables 
# - The Land cover variables 
# - The land use variables 
# - The LCZ variables 
# - The base variables 
#
# All the csv files are finally joined into a csv containing all covariates (full_cov.csv)
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities                 List with one vector c(city name, city short name, city epsg of projected crs) per city 
#' @param crs_REACT              Geographic coordinate reference system  
#' @param aim                    Aim of covariates extraction, either "training" or "prediction" 
#' @param n.cores                Number of cores used for parallelelization 
#' @param Dir_input_var          Path to input geospatial variables  
#' @param Dir_input_malaria_data Path to malaria database 
#' @param Dir_out_cov            Path to output covariates  
#' @param Dir_input_pts          Path to input points (pixel centroids) of prediction grid (Default NULL)
#
#' @return                       Nothing 
#' @export                       # Per city: A csv file containing extracted COSMO covariates 
#'                                           A csv file containing extracted LC covariates 
#'                                           A csv file containing extracted LU covariates 
#'                                           A csv file containing extracted LCZ covariates 
#'                                           A csv file containing extracted Base covariates 
#'                                           A csv file containing ALL extracted covariates (full_cov.csv)
#--------------------------------------------------------------------------------------------------------------------------


covariates.extraction = function(cities, crs_REACT, aim, n.cores, Dir_input_var, Dir_input_malaria_data, Dir_out_cov, Dir_input_pts = NULL){
  
  for (city_vect in cities){
    
    city = city_vect[1]
    city_short = city_vect[2]
    myepsg = city_vect[3] 
    if (aim == "prediction"){
      pred_pts = st_read(paste(Dir_input_pts, city, "/", city_short, '_1km_pts/', city_short, "_1km_pts.shp", sep=""))
      cov_path = paste(Dir_out_cov, city, '/prediction/', sep="") ; cov_path
      dir.create(file.path(cov_path), showWarnings = F, recursive = T)
      pred_pts = dplyr::rename(pred_pts, ID = idx)
      ref_grid = raster(paste(Dir_input_pts,  city, "/", city_short, "_1km_grid.tif", sep=""))
      
      #read rectangular buffer shapefile  
      rect_buffers = st_read(paste(Dir_input_pts, city, "/", city_short, '_1km_vect/', city_short, "_1km_vect.shp", sep=""))
      rect_buffers = st_as_sf(rect_buffers)
      rect_buffers = dplyr::rename(rect_buffers, ID = idx)
      
      proj_rect_buffers = st_transform(rect_buffers, crs = crs_REACT)
      pts_with_ID_transf = pred_pts
      #sometimes do not work in CRS react but projected CRS 
      pts_with_ID = st_transform(pts_with_ID_transf, crs=crs_REACT)
  
      
    } else if (aim == "training"){
      city_name = paste("malaria_data_", as.character(city_short), sep="")
      
      malaria_data_sf=st_read(paste(Dir_input_malaria_data, city, '/', city_name, '.shp', sep=""))
      cov_path = paste(Dir_out_cov, city, '/', sep="")
      dir.create(file.path(cov_path), showWarnings = F, recursive = T)
      
      #read buffer shapefile
      rect_buffers = st_read(paste(Dir_input_malaria_data, city, '/', city_name, '_1km_round.shp', sep=""))
      rect_buffers = st_as_sf(rect_buffers)
      rect_buffers = rect_buffers[,c("ID", "geometry")]
      proj_rect_buffers = st_transform(rect_buffers, crs = crs_REACT )
      pts_with_ID = malaria_data_sf[,c("ID")]
      
      #sometimes do not work in CRS react but projected CRS 
      pts_with_ID_transf = st_transform(pts_with_ID, crs=CRS(paste("+init=epsg:", myepsg, sep="")))
      
    }
    
  
    # -------------------------------
    #  COSMO
    # -------------------------------

    #initialize df with cosmo variables
    cosmo_cov = data.frame("ID" = pts_with_ID$ID)
    cosmo_stack = stack()
    setwd(paste(sep="", Dir_input_var, "COSMO/", city, "/"))
    #get all files in the cosmo directory
    files <- list.files(path = paste(sep="", Dir_input_var, "COSMO/", city, "/"), pattern="*.tif", full.names=FALSE, recursive=FALSE, ignore.case = FALSE)
    #iterate over each file
    for (file in files){
      if (endsWith(file, ".tif")){
        #get variable name and raster
        layer_name = sub(pattern=".tif.*", "", x=file)
        variable_name = sub(pattern = ".*DRY", "", x=layer_name)
        assign(x=variable_name, value = raster(file))
        raster_variable = get(variable_name)

        #assign crs
        crs(raster_variable) = crs_REACT

        #extract cosmo variable at malaria points
        extracted_cov <- cbind(st_drop_geometry(pts_with_ID), raster::extract(raster_variable,pts_with_ID))#res is 1km

        #add this newly extracted variable to big cosmo df
        cosmo_cov = join(cosmo_cov, extracted_cov, by = "ID")
        names(cosmo_cov)[length(names(cosmo_cov))]<-variable_name
      }
    }
    write.table(cosmo_cov, file=paste(cov_path, city_short, "_cosmo_cov.csv", sep=""), sep=";", dec=",", row.names = F)

    # -------------------------------
    #  LCZ
    # -------------------------------

    #create output directories
    binary_maps_path = paste(sep="", Dir_input_var, "LCZ/", city, "/binary_maps/")
    dist_path = paste(sep="", Dir_input_var, "LCZ/", city, "/distance_maps/")

    legend_classes = list(c("LCZ_compact", list(1,2,3)),
                          c("LCZ_open", list(4,5,6)),
                          c("LCZ_indu", list(8,10)),
                          c("LCZ_trees", list(11,12)),
                          c("LCZ_lowland", list(13,14)),
                          c("LCZ_informal", list(7)),
                          c("LCZ_sparse", list(9)),
                          c("LCZ_water", list(17)),
                          c("LCZ_mineral", list(15,16)),
                          c("LCZ_wetlands", list(18)))

    lcz_cov = data.frame("ID" = pts_with_ID$ID)

    for (class in legend_classes){
      #name of aggregated lcz class
      lcz_name = class[1]

      #get binary map
      binary_map = raster(paste0(binary_maps_path, lcz_name, ".tif"))
      if (st_crs(crs(binary_map)) != st_crs(CRS(paste("+init=epsg:", myepsg, sep="")))){
        binary_map = projectRaster(binary_map, crs = CRS(paste("+init=epsg:", myepsg, sep="")))
      }

      df_binary = data.frame(bool = values(binary_map)>0)
      len = length(df_binary[df_binary$bool == "TRUE",])

      if (len != 0){

        #extract proportion variables at malaria points
        lcz_prop_cov = cbind(st_drop_geometry(rect_buffers), exact_extract(binary_map, rect_buffers, fun = 'mean', progress = T))
        colnames(lcz_prop_cov) = c("ID", lcz_name)

        #get distance map
        lcz_dist = raster(paste0(dist_path, "Dist_", substring(lcz_name, 5), ".tif "))

        #extract distances variable at malaria points
        lcz_dist_cov = cbind(st_drop_geometry(rect_buffers), exact_extract(lcz_dist, rect_buffers, fun = 'mean', progress = T))
        colnames(lcz_dist_cov) = c("ID", paste0("Dist_", substring(lcz_name, 5)))

        lcz_cov$ID = as.numeric(lcz_cov$ID)
        lcz_prop_cov$ID = as.numeric(lcz_prop_cov$ID)
        lcz_dist_cov$ID = as.numeric(lcz_dist_cov$ID)

        lcz_cov = left_join(lcz_cov, lcz_prop_cov, by = c("ID"= "ID"))
        lcz_cov = left_join(lcz_cov, lcz_dist_cov, by = c("ID"= "ID"))

      }

    }

    write.table(lcz_cov, file=paste(cov_path, city_short, "_LCZ_cov.csv", sep=""), sep=";", dec=",", row.names=F)

    # -------------------------------
    #  LC
    # -------------------------------

    #NOTE HERE WE DONT WORK IN CRS REACT but in projected crs of the country (as distances are computed)

    #create output directories
    binary_maps_path = paste(sep="", Dir_input_var, "LULC/", city, "/LC/binary_maps/")
    dir.create(file.path(binary_maps_path), showWarnings = F)

    #get input LC map
    setwd(paste(sep="", Dir_input_var, "LULC/", city, "/LC/"))
    LC_raster = raster(paste(city_short, "_LC.tif", sep=""))

    if (city == "Kampala" | city == "kampala" | city == "Kamp") {legend_classes = list(c("LC_buildings", c(7)),
                                                         c("LC_trees", c(3)),
                                                         c("LC_low_veg", c(4)),
                                                         c("LC_water", c(2)),
                                                         c("LC_bare_ground", c(5)))
                                                         #c("LC_artificial_ground", c(6)) )

    } else {legend_classes =  list(c("LC_buildings", c(1, 111, 112, 113)),
                                                                c("LC_trees",  c(5)),
                                                                c("LC_low_veg", c(4)),
                                                                c("LC_water", c(2)),
                                                                c("LC_bare_ground", c(7)))
  }

    if (length(list.files(binary_maps_path)) == 0){
      empty_dir = T
    } else { empty_dir = F }

    lc_cov = foreach(class=iter(legend_classes), .combine = "join", .packages=c("raster", "exactextractr", "sf", "rgdal", "dplyr")) %dopar% {

      #name of aggregated lc class
      lc_name = class[1]
      #values of lcz sub classes
      lc_numeric = tail(class, -1)

      if(empty_dir){
        #create BINARY MAP
        binary_map= (LC_raster %in% lc_numeric)
        setwd(paste(sep="", Dir_input_var, "LULC/", city, "/LC/binary_maps/"))
        raster::writeRaster(x = binary_map, filename = lc_name, format="GTiff", overwrite = TRUE)

      } else {
        binary_map = raster(paste0(binary_maps_path, lc_name, ".tif"))
      }

      #create PROPORTIONS MAP
      lc_prop_cov = cbind(data.frame(sf::st_drop_geometry(rect_buffers)), data.frame(exactextractr::exact_extract(binary_map, rect_buffers, fun = 'mean', progress = T)))
      colnames(lc_prop_cov) = c("ID", lc_name)
      lc_prop_cov

  }
    write.table(lc_cov, file=paste(cov_path, city_short, "_LC_cov.csv", sep=""), sep=";", dec=",", row.names=F)

    # -------------------------------
    #  LU
    # -------------------------------

    #NOTE HERE WE DONT WORK IN CRS REACT but in projected crs of the country (as distances are computed)
    #get input LU map
    setwd(paste(sep="", Dir_input_var, "LULC/", city, "/LU/"))
    LU_map = st_read(paste(city_short, '_LU.shp', sep=""))

    if (st_crs(LU_map) != st_crs(CRS(paste("+init=epsg:", myepsg, sep="")))){
      LU_map = st_transform(LU_map, crs = myepsg)
    }

    legend_classes = list(c("LU_ACS", c("ACS")),
                          c("LU_planned", c("PLAN", "PLAN_HD", "PLAN_MD", "PLAN_LD")),
                          c("LU_informal", c("UNPLAN", "UNPLAN_HD", "UNPLAN_MD", "UNPLAN_LD", "DEPR")),
                          c("LU_wetlands", c("WET")))

    #When label is uncertain, replace with first label found by OBIA classif
    LU_map$MAP_LABEL = ifelse(LU_map$MAP_LABEL == "UNCERT", LU_map$FIRST_LABE, LU_map$MAP_LABEL)

    lu_cov = foreach(class=iter(legend_classes), .combine = "join", .packages=c("raster", "units", "sf", "sp", "dplyr")) %dopar% {

      #name of aggregated lc class
      lu_name = class[1]

      #values of lcz sub classes
      lu_values = tail(class, -1)

      #select lu classes
      selected_lu_map = st_make_valid(LU_map[LU_map$MAP_LABEL %in% lu_values,])
      selected_lu_map = st_as_sf(st_union(selected_lu_map,by_feature = F ))

      #select points in LU map extent
      LU_map_labels = LU_map[,c("MAP_LABEL")]
      pts_in_extent = pts_with_ID_transf[!is.na(over(as_Spatial(pts_with_ID_transf), as_Spatial(LU_map_labels))), ]
      rect_buffers_in_extent = rect_buffers[rect_buffers$ID %in% pts_in_extent$ID,]
      rect_buffers_in_extent$area = st_area(rect_buffers_in_extent)

      pts_not_in_extent = pts_with_ID_transf[is.na(over(as_Spatial(pts_with_ID_transf), as_Spatial(LU_map_labels))), ]
      rect_buffers_not_in_extent = st_drop_geometry(rect_buffers[rect_buffers$ID %in% pts_not_in_extent$ID,])
      if (nrow(rect_buffers_not_in_extent) != 0){rect_buffers_not_in_extent$proportion = NA}

      #intersect lu map with rect buffers
      lu_in_rect = st_intersection(selected_lu_map, rect_buffers_in_extent)
      lu_in_rect$intersect_area = st_area(lu_in_rect)

      #compute proportions
      lu_in_rect = st_drop_geometry(lu_in_rect)
      lu_in_rect$proportion = lu_in_rect$intersect_area/lu_in_rect$area
      pts_in_extent$ID = as.numeric(pts_in_extent$ID)
      lu_in_rect$ID = as.numeric(lu_in_rect$ID)
      lu_prop_cov = left_join(st_drop_geometry(pts_in_extent), lu_in_rect, by=c("ID" = "ID"))
      lu_prop_cov= lu_prop_cov[,c("ID", "proportion")]
      lu_prop_cov[is.na(lu_prop_cov)] = 0
      lu_prop_cov = rbind(lu_prop_cov, rect_buffers_not_in_extent)
      colnames(lu_prop_cov) = c("ID", lu_name)
      lu_prop_cov = drop_units(lu_prop_cov)
      lu_prop_cov

    }

     write.table(lu_cov, file=paste(cov_path, city_short, "_LU_cov.csv", sep=""), sep=";", dec=",", row.names=F)

    # -------------------------------
    #  BASE MODEL
    # -------------------------------

    setwd(paste(sep="", Dir_input_var, "BASE/", city, "/"))
    rect_buffers = rect_buffers[, c("ID", "geometry")]

    #get all files in the base directory
    files <- list.files(path = paste(sep="", Dir_input_var, "BASE/", city, "/"), pattern="*.tif", full.names=FALSE, recursive=FALSE, ignore.case = FALSE)
    base_cov = data.frame("ID" = pts_with_ID_transf$ID)

    for (file in files){
      if (endsWith(file, ".tif") & !startsWith(file, "Landsat_STD")){
        #get variable name and raster
        var_name = sub(pattern=paste("_", city_short, ".*", sep = ""), "", x=file)
        var_name = sub(pattern=paste(".*", "Landsat_", sep = ""), "", x=var_name)
        assign(x=var_name, value = raster(file))
        raster_variable = get(var_name)
        projected_raster = projectRaster(raster_variable, crs = CRS(paste("+init=epsg:", myepsg, sep="")))

        #extract in buffers
        base_cov_extracted = cbind(st_drop_geometry(rect_buffers), exact_extract(projected_raster, rect_buffers, fun = 'mean', progress = T)) #res is 100m
        base_cov = join(base_cov, base_cov_extracted, by = "ID")
        names(base_cov)[length(names(base_cov))]<-var_name
      }
    }

    write.table(base_cov, file=paste(cov_path, city_short, "_base_cov.csv", sep=""), sep=";", dec=",", row.names=F)

    # -------------------------------
    #  MERGE ALL CSV OF COVARIATES 
    # -------------------------------
    #merge all covariates in one DF
    
    files = list.files(path = cov_path, pattern="*.csv", full.names=FALSE, recursive=FALSE, ignore.case = FALSE)
      
    full_cov = data.frame("ID" = pts_with_ID_transf$ID)
    full_cov$ID = as.numeric(full_cov$ID)
    
    for (file in files) {
      if (! endsWith(file, "full_cov.csv")){
        cov_df = read.csv(paste0(cov_path, file), sep = ";") ; head(cov_df)
        cov_df$ID = as.numeric(cov_df$ID)
        full_cov = join(full_cov, cov_df, by = "ID")
        }
    }
     write.table(full_cov, file=paste(cov_path, city_short, "_full_cov.csv", sep=""), sep=";", dec=",", row.names=F)
  }
  
}
  

