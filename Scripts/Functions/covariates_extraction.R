#-------------------------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# This function extracts covariates using the coordinates of the malaria surveys for both the TRAINING and PREDICTION models
# It creates, for each city, a csv file for each group of covariates :
# - The cosmo variables
# - The Land cover variables
# - The land use variables
# - The LCZ variables
# - The base/ancillary variables
# - The climatic variables
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
#' @param Dir_pred_grid          Path to input points (pixel centroids) of prediction grid (Default NULL)
#' @param method                 Method to implement, either "m0", "m1", "m2" or "TM" (useful for simulation study) (Default "m0")
#' @param simulation             Boolean indicates whether simulation study is conducted (Default F)
#
#' @return                       Nothing
#' @export                       # Per city: A separate csv file containing extracted covariates for each set of covariates
#'                                           A csv file containing ALL extracted covariates (full_cov.csv)
#--------------------------------------------------------------------------------------------------------------------------


covariates.extraction = function(cities, crs_REACT, aim, n.cores, Dir_input_var, Dir_input_malaria_data, Dir_out_cov, Dir_pred_grid = NULL, method = "m0", simulation = F){

  for (city_vect in cities){

    city = city_vect[1]
    city_short = city_vect[2]
    myepsg = city_vect[3]

    if (aim == "training"){

      city_name = paste0("malaria_data_", as.character(city_short))
      subdir = '/'
      save_suffix = ''

     if (simulation) {

        subdir = 'simul_displ/'

        if (method == "m0") { #use displaced data
          city_name = paste0(city_name, "_displ")
          save_suffix = "_displ"
        }

      }

      if (method == "m1" | method == "m2"){

        subdir = paste0(subdir, 'duplicates/')

      }

      data_pts_4326 = paste0(Dir_input_malaria_data, city, '/', subdir, city_name, '.shp') %>%
        st_read() %>%
        dplyr::select(ID)

      data_pts = st_transform(data_pts_4326, crs=CRS(paste("+init=epsg:", myepsg, sep="")))

      cov_path = paste0(Dir_out_cov, city, '/', subdir)
      dir.create(file.path(cov_path), showWarnings = F)

      #read buffer shapefile
      buffer = paste0(Dir_input_malaria_data, city, '/', subdir, city_name, '_1km_round.shp') %>%
        st_read() %>%
        st_as_sf() %>%
        dplyr::select(ID, geometry)

      buffer_4326 = st_transform(buffer, crs = crs_REACT)

    } else if (aim == "prediction"){

      data_pts_4326 = paste0(Dir_pred_grid, city, "/", city_short, '_1km_pts/', city_short, "_1km_pts.shp") %>%
        st_read() %>%
        dplyr::rename(ID = idx)

      data_pts = st_transform(data_pts_4326, crs=CRS(paste("+init=epsg:", myepsg, sep="")))

      cov_path = paste0(Dir_out_cov, city, '/prediction/')
      dir.create(cov_path, showWarnings = F, recursive = T)
      ref_grid = paste0(Dir_pred_grid, city, "/", city_short, "_1km_grid.tif") %>% raster()

      #read rectangular buffer shapefile
      buffer = paste0(Dir_pred_grid, city, "/", city_short, '_1km_vect/', city_short, "_1km_vect.shp") %>%
        st_read() %>%
        st_as_sf() %>%
        dplyr::rename(ID = idx)

      buffer_4326 = buffer %>% st_transform(crs = crs_REACT)

    }

    # -------------------------------
    #  COSMO
    # -------------------------------

    if (city == "Kampala"||city=="Dar es salaam"){

      #initialize df with cosmo variables
      cosmo_cov = data.frame("ID" = data_pts_4326$ID)
      cosmo_stack = stack()

      #get all files in the cosmo directory
      files <- list.files(path = paste0(Dir_input_var, "COSMO/", city, "/"), pattern = ".tif$", full.names=T)

      #iterate over each file
      for (file in files){

        #get variable name and raster
        layer_name = sub(pattern=".tif.*", "", x = file)
        variable_name = sub(pattern = ".*DRY", "", x = layer_name)

        assign(x=variable_name, value = raster(file))
        raster_variable = get(variable_name)

        #assign crs
        crs(raster_variable) = crs_REACT

        #extract cosmo variable at malaria points
        extracted_cov <- data_pts_4326 %>%
            st_drop_geometry() %>%
            cbind(raster::extract(raster_variable, data_pts_4326))

        #add this newly extracted variable to big cosmo df
        cosmo_cov = join(cosmo_cov, extracted_cov, by = "ID")
        names(cosmo_cov)[length(names(cosmo_cov))]<-variable_name

      }

      write.table(cosmo_cov, file = paste0(cov_path, city_short, "_cosmo_cov", save_suffix, ".csv"), sep=";", dec=",", row.names = F)

    }

    # -------------------------------
    #  LCZ
    # -------------------------------

    #output directories
    binary_maps_path = paste0(Dir_input_var, "LCZ/", city, "/binary_maps/")
    dist_path = paste0(Dir_input_var, "LCZ/", city, "/distance_maps/")

    legend_classes = list(c("LCZ_compact", list(1,2,3)),
                          c("LCZ_open", list(4,5,6)),
                          c("LCZ_indu", list(8,10)),
                          c("LCZ_trees", list(11,12)),
                          c("LCZ_lowland", list(13,14)),
                          c("LCZ_water", list(17)))
    # c("LCZ_informal", list(7)),
    # c("LCZ_sparse", list(9)),
    # c("LCZ_water", list(17)),
    # c("LCZ_mineral", list(15,16)),
    # c("LCZ_wetlands", list(18)))

    lcz_cov = data.frame("ID" = as.numeric(data_pts_4326$ID))

    for (class in legend_classes){

      #name of aggregated lcz class
      lcz_name = class[1]

      #get binary map
      binary_map = paste0(binary_maps_path, lcz_name, ".tif") %>% raster()

      #extract proportion variables at malaria points
      lcz_prop_cov = buffer_4326 %>%
        st_drop_geometry() %>%
        mutate(ID = as.numeric(ID)) %>%
        cbind(exact_extract(binary_map, buffer_4326, fun = 'mean', progress = T))

      colnames(lcz_prop_cov) = c("ID", lcz_name)

      #get distance map
      lcz_dist = paste0(dist_path, "dist_", lcz_name, ".tif ") %>% raster()

      #extract distances variable at malaria points
      lcz_dist_cov = buffer %>%
        st_drop_geometry %>%
        mutate(ID = as.numeric(ID)) %>%
        cbind(exact_extract(lcz_dist, buffer, fun = 'mean', progress = T))

      colnames(lcz_dist_cov) = c("ID", paste0("dist_", lcz_name))

      lcz_cov = lcz_cov %>%
        left_join(lcz_prop_cov, by = c("ID"= "ID")) %>%
        left_join(lcz_dist_cov, by = c("ID"= "ID"))

    }

    write.table(lcz_cov, file=paste(cov_path, city_short, "_LCZ_cov", save_suffix, ".csv", sep=""), sep=";", dec=",", row.names=F)

    # -------------------------------
    #  LC
    # -------------------------------

    #NOTE HERE WE DONT WORK IN CRS REACT but in projected crs of the country (as distances are computed)

    #create output directories
    binary_maps_path = paste0(Dir_input_var, "LULC/", city, "/LC/binary_maps/")
    dir.create(binary_maps_path, showWarnings = F)

    if (city == "Dakar") {
      legend_classes = list(c("LC_buildings", c(12, 111, 112, 113)),
                            c("LC_trees", c(22)),
                            c("LC_low_veg", c(23)),
                            c("LC_water", c(33)),
                            c("LC_bare_ground", c(45)))

    } else if (city == "Kampala") {
      legend_classes = list(c("LC_buildings", c(7)),
                            c("LC_trees", c(3)),
                            c("LC_low_veg", c(4)),
                            c("LC_water", c(2)),
                            c("LC_bare_ground", c(5)))

    } else if (city == "Dar es salaam") {
      legend_classes =  list(c("LC_buildings", c(1, 111, 112, 113)),
                             c("LC_trees",  c(5)),
                             c("LC_low_veg", c(4)),
                             c("LC_water", c(2)),
                             c("LC_bare_ground", c(7)))

    } else {
      legend_classes =  list(c("LC_buildings", c(11, 112, 111)),
                             c("LC_trees", c(31)),
                             c("LC_low_veg", c(30)),
                             c("LC_water", c(41)),
                             c("LC_bare_ground", c(20)))}

    lc_cov = foreach(class=iter(legend_classes), .combine = "join", .packages=c("raster", "exactextractr", "sf", "rgdal", "dplyr")) %dopar% {

      #name of aggregated lc class
      lc_name = class[1]
      print(lc_name)

      #get binary map
      binary_map = paste0(binary_maps_path, lc_name, ".tif") %>% raster()
      print('starts extraction')
      #create PROPORTIONS MAP
      lc_prop_cov = buffer %>%
        st_drop_geometry() %>%
        cbind(exact_extract(binary_map, buffer, fun = 'mean', progress = T))

      colnames(lc_prop_cov) = c("ID", lc_name)

      lc_prop_cov

    }

    write.table(lc_cov, file = paste0(cov_path, city_short, "_LC_cov", save_suffix, ".csv"), sep=";", dec=",", row.names=F)

    # -------------------------------
    #  LU
    # -------------------------------
    #NOTE HERE WE DONT WORK IN CRS REACT but in projected crs of the country (as distances are computed)
    #get input LU map
    LU_map = paste0(Dir_input_var, "LULC/", city, "/LU/", city_short, '_LU.shp') %>%
      st_read()

    legend_classes = list(c("LU_ACS", c("ACS")),
                          c("LU_planned", c("PLAN", "PLAN_HD", "PLAN_MD", "PLAN_LD")),
                          class = c("LU_informal", c("UNPLAN", "UNPLAN_HD", "UNPLAN_MD", "UNPLAN_LD", "DEPR")),
                          c("LU_wetlands", c("WET")))

    #When label is uncertain, replace with first label found by OBIA classif
    LU_map$MAP_LABEL = ifelse(LU_map$MAP_LABEL == "UNCERT", LU_map$FIRST_LABE, LU_map$MAP_LABEL)

    lu_cov = foreach(class = iter(legend_classes), .combine = "join", .packages=c("raster", "units", "sf", "sp", "dplyr")) %dopar% {

      #name of aggregated lc class
      lu_name = class[1]

      #values of lcz sub classes
      lu_values = tail(class, -1)

      #select lu classes
      selected_lu_map = LU_map %>%
        st_make_valid() %>%
        filter(MAP_LABEL %in% lu_values) %>%
        st_union(by_feature = F ) %>%
        st_make_valid() %>%
        st_as_sf()

      #select points in LU map extent
      LU_map_labels = LU_map %>%
        dplyr::select(MAP_LABEL) %>%
        as_Spatial()

      is_pt_in_extent = data_pts %>%
        as_Spatial() %>%
        over(LU_map_labels) %>%
        is.na()

      pts_in_extent = data_pts[! is_pt_in_extent, ]

      buffer_in_extent = buffer %>%
        filter(ID %in% pts_in_extent$ID)
      buffer_in_extent$area = st_area(buffer_in_extent)

      pts_not_in_extent = data_pts[is_pt_in_extent, ]
      buffer_not_in_extent = buffer %>%
        filter(ID %in% pts_not_in_extent$ID) %>%
        st_drop_geometry()

      if (nrow(buffer_not_in_extent) != 0) {
        buffer_not_in_extent$proportion = NA
      }

      #intersect lu map with rect buffers
      lu_in_rect = selected_lu_map %>%
        st_buffer(0) %>%
        st_intersection(buffer_in_extent)

      lu_in_rect$intersect_area = st_area(lu_in_rect)

      #compute proportions
      lu_in_rect = lu_in_rect %>%
        st_drop_geometry() %>%
        mutate(proportion = intersect_area/area, ID = as.numeric(ID))

      lu_prop_cov = pts_in_extent %>%
        mutate(ID = as.numeric(ID)) %>%
        st_drop_geometry() %>%
        left_join(lu_in_rect, by=c("ID" = "ID")) %>%
        dplyr::select(ID, proportion)

      lu_prop_cov[is.na(lu_prop_cov)] = 0
      lu_prop_cov = lu_prop_cov %>%
        rbind(buffer_not_in_extent) %>%
        units::drop_units()

      colnames(lu_prop_cov) = c("ID", lu_name)

      lu_prop_cov

    }

    write.table(lu_cov, file=paste0(cov_path, city_short, "_LU_cov", save_suffix, ".csv"), sep=";", dec=",", row.names=F)

    # -------------------------------
    #  BASE MODEL
    # -------------------------------
    buffer = buffer %>%
      dplyr::select(ID, geometry)

    #get all files in the cosmo directory
    files <- list.files(path = paste0(Dir_input_var, "BASE/", city, "/"), pattern="*.tif$", full.names=F, recursive=FALSE, ignore.case = FALSE)

    base_cov = data.frame("ID" = data_pts$ID)


    for (file in files){
      if (!grepl("STD", file)){
        #get variable name and raster
        var_name = sub(pattern=paste("_", city_short, ".*", sep = ""), "", x=file)
        var_name = sub(pattern=paste(".*", "Landsat_", sep = ""), "", x=var_name)
        assign(x=var_name, value = raster(paste0(Dir_input_var, "BASE/", city, "/", file)))
        raster_variable = get(var_name)

        projected_raster = projectRaster(raster_variable, crs = CRS(paste("+init=epsg:", myepsg, sep="")))

        #extract in buffers
        base_cov_extracted = buffer %>%
          st_drop_geometry() %>%
          cbind(exact_extract(projected_raster, buffer, fun = 'mean', progress = T))

        base_cov = join(base_cov, base_cov_extracted, by = "ID")
        names(base_cov)[length(names(base_cov))] <- var_name

      }
    }

    write.table(base_cov, file=paste(cov_path, city_short, "_base_cov", save_suffix, ".csv", sep=""), sep=";", dec=",", row.names=F)

    # -------------------------------
    #  CLIMATIC
    # -------------------------------

    clim_cov = data.frame("ID" = data_pts_4326$ID)
    clim_dir = paste0(Dir_input_var, "CLIMATIC/", city, '/')
    files = list.files(clim_dir, recursive = TRUE, full.names = F, pattern = ".tif$")

    for (file in files){

      if (startsWith(file, "LST_Day") || startsWith(file, "LST_Night")){

        if (grepl("Day", file, fixed=TRUE)) {
          timing = "LST_Day"
        } else {
          timing = "LST_Night"
        }

        r = stack(paste0(clim_dir, file))

        cov = data_pts_4326 %>%
          st_drop_geometry() %>%
          cbind(raster::extract(r,data_pts_4326)) #res is 1km

        var_name =  sub(pattern=".*Celsius_mean_mean_", "", x=names(r))

        colnames(cov) = c("ID", paste(timing, "_", var_name, sep=""))


      } else {
        var_name = sub(pattern="_2005.*|_1970.*", "", x=file)

        r = raster(paste0(clim_dir, file)) %>%
          projectRaster(crs=CRS(paste("+init=epsg:", myepsg, sep="")))

        cov = cbind(st_drop_geometry(data_pts), exact_extract(r, buffer, fun = 'mean', progress = T)) #res is 100m

        colnames(cov) = c("ID", var_name)

      }

      clim_cov = join(clim_cov, cov, by ="ID")

    }

    write.table(clim_cov, file=paste(cov_path, city_short, "_clim_cov", save_suffix, ".csv", sep=""), sep=";", dec=",", row.names = F)

    # -------------------------------
    # MERGE ALL COVARIATES
    # -------------------------------
    #merge all covariates in one DF

    files = list.files(path = cov_path, pattern="*.csv", full.names= T, recursive=FALSE, ignore.case = FALSE)
    full_cov = data.frame("ID" = data_pts$ID)

    for (file in files) {
      if (! grepl("full_cov", file) & grepl(paste0("cov", save_suffix, ".csv"), file)){

        cov_df = read.csv(file, sep = ";")
        full_cov = join(full_cov, cov_df, by = "ID")

      }
    }

    write.table(full_cov, file=paste0(cov_path, city_short, "_full_cov", save_suffix, ".csv"), sep=";", dec=",", row.names=F)

  }
}

