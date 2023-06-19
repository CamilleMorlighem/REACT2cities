#----------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# This function simulates a spatial displacement on non-DHS surveys, following the displacement procedure of DHS.
#
# It creates
# - A csv file and a shapefile with displaced surveys (for each city)
# - A shapefile with buffers of 1 km around each survey data point (for each city)
# - A training and testing sets on the original non-displaced data
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities              List with one vector c(city name, city short name, city epsg of projected crs) per city
#' @param crs_REACT           Geographic coordinate reference system
#' @param malaria_data_db     Path to malaria database
#' @param DirOut              Directory for malaria outputs
#' @param Max_up_age          Maximum age of survey participants (integer)
#' @param survey_period       Vector with minimum and maximum survey dates (in years)
#' @param extents             Vector with path to Shapefiles whose extent define the city extent
#' @param admin_files         Path to admin level 2 administrative boundaries shapefile

#'
#' @return                    Nothing
#' @export                    # Per city: A csv and shp for displaced surveys
#'                                        A shp with 1km buffers around displaced data point
#--------------------------------------------------------------------------------------------------------------------------


simulate.displacement = function(cities, crs_REACT, malaria_data_db, DirOut, Max_up_age, survey_period, extents, admin_files){

  # -------------------------------
  #  IMPORT AND PROCESS RAW DB
  # -------------------------------

  xlsx_wb<-loadWorkbook(malaria_data_db)
  malaria_data <- readWorksheet(xlsx_wb, sheet = "REACT_data")

  ### _Selecting the relevant variables and rows ###
  malaria_data <- malaria_data[, c("No",
                                   "ID",
                                   "COUNTRY",
                                   "CITY",
                                   "FULL_NAME",
                                   "LAT",
                                   "LONG",
                                   "MM",
                                   "YY",
                                   "LO_AGE",
                                   "UP_AGE",
                                   "NUM_EX",
                                   "NUM_Pf",
                                   "PfPR2_10",
                                   "SPATIAL_RELIABILITY",
                                   "SAMPLE_ORIGIN",
                                   "SEASON")]

  ### Selection of malaria data upon criteria in paper ###
  malaria_data = malaria_data %>%
    filter(CITY %in%  c(cities[[1]][1])) %>%
    filter(YY<=survey_period[[2]] & YY>=survey_period[[1]]) %>%
    filter(UP_AGE <= Max_up_age)

  malaria_data <- malaria_data[malaria_data$SPATIAL_RELIABILITY!="DHS",]

  ### _[Assigning short names] ###
  colnames(malaria_data) <-      c("No",
                                   "ID",
                                   "COUNTRY",
                                   "CITY",
                                   "NAME",
                                   "LAT",
                                   "LONG",
                                   "MM",
                                   "YY",
                                   "LO_AGE",
                                   "UP_AGE",
                                   "NUM_EX",
                                   "NUM_Pf",
                                   "PfPR2_10",
                                   "SPAT_REL",
                                   "SAMPLE_O",
                                   "SEASON")


  # -------------------------------
  #  SIMULATE DISPLACEMENT
  # -------------------------------

  displaceData <- function(dataset.sf.unflag, admin_shp) {

    processed_df = data.frame()

    na_count = dataset.sf.unflag %>%
      dplyr::select(flag) %>%
      is.na() %>%
      sum()

    while (na_count > 0) {

      dataset.sf = dataset.sf.unflag %>%
        filter(is.na(flag))

      # 2. Generate a random direction by generating angle between 0 and 360, and converting the angle from degrees to radians.
      dataset.sf$angle_rad <- deg2rad(runif(nrow(dataset.sf), 0, 360))

      # 3. Generate a random distance in meters of 0-5,000 meters for Rural points with 1% of rural points being given 0-10,000 meter distance.
      # # get number of rural clusters
      # nRuralClusters <- sum(dataset.sf$URBAN_RURA == "R")
      # #   Split urban and rural
      # dataset_R.sf <- subset(dataset.sf, URBAN_RURA == "R")
      # dataset_U.sf <- subset(dataset.sf, URBAN_RURA == "U")
      #   For rural, assign 5000 meter displacement with 1% randomly assigned up to 10000 meter displacement
      # dataset_R.sf$m_displaced <- ifelse({runif(nRuralClusters) |> rank(ties.method = "random") <= floor(nRuralClusters*0.01)}, 10000, 5000)
      #   For urban, assign 2000 meter displacement
      dataset.sf$m_displaced <- 2000
      # #   Combine them back together
      # dataset.sf <- rbind(dataset_R.sf, dataset_U.sf)

      dataset.sf$random_meters <- runif(nrow(dataset.sf), min = 0, max = dataset.sf$m_displaced)

      # 4. Generate the offset by applying trigonometry formulas (law of cosines) using the distance as the hypotenuse and the radians calculated in step 2.
      dataset.sf$xOffset <- sin(dataset.sf$angle_rad)*dataset.sf$random_meters
      dataset.sf$yOffset <- cos(dataset.sf$angle_rad)*dataset.sf$random_meters

      # 5. Add the offset to the original coordinate (in meters) to return the displaced coordinates.
      dataset.sf$newX <- st_coordinates(dataset.sf)[,1] + dataset.sf$xOffset
      dataset.sf$newY <- st_coordinates(dataset.sf)[,2] + dataset.sf$yOffset

      # Remove geometry
      dataset_C <- st_set_geometry(dataset.sf, NULL)
      dataset_f.sf <- st_as_sf(dataset_C, coords = c("newX", "newY"), crs = st_crs(admin_shp))

      # Check if points are still in same admin unit level 2
      check_admin = st_intersection(dataset_f.sf, admin_shp) %>%
        filter(orig_admin_ID == admin_ID)

      dataset.sf.unflag = dataset.sf.unflag %>%
        filter(ID %in% check_admin$ID) %>%
        mutate(flag = "flagged") %>%
        rbind(filter(dataset.sf.unflag, !ID %in% check_admin$ID))

      na_count = dataset.sf.unflag %>%
        dplyr::select(flag) %>%
        is.na() %>%
        sum()

      processed_df = dataset.sf %>%
        st_drop_geometry() %>%
        filter(ID %in% check_admin$ID) %>%
        #dplyr::select(ID, newX, newY, xOffset, yOffset) %>%
        rbind(processed_df)

    }

    processed_df <- st_as_sf(processed_df, coords = c("newX", "newY"), crs = st_crs(admin_shp))
    return(processed_df)
  }

  malaria_data_sf = st_as_sf(malaria_data,
                             coords = c("LONG", "LAT"),
                             crs = crs_REACT)

  lapply(cities, FUN = function(city_lst){

    city = city_lst[1]
    city_short = city_lst[2]
    myepsg = city_lst[3]
    city_i = match(list(city_lst), cities)

    # convert to proj. crs
    malar_sf_4326 = malaria_data_sf %>%
      filter(CITY == city)

    admin = admin_files[city_i] %>%
      st_read() %>%
      st_transform(as.numeric(myepsg)) %>%
      rownames_to_column("admin_ID") %>%
      dplyr::select(admin_ID)

    malar_sf = malar_sf_4326 %>%
      st_transform(as.numeric(myepsg)) %>%
      st_intersection(admin) %>%
      dplyr::rename(orig_admin_ID = admin_ID) %>%
      mutate(flag = NA)

    # displace points
    malar_sf_displ_4326 = displaceData(malar_sf, admin) %>%
      st_transform(crs_REACT) %>%
      dplyr::select(-c(flag, orig_admin_ID))

    # # visualize results
    # plot(malar_sf_displ$geometry)
    # plot(malar_sf$geometry, add = T, col = "red")

    # -------------------------------
    #  CREATE SHP AND CSV FILES
    # -------------------------------

    city_DirOut = paste0(DirOut, as.character(city), '/simul_displ')
    dir.create(file.path(city_DirOut), showWarnings = F, recursive = T)
    city_name = paste("malaria_data_", as.character(city_short), sep="")

    ### Selection of data points in smaller covariate extent (LU)
    extent = st_read(extents[city_i]) %>%
      st_make_valid() %>%
      st_buffer(0) %>%
      st_union() %>%
      st_transform(crs_REACT) %>%
      st_make_valid()

    df = list(malar_sf_4326, malar_sf_displ_4326)

    i = 1

    for (m_sf in df){

      in_extent = st_intersects(m_sf, extent, sparse=F)
      m_sf = m_sf[in_extent,] ; nrow(m_sf)

      if (i == 2){
        city_name = paste0(city_name, "_displ")
      }

      ### _Build one csv file ###
      write.table(m_sf, file = paste0(city_DirOut, "/", city_name, ".csv"), quote=F, row.names=F, sep=";", dec=",")

      ### _Building the shapefile ###
      st_write(m_sf,
               dsn=city_DirOut,
               layer=city_name,
               layer_options = "ENCODING=UTF-8",
               driver="ESRI Shapefile",
               delete_layer = TRUE)

      ### Building one km buffer ###
      m_sf %>%
        st_transform(as.numeric(myepsg)) %>%
        st_buffer(1000) %>%
        st_write(dsn=city_DirOut,
                 layer=paste0(city_name, "_1km_round"),
                 layer_options = "ENCODING=UTF-8",
                 driver="ESRI Shapefile",
                 delete_layer = TRUE)

      i = i + 1

    }

    # -------------------------------
    #  Generate test-training sets
    # -------------------------------

    # Load shapefile of points with malaria prevalence, as 'sf' object
    malaria_data = malar_sf_4326 %>%
      mutate(ID = as.numeric(ID),
             PfPR2_10 = as.numeric(PfPR2_10)) %>%
      dplyr::select(-geometry)

    resampled = malaria_data %>%
      sample_frac(0.77)

    write.csv2(resampled, paste0(city_DirOut, "/malaria_data_", city_short, "_training_set.csv"))

    test = malaria_data %>%
      filter(!ID %in% resampled$ID) ; nrow(test)

    write.csv2(test, paste0(city_DirOut, "/malaria_data_", city_short, '_test_set.csv'))

  })

}
