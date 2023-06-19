#----------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# This function loads the input malaria database and allows to select data based on three criteria:
# - Surveys that are or not DHS
# - Surveys that include or not adults (<= 16 years old)
# - Surveys from a certain time period
#
# This function creates
# - A csv file and a shapefile with selected surveys (for each city)
# - A shapefile with buffers of 1 km around each survey data point (for each city)
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities              List with one vector c(city name, city short name, city epsg of projected crs) per city
#' @param crs_REACT           Geographic coordinate reference system
#' @param malaria_data_db     Path to malaria database
#' @param DirOut              Directory for malaria outputs
#' @param Max_up_age          Maximum age of survey participants (integer)
#' @param survey_period       Vector with minimum and maximum survey dates (in years)
#' @param use_DHS             Boolean. If T use DHS data; if F don't use DHS. F is default
#' @param extents             Vector with path to Shapefiles whose extent define the city extent
#'
#' @return                    Nothing
#' @export                    # Per city: A csv and shp with selected surveys
#'                                        A shp with 1km buffers around survey data point
#--------------------------------------------------------------------------------------------------------------------------


select.malaria.data = function(cities, crs_REACT, malaria_data_db, DirOut, Max_up_age, survey_period, use_DHS=F, extents){

  # -------------------------------
  #  IMPORT AND PROCESS RAW DB
  # -------------------------------

  malaria_data <- malaria_data_db %>%
    loadWorkbook() %>%
    readWorksheet(sheet = "REACT_data") %>%
    select(c("No", "ID", "COUNTRY", "CITY", "FULL_NAME", "LAT", "LONG", "YY",
             "LO_AGE","UP_AGE", "NUM_EX", "NUM_Pf", "PfPR2_10", "SPATIAL_RELIABILITY",
             "SAMPLE_ORIGIN", "SEASON"))

  ### _[Assigning short names] ###
  colnames(malaria_data) <-      c("No",
                                   "ID",
                                   "COUNTRY",
                                   "CITY",
                                   "NAME",
                                   "LAT",
                                   "LONG",
                                   "YY",
                                   "LO_AGE",
                                   "UP_AGE",
                                   "NUM_EX",
                                   "NUM_Pf",
                                   "PfPR2_10",
                                   "SPAT_REL",
                                   "SAMPLE_O",
                                   "SEASON")

  ### Selection of malaria data upon criteria in paper ###
  malaria_data = malaria_data %>%
    filter(CITY %in%  c(cities[[1]][1], cities[[2]][1], cities[[3]][1], cities[[4]][1])) %>%
    filter(YY<=survey_period[[2]] & YY>=survey_period[[1]]) %>%
    filter(UP_AGE <= Max_up_age)

  if (use_DHS == F){
    malaria_data <- malaria_data %>%
      filter(SPAT_REL!="DHS")
  }

  malaria_data_sf = st_as_sf(malaria_data,
                             coords = c("LONG", "LAT"),
                             crs = crs_REACT)

  # -------------------------------
  #  CREATE SHP AND CSV FILES
  # -------------------------------

  lapply(cities, FUN = function(city_lst){

    city = city_lst[1]
    city_short = city_lst[2]
    myepsg = city_lst[3]
    city_i = match(list(city_lst), cities)

    # convert to proj. crs
    malar_sf_4326 = malaria_data_sf %>%
      filter(CITY == city)

    city_DirOut = paste0(DirOut, as.character(city))
    dir.create(file.path(city_DirOut), showWarnings = F, recursive = T)
    city_name = paste("malaria_data_", as.character(city_short), sep="")

    ### Selection of data points in smaller covariate extent (LU)
    extent = st_read(extents[city_i]) %>%
      st_make_valid() %>%
      st_buffer(0) %>%
      st_union() %>%
      st_transform(crs_REACT) %>%
      st_make_valid()

    in_extent = st_intersects(malar_sf_4326, extent, sparse=F)
    malar_sf_4326 =  malar_sf_4326[in_extent,]

    ### _Build one csv file ###
    write.table(malar_sf_4326, file = paste(city_DirOut, "/", city_name, ".csv", sep=""), quote=F, row.names=F, sep=";", dec=",")

    ### _Building the shapefile ###
    st_write(malar_sf_4326,
             dsn=city_DirOut,
             layer= city_name,
              #layer_options = "ENCODING=UTF-8",
             driver="ESRI Shapefile",
             delete_layer = TRUE)

    ### Building one km buffer ###
    malar_sf_4326 %>%
      st_transform(as.numeric(myepsg)) %>%
      st_buffer(1000) %>%
      st_write(dsn=city_DirOut,
               layer=paste0(city_name, "_1km_round"),
               layer_options = "ENCODING=UTF-8",
               driver="ESRI Shapefile",
               delete_layer = TRUE)

  })
}
