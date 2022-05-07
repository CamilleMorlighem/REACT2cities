#----------------------------------------------------------------------------------------------------------
# Multi-satellite environmental and socio-economic predictors of vector-borne diseases in African cities: 
# malaria as an example
#
# This function loads the input malaria database and selects data based on three criteria: 
# - Surveys which are not DHS 
# - Surveys which do not include adults (<= 16 years old)
# - Surveys between 2005 and 2016
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
  malaria_data <- malaria_data[malaria_data$CITY %in% c(cities[[1]][1], cities[[2]][1]),]
  malaria_data <- malaria_data[malaria_data$YY<=survey_period[[2]] & malaria_data$YY>=survey_period[[1]],]
  malaria_data <- malaria_data[malaria_data$UP_AGE <= Max_up_age,]
  if (use_DHS == F){
    malaria_data <- malaria_data[malaria_data$SPATIAL_RELIABILITY!="DHS",]
  }

  
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
  
  city_i = 1
  
  for (city_vect in cities){
    
    city = city_vect[1]
    city_short = city_vect[2]
    city_DirOut = paste(DirOut, as.character(city), sep="")
    dir.create(file.path(city_DirOut), showWarnings = F, recursive = T)
    
    
    # -------------------------------
    #  CREATE SHP AND CSV FILES
    # -------------------------------
    
    #create spatial df
    city_name = paste("malaria_data_", as.character(city_short), sep="")
    assign(x = city_name, value=malaria_data[malaria_data$CITY==city,])
    city_df = eval(parse(text = city_name))
    malaria_data_sf = st_as_sf(city_df,
                               coords = c("LONG", "LAT"),
                               crs = crs_REACT)
    
    
    
    ### Selection of data points in smaller covariate extent (LU)
    extent = st_read(extents[city_i])
    extent=st_union(extent)
    extent = st_transform(extent, crs_REACT)
    in_extent=st_intersects(malaria_data_sf, extent, sparse=F)
    malaria_data_sf=malaria_data_sf[in_extent,] ; nrow(malaria_data_sf)
    
    ### _Build one csv file ###
    write.table(city_df, file = paste(city_DirOut, "/", city_name, ".csv", sep=""), quote=F, row.names=F, sep=";", dec=",")
    
    
    ### _Building the shapefile ###
    st_write(malaria_data_sf,
             dsn=city_DirOut,
             layer=city_name,
             layer_options = "ENCODING=UTF-8",
             #encoding="UTF-8",
             driver="ESRI Shapefile",
             #delete_dsn = TRUE,
             delete_layer = TRUE)
    
    
    # -------------------------------
    #  CREATE ONE KM BUFFERS
    # -------------------------------
    
    myepsg = city_vect[3]
    malaria_data_sf_meter <- st_transform(malaria_data_sf, crs = as.numeric(myepsg)) # requires first a conversion from decimal degreses to meters
    
    #ROUND BUFFER
    malaria_data_sf_1K <- st_buffer(malaria_data_sf_meter, 1000) # Set the 1k buffer
    
    st_write(malaria_data_sf_1K,
             dsn=city_DirOut,
             layer=paste(city_name, "_1km_round", sep=""),
             layer_options = "ENCODING=UTF-8",
             driver="ESRI Shapefile",
             delete_layer = TRUE)
    
    city_i = city_i + 1
  }
}
