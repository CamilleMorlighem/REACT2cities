#----------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
#
# This function creates a grid to be used for predicting Pfpr2-10. More specifically it creates:
# - a 1km resolution grid stored as a raster
# - a 1km resolution grid stored as a vector (shp)
# - a point shp with the centroids of the grid cells
# - a csv that stores the index of the cells and their xy coordinates
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities                 List with one vector c(city name, city short name, city epsg of projected crs) per city
#' @param crs_REACT              Geographic coordinate reference system
#' @param Dir_input_var          Path to input geospatial variables
#' @param Dir_input_malaria_data Path to input malaria data-points
#' @param Dir_output_data        Path to output prediction grid
#
#' @return                       Nothing
#' @export                       # Per city: A 1km resolution grid stored as a raster
#                                            A 1km resolution grid stored as a vector (shp)
#                                            A point shp with the centroids of the grid cells
#                                            A csv that stores the index of the cells and their xy coordinates
#--------------------------------------------------------------------------------------------------------------------------


prediction.grid = function(cities, crs_REACT, Dir_input_var, Dir_input_malaria_data, Dir_output_data){

  for (city_vect in cities){

    city = city_vect[1]
    city_short = city_vect[2]

    if(!dir.exists(file.path(paste0(Dir_output_data, city, "/")))){ dir.create(file.path(paste0(Dir_output_data, city, "/")), recursive = TRUE)   }

    # -------------------------------
    #  Create 1km raster grid
    # -------------------------------
    smaller_extent = raster(paste0(Dir_output_data, city, "/GRID_1KM_", city_short))

    r = raster(ext = extent(smaller_extent), resolution =1000, crs = crs(smaller_extent))
    ref_grid = projectRaster(smaller_extent, r, crs = crs_REACT, method = "ngb")

    #Define index values
    idx <- Which(!is.na(ref_grid), cells=TRUE)
    raster::values(ref_grid)[idx] <- 1
    writeRaster(ref_grid, paste0(Dir_output_data, city, "/", city_short, "_1km_grid"), format="GTiff", overwrite=TRUE)

    # -------------------------------
    #  Convert raster to polygons
    # -------------------------------
    poly_ref <- spex::polygonize(ref_grid)
    st_crs(poly_ref) = crs(smaller_extent)
    poly_ref <- subset(poly_ref, select = - 1)
    poly_ref = tibble::rownames_to_column(poly_ref, "idx")
    st_agr(poly_ref) = "constant" # for avoind warning message "assumes attributes are constant over geometries of x"

    st_write(poly_ref, dsn=paste0(Dir_output_data,city, "/", city_short, "_1km_vect"), driver = "ESRI Shapefile",  delete_layer = TRUE)

    # -------------------------------
    #  Extract ref points
    # -------------------------------
    pt_ref <- st_centroid(poly_ref)
    st_write(pt_ref, dsn=paste0(Dir_output_data,city, "/", city_short, "_1km_pts"), driver = "ESRI Shapefile", delete_layer = TRUE)

    # -------------------------------
    #  Get ref table
    # -------------------------------
    table_ref <- data.frame(st_coordinates(pt_ref))
    table_ref$idx <- idx # !
    table_ref <- table_ref[,c("idx", "X","Y")]

    write.csv(table_ref,
              paste0(Dir_output_data, city, "/", city_short, "_1km_ref_table.csv"),
              row.names=FALSE)

  }

}
