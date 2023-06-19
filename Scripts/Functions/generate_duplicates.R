# ------------------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# This code generates duplicates to DHS data-points in the fourth cardinal directions at 500 m and 1 km, which can
# later be used to implement spatial optimisation methods
#
# This function creates
# - A csv file and a shapefile with duplicated data-points
# - A shapefile with buffers of 1 km around each survey data point
#--------------------------------------------------------------------------------------------------------------------------
#' @param cities                   List with one vector c(city name, city short name, city epsg of projected crs) per city
#' @param crs_REACT                Geographic coordinate reference system
#' @param Dir_input_malaria_data   Path to malaria database
#' @param admin_files              Path to admin level 2 administrative boundaries shapefile
#' @param simulation               Boolean indicates whether simulation study is conducted (Default F)

#' @return                    Nothing
#' @export                    # Per city: A csv and shp with duplicates
#'                                        A shp with 1km buffers around duplicates
#--------------------------------------------------------------------------------------------------------------------------


generate.duplicates = function(cities, crs_REACT, Dir_input_malaria_data, admin_files, simulation = F){

    for (city_vect in cities){

        city = city_vect[1]
        city_short = city_vect[2]
        myepsg = city_vect[3]
        city_i = match(list(city_vect), cities)

        subdir = '/'
        save_suffix = ""

        if (simulation) {

            subdir = 'simul_displ/'
            save_suffix = "_displ"

        }

        city_name = paste0("malaria_data_", city_short)

        dhs_data = paste0(Dir_input_malaria_data, city, '/', subdir, city_name, save_suffix, '.csv') %>%
            read.csv2() %>%
            filter(UP_AGE <= 16) %>%
            dplyr::select(-geometry)

        dhs_pts = paste0(Dir_input_malaria_data, city, '/', subdir, city_name, save_suffix, '.shp') %>%
            st_read() %>%
            # transform to proj crs
            st_transform(crs = as.numeric(myepsg)) %>%
            dplyr::select(ID) %>%
            mutate(ID = as.numeric(ID)) %>%
            right_join(dhs_data)

        # get admin id
        admin = admin_files[city_i] %>%
            st_read() %>%
            st_transform(as.numeric(myepsg)) %>%
            rownames_to_column("admin_ID") %>%
            dplyr::select(admin_ID)

       dhs_pts = dhs_pts %>%
            st_intersection(admin) %>%
            dplyr::rename(orig_admin_ID = admin_ID)

        if (! simulation) {

            #select DHS survey points
            dhs_pts = dhs_pts %>%
                filter(SPAT_REL == "DHS")

        }

        #create displaced directory
        out_displaced_dir = paste0(Dir_input_malaria_data, city, '/', subdir, "duplicates")
        dir.create(out_displaced_dir, showWarnings = F)

        #convert spatial df to df with x,y columns
        separated_coord <- dhs_pts %>%
            mutate(northing = unlist(map(dhs_pts$geometry,2)),
                   easting = unlist(map(dhs_pts$geometry,1))) %>%
            dplyr::select(ID, northing, easting, PfPR2_10, orig_admin_ID) %>%
            st_drop_geometry()

        #create duplicates in the North, South, West and East directions for two radii: 500 and 1000 m
        radius = c(500,1000)
        x = separated_coord$easting
        y = separated_coord$northing
        PfPR2_10 = as.numeric(separated_coord$PfPR2_10)
        or_ID = as.numeric(separated_coord$ID)
        duplicates = data.frame()

        for (rad in radius){
            xPlus <- separated_coord$easting+rad
            xMinus <- separated_coord$easting-rad
            yPlus <- separated_coord$northing+rad
            yMinus <- separated_coord$northing-rad

            up = data.frame(cbind("PfPR2_10" = PfPR2_10 , "orig_admin_ID" =  separated_coord$orig_admin_ID ,"Orig_ID" = or_ID, "x" = x, "y"= yPlus))
            down = data.frame(cbind("PfPR2_10" = PfPR2_10 , "orig_admin_ID" =  separated_coord$orig_admin_ID ,"Orig_ID" = or_ID, "x" = x, "y" = yMinus))
            east = data.frame(cbind("PfPR2_10" = PfPR2_10 ,"orig_admin_ID" =  separated_coord$orig_admin_ID , "Orig_ID" = or_ID, "x" = xPlus, y))
            west = data.frame(cbind("PfPR2_10" = PfPR2_10 , "orig_admin_ID" =  separated_coord$orig_admin_ID ,"Orig_ID" = or_ID, "x" = xMinus, y))

            duplicates = rbind(duplicates, up, down, east, west)
        }

        #convert df with x,y columns to spatial df
        duplicates = st_as_sf(duplicates, coords = c('x', 'y'), crs = as.numeric(myepsg))

        #set ID
        duplicates <- tibble::rownames_to_column(duplicates, "ID") %>%
            mutate(ID = as.numeric(ID))

        #check admin area in which duplicates lie
        duplicates = duplicates %>%
            st_intersection(admin) %>%
            filter(orig_admin_ID == admin_ID) %>%
            dplyr::select(-admin_ID, -orig_admin_ID)

        #add original DHS points & non-DHS points to spatial df
        duplicates = dhs_pts %>%
            dplyr::select(ID, PfPR2_10, geometry) %>%
            mutate(Orig_ID = as.numeric(ID)) %>%
            rbind(duplicates)

        #transform to geographic CRS
        duplicates_4326 <- duplicates %>%
            st_transform(crs = crs_REACT) %>%
            mutate(ID = as.numeric(ID),
                   PfPR2_10 = as.numeric(PfPR2_10)) %>%
            dplyr::select(ID, Orig_ID, PfPR2_10)

        #save to csv and shp
        duplicates_4326 %>%
            write.csv2(file = paste(out_displaced_dir, '/', city_name, ".csv", sep=""), quote=F, row.names=F, sep=";")

        st_write(duplicates_4326,
                 dsn=out_displaced_dir,
                 layer=city_name,
                 driver="ESRI Shapefile",
                 delete_layer = TRUE)

        #ROUND BUFFER
        duplicates %>%
            st_buffer(1000) %>%
            st_write(dsn=out_displaced_dir,
                     layer=paste0(city_name, "_1km_round"),
                     layer_options = "ENCODING=UTF-8",
                     driver="ESRI Shapefile",
                     delete_layer = TRUE)

    }
}
