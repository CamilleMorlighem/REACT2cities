#----------------------------------------------------------------------------------------------------------
# Malaria risk mapping in sub-Saharan African cities using environmental and socio-economic predictors
#
# This code does a descriptive analysis of the input covariates
#
# It creates
# - Correlation matrices (pearson and spearman) between the covariates
# - Scatterplots showing the pfpr2-10 in y and covariates in x
#----------------------------------------------------------------------------------------------------------
#' @param cities              List with one vector c(city name, city short name, city epsg of projected crs) per city
#' @param crs_REACT           Geographic coordinate reference system
#' @param input_cov           Path to input covariates
#' @param input_malaria       Path to input malaria data-points

#' @return                    Nothing
#' @export                    # Per city: Correlation matrices (pearson and spearman) between the covariates (pdf)
#'                                        Scatterplots showing the pfpr2-10 in y and covariates in x (pdf)
#--------------------------------------------------------------------------------------------------------------------------


check.correlation = function(cities, crs_REACT, input_cov, input_malaria) {

  for (city_vect in cities){

    city = city_vect[1]
    city_short = city_vect[2]
    city_name = paste("malaria_data_", as.character(city_short), sep="")

    cov_df = paste0(input_cov, city, "/", city_short, "_full_cov.csv") %>%
      read.csv(sep = ";", dec=",") %>%
      mutate(ID = as.numeric(ID))

    malaria_data_sf = paste0(input_malaria, city, '/', city_name, '.shp') %>%
      st_read() %>%
      mutate(ID = as.numeric(ID)) %>%
      dplyr::select(ID, PfPR2_10)

    full_cov = dplyr::left_join(malaria_data_sf, cov_df, by=c("ID"= "ID"))

    #remove NaN values from df
    covariates = full_cov %>%
      na.omit %>%
      subset(select = -c(ID)) %>%
      st_drop_geometry()

    dir.create(file.path(paste(input_cov,  '/', city , "/Correlation/", sep="")), showWarnings = F)

    # -------------------------------
    #  CORRELATION
    # -------------------------------
    #spearman checks whether there is a correlation (does not have to be linear)
    #pearson checks for LINEAR correlation
    #pearson correlated ==> spearman correlated but reverse is not true
    #SPEARMAN
    par(mar=c(7,4,8,2))

    M = stats::cor(covariates, use = "all", method = "spearman")
    pdf(file = paste0(input_cov,  '/', city ,'/Correlation/', city_short, "_spearman_corr.pdf", sep=""))
    p_values <- cor.mtest(covariates, conf.level = 0.95)
    try(corrplot(M, p.mat = p_values$p, 'circle', sig.level = 0.05,  insig = "blank", tl.cex = 0.5))
    dev.off()

    #PEARSON
    M2 = stats::cor(covariates, use = "all", method = "pearson")
    pdf(file = paste0(input_cov,  '/', city ,'/Correlation/', city_short, "_pearson_corr.pdf", sep=""))
    p_values <- cor.mtest(covariates, conf.level = 0.95)
    try(corrplot(M2, p.mat = p_values$p, 'circle', sig.level = 0.05,  insig = "blank", tl.cex = 0.5))
    dev.off()

    # -------------------------------
    #  SCATTERPLOTS
    # -------------------------------
    dir.create(file.path(paste(input_cov,  '/', city , "/Correlation/Scatterplot/", sep="")), showWarnings = F)

    for (i in 1:(ncol(covariates)-1)){

      tryCatch({

      par(mar = c(5,6,4,4))
      pdf(file = paste0(input_cov,  '/', city , "/Correlation/Scatterplot/" ,city_short,"_plot_PfPR_",colnames(covariates[i+1]),".pdf"))
      plot(y = covariates$PfPR2_10, x = covariates[,i+1], xlab = paste0(colnames(covariates[i+1])), ylab = "PfPR2_10")
      abline(lm(covariates$PfPR2_10 ~ covariates[,i+1]), col = "grey")
      model <- lm(covariates$PfPR2_10 ~ covariates[,i+1])

      dev.off()

      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    }


    pdf(file = paste0(input_cov,  '/', city , "/Correlation/Scatterplot/" ,city_short,"_PfPR_allCov.pdf"))
    par(mar = c(5,6,4,4), mfrow = c(3,3))

    for (i in 1:(ncol(covariates)-1)){

      tryCatch({

      plot(y = covariates$PfPR2_10, x = covariates[,i+1], xlab = paste0(colnames(covariates[i+1])), ylab = "PfPR2_10")
      abline(lm(covariates$PfPR2_10 ~ covariates[,i+1]), col = "grey")

      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    }

    dev.off()

  }

}




