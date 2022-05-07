#----------------------------------------------------------------------------------------------------------
# Multi-satellite environmental and socio-economic predictors of vector-borne diseases in African cities: 
# malaria as an example
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
    
    cov_df = read.csv(paste(input_cov, city, "/", city_short, "_full_cov.csv", sep = "") , sep = ";", dec=",")
    malaria_data_sf=st_read(paste(input_malaria, city, '/', city_name, '.shp', sep=""))
    cov_df$ID = as.numeric(cov_df$ID) ; malaria_data_sf$ID = as.numeric(malaria_data_sf$ID)
    full_cov = dplyr::left_join(malaria_data_sf[,c("ID", "PfPR2_10")], cov_df, by=c("ID"= "ID"))
    
    #remove NaN values from df 
    covariates = na.omit(full_cov)
    covariates = subset(covariates, select = -c(ID))
    
    covariates = st_drop_geometry(covariates)
    
    dir.create(file.path(paste(input_cov,  '/', city , "/Correlation/", sep="")), showWarnings = F)
    
    # -------------------------------
    #  CORRELATION  
    # -------------------------------
    #check covariance
    #cov(covariates)
    
    #check correlation 
    #spearman checks whether there is a correlation (does not have to be linear)
    #pearson checks for LINEAR correlation 
    #pearson correlated ==> spearman correlated but reverse is not true 
    #SPEARMAN 
    par(mar=c(7,4,8,2))
    
    
    M = stats::cor(covariates, use = "all", method = "spearman")
    try(corrplot(M, method = 'square', order = 'FPC', type = 'lower', diag=F, tl.cex = 0.5))
    try(corrplot(M, method = 'number', order = 'FPC', type = 'lower', diag=F))
    
    pdf(file = paste0(input_cov,  '/', city ,'/Correlation/', city_short, "_spearman_corr.pdf", sep=""))
    p_values <- cor.mtest(covariates, conf.level = 0.95)
    try(corrplot(M, p.mat = p_values$p, 'circle', sig.level = 0.05,  insig = "blank", tl.cex = 0.5))
    dev.off()
    
    #PEARSON 
    M2 = stats::cor(covariates, use = "all", method = "pearson")
    try(corrplot(M2, method = 'square', order = 'FPC', type = 'lower'))
    try(corrplot(M2, method = 'number', order = 'FPC', type = 'lower'))
    
    pdf(file = paste0(input_cov,  '/', city ,'/Correlation/', city_short, "_pearson_corr.pdf", sep=""))
    p_values <- cor.mtest(covariates, conf.level = 0.95)
    try(corrplot(M2, p.mat = p_values$p, 'circle', sig.level = 0.05,  insig = "blank", tl.cex = 0.5))
    dev.off()
    
    par(mar = c(5, 4, 4, 2) + 0.1) # default
    
    # -------------------------------
    #  SCATTERPLOTS
    # -------------------------------
    dir.create(file.path(paste(input_cov,  '/', city , "/Correlation/Scatterplot/", sep="")), showWarnings = F)
    
    i = 1
    for (i in 1:(ncol(covariates)-1)){
      par(mar = c(5,6,4,4))
      pdf(file = paste(input_cov,  '/', city , "/Correlation/Scatterplot/" ,city_short,"_plot_PfPR_",colnames(covariates[i+1]),".pdf", sep=""))
      plot(y = covariates$PfPR2_10, x = covariates[,i+1], xlab = paste0(colnames(covariates[i+1])), ylab = "PfPR2_10")
      abline(lm(covariates$PfPR2_10 ~ covariates[,i+1]), col = "grey")
      model <- lm(covariates$PfPR2_10 ~ covariates[,i+1])
      
      dev.off()
    }
    
    
    pdf(file = paste0(input_cov,  '/', city , "/Correlation/Scatterplot/" ,city_short,"_PfPR_allCov.pdf"))
    par(mar = c(5,6,4,4), mfrow = c(3,3))
    for (i in 1:(ncol(covariates)-1)){
      
      plot(y = covariates$PfPR2_10, x = covariates[,i+1], xlab = paste0(colnames(covariates[i+1])), ylab = "PfPR2_10")
      abline(lm(covariates$PfPR2_10 ~ covariates[,i+1]), col = "grey")
      
    } ; 
    dev.off()
    
    par(mar = c(5, 4, 4, 2) + 0.1, mfrow = c(1,1))
    
    pdf(file = paste0(input_cov,  '/', city , "/Correlation/Scatterplot/" ,city_short,"_OtherVar.pdf"))
    par(mar = c(5,6,4,4), mfrow = c(3,3))
    
    plot(y = covariates$LCZ_compact, x = covariates$LU_planned, xlab = "LU_Planned", ylab = "LCZ_compact")
    abline(lm(covariates$LCZ_compact ~ covariates$LU_planned), col = "grey")
    
    plot(y=covariates$LC_trees, x=covariates$LC_bare_ground)
    abline(lm(covariates$LC_trees ~ covariates$LC_bare_ground), col = "grey")
    
    plot(y=covariates$LU_ACS, x=covariates$LC_buildings)
    abline(lm(covariates$LU_ACS ~ covariates$LC_buildings), col = "grey")
    
    plot(y=covariates$LC_bare_ground, x=covariates$LC_trees)
    abline(lm(covariates$LC_bare_ground ~ covariates$LC_trees), col = "grey")
    
    plot(y=covariates$LC_buildings, x=covariates$LC_trees)
    abline(lm(covariates$LCZ_compact ~ covariates$LC_trees), col = "grey")
    
    
    if (city_short == "Kamp" || city_short == "DES"){
      plot(y=covariates$AVG_TS, x=covariates$AVG_RH2M)
      abline(lm(covariates$AVG_TS ~ covariates$AVG_RH2M), col = "grey")
      
      plot(y=covariates$AVG_TS, x=covariates$LU_informal)
      abline(lm(covariates$AVG_TS ~ covariates$LU_informal), col = "grey")
      
      plot(y=covariates$LC_trees, x=covariates$MAX_PREC)
      abline(lm(covariates$LC_trees ~ covariates$MAX_PREC), col = "grey")
      
      plot(y=covariates$LC_buildings, x=covariates$AVG_T2M)
      abline(lm(covariates$LC_buildings ~ covariates$AVG_T2M), col = "grey")
      
      
      plot(y=covariates$LC_buildings, x=covariates$AVG_TS)
      abline(lm(covariates$LC_buildings ~ covariates$AVG_TS), col = "grey")
      
    }
    
    dev.off()
    
  }
  
}




