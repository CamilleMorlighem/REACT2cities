# ------------------------------------------------------------------------------------------------------------------
#' Random forest model, cov importance & partial dependence plot computation: 
#' 
# - Run a random forest (RF) model in spatial cross validation with hyperparameter tuning 
# - store RF results, model characteristics and goodness-of-fit indices in a csv file     
# - Extract covariates importance and store these values store in a csv file
# - Create a barplot figure with the importance of each covariate - in pdf
# - Create a dependence plot for each covariate of the model - in pdf 
# ------------------------------------------------------------------------------------------------------------------
#' @param variables         A dataframe containing all variables of the model (dependent variables + predictor variables)
#' @param coordinates       A dataframe containing the coordinates of the dependent variable data 
#' @param dependent_var     A chr defining the name of the column containing the dependent variable
#' @param var_group         A chr defining the group of variables being assessed: "all3Geo", "all_predictors", "GeoSpDt"
#' @param city              A chr defining the study area : here a city among 'DES' or 'Kampala'
#' @param mode              A chr defining the modelling mode: "ALL" or "RFE"
#' @param nb_test           A chr defining the step in RFE variable selection
#' @param folds_rf          A int defining the number of folds (see mlr) in the spatial cross-validation 
#' @param reps_rf           A int defining the number of repetitions in the spatial cross-validation   
#' @param iters_rf          A int defining the number of iterations in the spatial cross-validation 
#' @param nb_cores          A int defining the number of cores which can be used for the paralellisation with the computer 
#' @param output_file       A chr defining the output file path 
#' @param output_file_name  A chr defining the output file name 
#' @param save_result       A chr defining whether results have to be saved or not: "yes" or "not
#' @param HYt_mode          An int defining how to perform the hyperparameter tuning 
#' @param save_RF_object    A chr defining whether RF object has to be saved or not: "yes" or "not
#'
#' @return                  A dataframe with a summary of the covariate importance 
#' @export                  # RF model results and GoF csv file, Covariate importance csv file, Barplot of cov imp pdf file and dependence plot pdf file
# ------------------------------------------------------------------------------------------------------------------


RF.CovImp.pdp = function(variables, coordinates, dependent_var,
                         var_group, city, mode, nb_test,
                         folds_rf, reps_rf, iters_rf, nb_cores, 
                         output_file, output_file_name, save_result, HYt_mode, save_RF_object = c("yes")){

  predictor_dataset = dplyr::select(variables,-dependent_var)
  
  # ------------------------------------------------------------------------------------------------------------------
  #       Spatial Cross-Validation + hyperparameters (HP) tuning with mlr R package & Ranger package 
  # ------------------------------------------------------------------------------------------------------------------
  
  # create task
  task = makeRegrTask(data=variables, target="PfPR2_10", coordinates=coordinates)
  
  # learner
  lrn_rf = makeLearner(cl="regr.ranger", predict.type = "response", importance = "permutation")  # Ranger package 
  
  # Performance estimation level - set with a spatial cross-validation and 5 folds
  perf_level = makeResampleDesc(method = "SpRepCV", folds = folds_rf, reps = reps_rf) #5*2*50
  
  # use 50 randomly selected hyperparameters
  ctrl = makeTuneControlRandom(maxit=50L) 
  
  # Define the outer limits of the randomly selected hyperparameters
  ps=makeParamSet(
    makeIntegerParam("mtry", lower=1, upper=(ncol(predictor_dataset)-1)),
    makeNumericParam("sample.fraction", lower = 0.2, upper = 0.9), # Pas fait par Camille
    makeIntegerParam("min.node.size", lower=1, upper=10) )
  
  # five spatially disjoint partitions
  tune_level = makeResampleDesc("SpCV", iters = iters_rf)
  
  # Modify the learner in accordance with all the characteristics defining the hyperparameter tuning
  wrapped_lrn_rf = makeTuneWrapper(learner = lrn_rf, 
                                   resampling = tune_level,
                                   par.set = ps,
                                   control = ctrl, 
                                   show.info = TRUE,
                                   measures = mlr::rmse)
  
  # Start paralellisation 
  parallelStartSocket(nb_cores)
  
  result = mlr::resample(learner = wrapped_lrn_rf, 
                         task = task,
                         resampling = perf_level,
                         extract = getTuneResult,
                         models = TRUE,
                         measures = list(mlr::rmse, mlr::mae, mlr::mse)) #GOF wanted
  # stop parallelization
  parallelStop()

  
  if(save_RF_object == "yes"){
    saveRDS(result, file = paste0(output_file, "RFoutput_",output_file_name,".rds"))
  }
  
  
  # --------------------------------------------------------------------------------------
  # Get performance results and summary results 
  # --------------------------------------------------------------------------------------
  
  # create table of goodness of fit (GoF) indicators and modelling parameters 
  tmp=lapply(result$models, function(i){ cbind(i$learner.model$next.model$learner.model$treetype, "PfPR2_10",
                                               i$learner.model$next.model$learner.model$mtry,
                                               i$learner.model$next.model$learner.model$min.node.size,
                                               i$learner.model$next.model$learner.model$num.trees, 
                                               i$learner.model$next.model$learner.model$prediction.error, 
                                               i$learner.model$next.model$learner.model$r.squared) }) 
  RF_result = as.data.frame(do.call(rbind,tmp))
  colnames(RF_result) =c("model", "dependent_var", 
                         "mtry", "nodesize", "ntree", 
                         "MSE", "R-squared")
  RF_result$RMSE = result$measures.test$rmse
  RF_result$MAE = result$measures.test$mae
  RF_result$MSE = result$measures.test$mse
  RF_result$sample_fraction = unlist(lapply(result$extract, function(i){c(i$x$sample.fraction)}))
  RF_result$covariates = paste(result$models[[1]]$learner.model$next.model$features, collapse = ",")
  RF_result$SpCV = perf_level$id
  RF_result$SpCV_details = paste("folds = ", folds_rf,", reps = ", reps_rf)
  RF_result$ctrl_HY = ctrl$extra.args$maxit
  if(var_group!= "BEST" & mode == 'RFE'){RF_result$nb_test = nb_test}
  RF_result$n = nrow(predictor_dataset)
  RF_result$n_cov = length(colnames(predictor_dataset))
  RF_result$var_group = paste(var_group)
  RF_result$city = paste(city)

  # transform column in numeric 
  RF_result$MSE = as.numeric(as.character(RF_result$MSE))
  RF_result$`R-squared` = as.numeric(as.character(RF_result$`R-squared`))
  RF_result$ntree = as.numeric(as.character(RF_result$ntree))
  RF_result$nodesize = as.numeric(as.character(RF_result$nodesize))
  RF_result$mtry = as.numeric(as.character(RF_result$mtry))
  
  
  if(save_result == "yes"){ write.csv2(RF_result, file = paste0(output_file,"RFResults_",output_file_name,".csv")) }
  
  
  # --------------------------------------------------------------------------------------
  # Importance of cov
  # --------------------------------------------------------------------------------------
  
  CovImp = lapply(result$models, function(x) {x$learner.model$next.model$learner.model$variable.importance})
  CovImp  = as.data.frame(do.call(rbind,CovImp))
  
  if(save_result == "yes"){  write.csv2(CovImp, file = paste0(output_file,"RFCovImp_",output_file_name,".csv")) }
  
  CovImp_summary = lapply(CovImp, function(x) { cbind(mean(x), sd(x)) } )
  names_tmp = names(CovImp_summary)
  CovImp_summary = as.data.frame(do.call(rbind,CovImp_summary))
  colnames(CovImp_summary) = c("mean","sd")
  rownames(CovImp_summary) = names_tmp
  CovImp_summary$var = names_tmp
  
  
  # Plot
  p<- ggplot(CovImp_summary, aes(x = reorder(var, -mean), y=mean)) + 
    geom_bar(stat="identity", color="black", fill="grey", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_grey() + # Use grey scale
    labs(y = "IncMSE", x = "")+
    ggtitle(paste0(city,"\n",var_group))+
    theme(axis.title.y = element_text(size = rel(1.8), margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
          axis.text.y = element_text(size = rel(1.8), margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
          axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = rel(1.8), margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)), 
          axis.text.x = element_text(size = rel(1.8), margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), 
                                     angle = 90),
          plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 0.5, l = 1.5, "cm"), 
          plot.title = element_text(size=17,hjust = 0.5,  margin = ggplot2::margin(t = 0, r = 0, b = 30, l = 0) ), 
          panel.background = element_rect(fill = "transparent") # bg of the panel
          , plot.background = element_rect(fill = "transparent") # bg of the plot
          , panel.grid.major = element_blank() # get rid of major grid
          , panel.grid.minor = element_line(color = "grey", linetype = 3)
          , panel.grid.major.y = element_line(color = "grey", linetype = 3)
    ) #p
  
  if(save_result == "yes"){
    ggsave(paste0("Barplot_CovImp_",output_file_name,".tiff"), 
           height=8, width=12, 
           units="in", dpi=300, plot= p, bg = "transparent",
           path = paste(output_file)) ; rm(p)          }
  
  
  # --------------------------------------------------------------------------------------
  # Partial dependence plot - 1 plot/country
  # --------------------------------------------------------------------------------------
  
  if (save_result == "yes"){
    
    # create dependence plot pdf and store values 
    valuesDepPlot_tot = dependency.plot(RF_mlr_output = result, predictor_dataset = predictor_dataset, variables = variables, 
                                        output_file = output_file ,output_file_name = output_file_name, HYt_mode = HYt_mode)
    
    # Save partial dependence plot values in rds file
    if(save_result == "yes"){
      names(valuesDepPlot_tot) = colnames(predictor_dataset)
      saveRDS(valuesDepPlot_tot, file = paste0(output_file, "DepPlotValues_",output_file_name,".rds"))  }
    
  }
  
  return (CovImp_summary)
}
