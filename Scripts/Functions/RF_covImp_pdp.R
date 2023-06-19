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
#' @param city              A chr defining the city
#' @param city_short        A chr defining the city short name
#' @param mode              A chr defining the modelling mode: "ALL" or "RFE"
#' @param nb_test           A chr defining the step in RFE variable selection
#' @param folds_rf          A int defining the number of folds (see mlr) in the spatial cross-validation
#' @param reps_rf           A int defining the number of repetitions in the spatial cross-validation
#' @param iters_rf          A int defining the number of iterations in the spatial cross-validation
#' @param nb_cores          A int defining the number of cores which can be used for the paralellisation with the computer
#' @param Dir_output_rf     A chr defining the output file path
#' @param save_result       A chr defining whether results have to be saved or not: "yes" or "not
#' @param save_RF_object    A chr defining whether RF object has to be saved or not: "yes" or "not
#'
#' @return                  A dataframe with a summary of the covariate importance
#' @export                  # RF model results and GoF csv file, Covariate importance csv file, Barplot of cov imp pdf file and dependence plot pdf file
# ------------------------------------------------------------------------------------------------------------------


RF.CovImp.pdp = function(variables, coordinates, city,
                         city_short, mode, nb_test,
                         folds_rf, reps_rf, iters_rf, nb_cores,
                         Dir_output_rf, save_RF_object, save_result){

  predictor_dataset = dplyr::select(variables,-"PfPR2_10")

  if (mode == "RFE"){
    Dir_output_rf = paste0(Dir_output_rf, 'RFE_', nb_test, "/")
    dir.create(file.path(Dir_output_rf), showWarnings = F)

  }

  # ------------------------------------------------------------------------------------------
  #  1. Spatial Cross-Validation + hyperparameters (HP) tuning with mlr & Ranger packages
  # ------------------------------------------------------------------------------------------

  # create task
  task = makeRegrTask(data=variables, target="PfPR2_10", coordinates=coordinates)

  # learner
  lrn_rf = makeLearner(cl="regr.ranger", predict.type = "response", importance = "permutation")  # Ranger package

  # Performance estimation level - set with a spatial cross-validation and 5 folds
  perf_level = makeResampleDesc(method = "SpRepCV", folds = folds_rf, reps = reps_rf)

  # use 50 randomly selected hyperparameters
  ctrl = makeTuneControlRandom(maxit=50L)

  # Define the outer limits of the randomly selected hyperparameters
  ps=makeParamSet(
    makeIntegerParam("mtry", lower=1, upper=(ncol(predictor_dataset)-1)),
    makeNumericParam("sample.fraction", lower = 0.2, upper = 0.9),# 0.2 et 0.9
    makeIntegerParam("min.node.size", lower=1, upper=10))

  # five spatially disjoint partitions
  tune_level = makeResampleDesc("SpCV", iters = iters_rf)

  # Modify the learner in accordance with all the characteristics defining the hyperparameter tuning
  wrapped_lrn_rf = makeTuneWrapper(learner = lrn_rf,
                                   resampling = tune_level,
                                   par.set = ps,
                                   control = ctrl,
                                   show.info = TRUE,
                                   measures = mlr::rmse)

  registerDoSEQ()
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)

  # Start paralellisation
  parallelStartSocket(nb_cores)

  result = mlr::resample(learner =wrapped_lrn_rf,
                         task = task,
                         resampling = perf_level,
                         extract = getTuneResult,
                         models = TRUE,
                         measures = list(mlr::rmse, mlr::mae, mlr::mse, mlr::rsq)) #GOF wanted
  # stop parallelization
  parallelStop()

  if(save_RF_object == "yes"){
    saveRDS(result, file = paste0(Dir_output_rf, "/RF_output_", city_short,".rds"))
  }

  # --------------------------------------------------------------------------------------
  #  2. Get performance results and summary results
  # --------------------------------------------------------------------------------------

  #create table of goodness of fit (GoF) indicators and modelling parameters
  tmp=lapply(result$models, function(i){ cbind(i$learner.model$next.model$learner.model$treetype, "PfPR2_10",
                                               i$learner.model$next.model$learner.model$mtry,
                                               i$learner.model$next.model$learner.model$min.node.size,
                                               i$learner.model$next.model$learner.model$num.trees,
                                               i$learner.model$next.model$learner.model$prediction.error ,
                                               i$learner.model$next.model$learner.model$r.squared)
    })

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

  if(mode == 'RFE'){RF_result$nb_test = nb_test}

  RF_result$n = nrow(predictor_dataset)
  RF_result$n_cov = length(colnames(predictor_dataset))
  RF_result$city = paste(city)

  # transform column in numeric
  RF_result$MSE = as.numeric(as.character(RF_result$MSE))
  RF_result$`R-squared` = as.numeric(as.character(RF_result$`R-squared`))
  RF_result$ntree = as.numeric(as.character(RF_result$ntree))
  RF_result$nodesize = as.numeric(as.character(RF_result$nodesize))
  RF_result$mtry = as.numeric(as.character(RF_result$mtry))

  if(save_result == "yes"){
    write.csv2(RF_result, file = paste0(Dir_output_rf,"/RF_Results_", city_short,".csv"))
  }


  # --------------------------------------------------------------------------------------
  # 3. Importance of cov
  # --------------------------------------------------------------------------------------

  #var imp for each cov for each model (permutation importance)
  CovImp = lapply(result$models, function(x) {x$learner.model$next.model$learner.model$variable.importance})
  CovImp  = as.data.frame(do.call(rbind,CovImp))

  if(save_result == "yes"){write.csv2(CovImp, file = paste0(Dir_output_rf,"/RF_VarImp_", city_short, ".csv")) }

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
    ggtitle(paste0(city,"\n"))+
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
    )

  if(save_result == "yes"){
    ggsave(paste0(Dir_output_rf, "/RF_VarImp_", city_short,".tiff"),
           height=8, width=12,
           units="in", dpi=300, plot= p,
           bg = "transparent") ; rm(p)          }



  # --------------------------------------------------------------------------------------
  # 4. Partial dependence plot - 1 plot/country
  # --------------------------------------------------------------------------------------

  if (mode=="All" || nb_test == 1){

    #create the cluster
    try(my.cluster <- parallel::makeCluster(
      nb_cores/2,
      type = "PSOCK"))

    #register it to be used by %dopar%
    try(doParallel::registerDoParallel(cl = my.cluster))

  }

  if (save_result == "yes"){

    # create dependence plot pdf and store values
    valuesDepPlot_tot = dependency.plot(RF_mlr_output = result, predictor_dataset = predictor_dataset, variables = variables,
                                        output_file = Dir_output_rf, output_file_name = city_short, nb_cores)

    # Save partial dependence plot values in rds file
    if(save_result == "yes"){
      names(valuesDepPlot_tot) = colnames(predictor_dataset)
      saveRDS(valuesDepPlot_tot, file = paste0(Dir_output_rf, "/RF_DepPlot_", city_short,".rds"))  }

  }


  return (CovImp_summary)


}


get_eval_metrics = function(results, nb_test){

  mean_R2 = mean(results$R.squared) ; mean_MAE = mean(results$MAE) ; mean_MSE = mean(results$RMSE)
  min_R2 = min(results$R.squared) ; min_MAE = min(results$MAE) ; min_MSE = min(results$RMSE)
  max_R2 = max(results$R.squared) ; max_MAE = max(results$MAE) ; max_MSE = max(results$RMSE)

  df = data.frame("nb_test" = nb_test,
                  "mean_R2" = mean_R2, "min_R2" = min_R2, 'max_R2' = max_R2,
                  "mean_MAE" = mean_MAE, "min_MAE" = min_MAE, "max_MAE" = max_MAE,
                  "mean_RMSE" = mean_MSE, "min_RMSE" = min_MSE, "max_RMSE" = max_MSE)

  return(df)
}
