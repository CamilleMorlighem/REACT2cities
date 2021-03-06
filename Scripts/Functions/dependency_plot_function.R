# ------------------------------------------------------------------------------------------------------------------
# Modelling of the intra-urban malaria risk using socioeconomic and environmental factors with a random forest (RF) model
# in two sub-Saharan cities: Kampala (Uganda) and Dar es Salaam (Tanzania)
# 
# This function creates dependency plots for covariates used in a given RF model 
#--------------------------------------------------------------------------------------------------------------------------
#' @param RF_mlr_output     RF model object 
#' @param variables         A dataframe containing all variables of the model (dependent variables + predictor variables)
#' @param predictor_dataset A dataframe containing all predictor variables of the model 
#' @param output_file       A chr defining the output file path 
#' @param output_file_name  A chr defining the output file name 
#' @param HYt_mode          An int defining how to perform the hyperparameter tuning 
#'
#' @return                  The dependency plot values stored into a list 
#' @export                  # A DP plot for each covariate in the RF model (pdf)
# ------------------------------------------------------------------------------------------------------------------


dependency.plot = function(RF_mlr_output, predictor_dataset, variables, output_file, output_file_name, HYt_mode){

  valuesDepPlot_tot = list()
  result = RF_mlr_output
  predictor_list = colnames(predictor_dataset) 
  
  # loop to store the values for the dependency plot for each of the covariate
  for (var_nb in 1:length(predictor_list)){
    
    var_name = predictor_list[var_nb]
    
    # get x values of the predictor 
    # ---------------------------------------------------------------------------------------------------------------------------------
    if(HYt_mode == "HYt_outsideCV"){
      x_values = pdp::partial(result$models[[1]]$learner.model, pred.var = var_name, plot = F, train = variables)[,var_name]
    } 
    if(HYt_mode == "HYt_insideCV"){ 
      x_values = pdp::partial(result$models[[1]]$learner.model$next.model$learner.model, pred.var = var_name, plot = F, train = variables)[,var_name] # with wrapped_lrn_rf
    }
    
    # get y values of the response variable 
    # ---------------------------------------------------------------------------------------------------------------------------------
    if(HYt_mode == "HYt_outsideCV"){
      y_values = lapply(result$models, function(x) { 
        pdp::partial(x$learner.model, plot = F, pred.var = var_name, train = variables)$yhat } )
    } 
    if(HYt_mode == "HYt_insideCV"){
      y_values = lapply(result$models, function(x) { 
      pdp::partial(x$learner.model$next.model$learner.model, plot = F, pred.var = var_name, train = variables)$yhat }  ) # with wrapped_lrn_rf
    }
    y_values = as.data.frame(do.call(cbind,y_values))
    
    # combine x and y values into a dataframe 
    # ---------------------------------------------------------------------------------------------------------------------------------
    valuesDepPlot =cbind(x_values, y_values)
    colnames(valuesDepPlot) =c(paste0(var_name,"_x"), paste0(paste0(paste0(var_name),"_yhat_"),1:length(result$models)))
    
    # start pdf plot 
    # ---------------------------------------------------------------------------------------------------------------------------------
    pdf(paste0(output_file,"DepPlot_",output_file_name,"_",var_name,".pdf"))
    
    # Define plot margins  
    par(mar=c(5,5,3,2))
    
    # Define ylim and xlim of the plot 
    plot(y = "", x = "", xlab = "", ylab="", xlim=range(x_values), ylim=range(y_values), 
         main = paste(var_name))
    
    # create rug to get the distribution of PfPR210 data on the plot  
    dec = quantile(predictor_dataset[,var_name], prob = seq(0, 1, length = 11), type = 5)
    dec = dec[-c(1,length(dec))]
    rug(x=dec, col="grey", lwd="1")
    
    # add a dependence plot line for each model simulation 
    for (p in 2:ncol(valuesDepPlot)){
      lines(x=valuesDepPlot[,1], y=valuesDepPlot[,p], col=rgb(0,0,0,0,alpha = 0.09))
    }
    
    # add a lines for the median of all simulations & median+sd & median-sd 
    summary_dpt = apply(valuesDepPlot[,2:length(result$models)],1,function(x) { cbind(mean(x), sd(x)) } )
    lines(x = valuesDepPlot[,1], y = summary_dpt[1,])
    lines(x = valuesDepPlot[,1], y = (summary_dpt[1,]+summary_dpt[2,]) , lty = 2)
    lines(x = valuesDepPlot[,1], y = (summary_dpt[1,]-summary_dpt[2,]) , lty = 2)
    
    # Stop plot - create final pdf 
    dev.off() 
    
    # store the values in a list 
    valuesDepPlot_tot = c(valuesDepPlot_tot,list(valuesDepPlot))
    
  }

return (valuesDepPlot_tot)

}