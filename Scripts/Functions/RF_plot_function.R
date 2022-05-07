# ------------------------------------------------------------------------------------------------------------------
# Modelling of the intra-urban malaria risk using socioeconomic and environmental factors with a random forest (RF) model
# in two sub-Saharan cities: Kampala (Uganda) and Dar es Salaam (Tanzania)
#  
# This function plots the evaluation metrics (R2, NMAE and NRMSE) of given RF models in horizontal or vertical mode  
#--------------------------------------------------------------------------------------------------------------------------
#' @param RFResult          Evaluation metric values of the given set of RF models built in spatial cross-validation 
#' @param stat_plot_list    Vector defining the statistic used to order the RF models in the boxplot (at most "Mean", "Median", "byFactor"). byFactor means models are ordered alphabetically by RF models name. 
#' @param factor            Either "nb_test" if used to compare RFE iterations or "GeoSpDt_modtype" if used to compare group of covariates 
#' @param output_file       A chr defining the output file path 
#' @param output_file_name  A chr defining the output file name 
 
#' @return                  Nothing 
#' @export                  # Boxplot (pdf) to compare RFE iterations or group of covariates 
#--------------------------------------------------------------------------------------------------------------------------


RF.plot.result = function(RFResult, factor, stat_plot_list, output_file, output_file_name){
 
  
  RFResult$factor = RFResult[,which(colnames(RFResult) == paste0(factor) )]
  
   
  GoF_list = c("R.squared", "MSE", "MAE", "RMSE")
  
  for (o in 1:length(GoF_list)){
    
    # o = 1
    GoF = GoF_list[o] 
    
    
    for (stat_i in 1:length(stat_plot_list)){
      
      # stat_i = 1
      my_stat = stat_plot_list[stat_i]
      
      
      if (my_stat == "Mean"){
        
        if(GoF == "R.squared"){
          p<-ggplot(RFResult, aes(x = reorder(factor,R2_mean), y=R.squared))+labs(y = "R.squared", x = "") }
        if(GoF == "MSE"){
          p<-ggplot(RFResult, aes(reorder(factor,MSE_mean), y=MSE))+labs(y = "MSE", x = "") }
        if(GoF == "MAE"){
          p<-ggplot(RFResult, aes(reorder(factor,MAE_mean), y=MAE))+labs(y = "MAE", x = "") }
        if(GoF == "RMSE"){
          p<-ggplot(RFResult, aes(reorder(factor,RMSE_mean), y=RMSE))+labs(y = "RMSE", x = "") }
             
      }
      
      
      if(my_stat=="Median"){
        
        if(GoF == "R.squared"){
          p<-ggplot(RFResult, aes(x = reorder(factor,R2_median), y=R.squared))+labs(y = "R.squared", x = "") }
        if(GoF == "MSE"){
          p<-ggplot(RFResult, aes(reorder(factor,MSE_median), y=MSE))+labs(y = "MSE", x = "") }
        if(GoF == "MAE"){
          p<-ggplot(RFResult, aes(reorder(factor,MAE_median), y=MAE))+labs(y = "MAE", x = "") }
        if(GoF == "RMSE"){
          p<-ggplot(RFResult, aes(reorder(factor,RMSE_median), y=RMSE))+labs(y = "RMSE", x = "") }
             
      }
      
      
      if(my_stat=="byFactor"){
        
        RFResult$nb_test = factor(RFResult$nb_test, levels = c(paste0("RFE-",1:length(levels(RFResult$nb_test)) )))
        
        GoF = GoF_list[o]
        if(GoF == "R.squared"){
          p<-ggplot(RFResult, aes(factor(nb_test), y=R.squared))+labs(y = "R.squared", x = "") }
        if(GoF == "MSE"){
          p<-ggplot(RFResult, aes(factor(nb_test), y=MSE))+labs(y = "MSE", x = "") }
        if(GoF == "MAE"){
          p<-ggplot(RFResult, aes(factor(nb_test), y=MAE))+labs(y = "MAE", x = "") }
        if(GoF == "RMSE"){
          p<-ggplot(RFResult, aes(factor(nb_test), y=RMSE))+labs(y = "RMSE", x = "") }
      
        }
        
    
      p<-p +
        geom_violin(adjust = 2, width=0.5)+
        geom_boxplot(width=0.1)+
        #geom_jitter(shape=16, position=position_jitter(0.2))+
        # geom_bar(stat="identity", color="black", fill="grey", 
        #          position=position_dodge()) +
        #reorder(CovImp_summary_sort,var)
        # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
        #               position=position_dodge(.9)) +
        # 
        #scale_fill_grey() + # # Use grey scale
        # scale_fill_manual(values=c("#E69F00", "#999999")) +
        #theme_classic()+
        ggtitle(paste0(RFResult$city[1]))+
        theme(axis.title.y = element_text(size = rel(1.8), margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
              axis.text.y = element_text(size = rel(1.8), margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
              
              axis.ticks.x = element_blank(),
              axis.title.x = element_text(size = rel(1.8), margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)), 
              axis.text.x = element_text(size = rel(1.8), margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0), 
                                         angle = 90),
              plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 0.5, l = 1.5, "cm"), 
              plot.title = element_text(size=17,hjust = 0.5,  margin = ggplot2::margin(t = 0, r = 0, b = 30, l = 0) ), 
              #legend.position="none",
              panel.background = element_rect(fill = "transparent") # bg of the panel
              , plot.background = element_rect(fill = "transparent") # bg of the plot
              , panel.grid.major = element_blank() # get rid of major grid
              #, panel.grid.minor = element_blank() # get rid of minor grid
              , panel.grid.minor = element_line(color = "grey", linetype = 3)
              , panel.grid.major.y = element_line(color = "grey", linetype = 3)
        )
      # p
      
      ggsave(paste0("Boxplot_",output_file_name,"_orderBy",my_stat,"_",GoF,".png"),
             height=8, width=12,
             units="in", dpi=350, plot= p, bg = "transparent",
             path = paste(output_file))
    }   
    
  }
}
  













