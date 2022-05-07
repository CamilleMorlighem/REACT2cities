# ------------------------------------------------------------------------------------------------------------------
# Modelling of the intra-urban malaria risk using socioeconomic and environmental factors with a random forest (RF) model
# in two sub-Saharan cities: Kampala (Uganda) and Dar es Salaam (Tanzania)
#  
# This function plots the evaluation metrics (R2, NMAE and NRMSE) of given RF models in horizontal or vertical mode  
#--------------------------------------------------------------------------------------------------------------------------
#' @param RFResult          Evaluation metric values of the given set of RF models built in spatial cross-validation 
#' @param stat_plot_list    Vector defining the statistic used to order the RF models in the boxplot (at most "Mean", "Median", "byFactor"). byFactor means models are ordered alphabetically by RF models name. 
#' @param orientation       Plot orientation ("vertical" or "horizontal")
#' @param GoF_list          Vector with evaluation metrics to be considered (at most "NMAE", "NRMSE", "R.squared", "MSE", "MAE")
#' @param output_file       A chr defining the output file path 
#' @param output_file_name  A chr defining the output file name 
#' @param rangeNRMSE        Range of NRMSE values  
#' @param rangeNMAE         Range of NMAE values 
#' @param rangeR2           Range of R2 values 
#' @param city              A chr defining the city of interest 
#' @param plot_mode         Either "ALL" (boxplots with models based on all covariates), "RFE" (boxplots based on RFE cov) or "ALL_RFE" (boxplots based on both)

#' @return                  Nothing 
#' @export                  # Boxplot (pdf) to compare different RF models
#--------------------------------------------------------------------------------------------------------------------------


RF.plot.comp.result = function(RFResult, stat_plot_list, orientation, GoF_list, output_file, output_file_name, 
                               rangeNRMSE, rangeNMAE, rangeR2, city, plot_mode){
 

  factor = "plot_name"
  RFResult$factor = RFResult[,which(colnames(RFResult) == paste0(factor) )]
  
  if (plot_mode == "RFE"){
    RFResult[RFResult$factor=="Base (ALL)",]$factor = "Base"
    RFResult[RFResult$factor=="LCZdp_LULC (RFE)",]$factor = "Base + LCZ + LULC"
    RFResult[RFResult$factor=="LCZdp_LULC_COSMO (RFE)",]$factor = "Base + LCZ + COSMO + LULC"
    RFResult[RFResult$factor=="LCZdp_COSMO (RFE)",]$factor = "Base + LCZ + COSMO"
    RFResult[RFResult$factor=="LCZdp (RFE)",]$factor = "Base + LCZ"
    
  } else if (plot_mode == "ALL_RFE"){
    RFResult[RFResult$factor=="Base (ALL)",]$factor = "Base (ALL)"
    RFResult[RFResult$factor=="LCZdp_LULC (RFE)",]$factor = "Base + LCZ + LULC (RFE)"
    RFResult[RFResult$factor=="LCZdp_LULC_COSMO (RFE)",]$factor = "Base + LCZ + COSMO + LULC (RFE)"
    RFResult[RFResult$factor=="LCZdp_COSMO (RFE)",]$factor = "Base + LCZ + COSMO (RFE)"
    RFResult[RFResult$factor=="LCZdp (RFE)",]$factor = "Base + LCZ (RFE)"
    RFResult[RFResult$factor=="LCZdp_LULC (ALL)",]$factor = "Base + LCZ + LULC (ALL)"
    RFResult[RFResult$factor=="LCZdp_LULC_COSMO (ALL)",]$factor = "Base + LCZ + COSMO + LULC (ALL)"
    RFResult[RFResult$factor=="LCZdp_COSMO (ALL)",]$factor = "Base + LCZ + COSMO (ALL)"
    RFResult[RFResult$factor=="LCZdp (ALL)",]$factor = "Base + LCZ (ALL)"
    
  }
  
  for (o in 1:length(GoF_list)){
    
    GoF = GoF_list[o] 
    
    for (stat_i in 1:length(stat_plot_list)){
      
      my_stat = stat_plot_list[stat_i]
      
      if (my_stat == "Mean"){
        
        if(GoF == "R.squared"){
          #p<-ggplot(RFResult, aes(x = reorder(factor,R2_mean), y=R.squared, fill = mode))+labs(y = "R.squared", x = "") }
          p<-ggplot(RFResult, aes(x = reorder(factor,R2_mean), y=R.squared))+labs(y = "R.squared", x = "") }
      
        if(GoF == "MSE"){
          #p<-ggplot(RFResult, aes(reorder(factor,MSE_mean), y=MSE, fill = mode))+labs(y = "MSE", x = "") }
          p<-ggplot(RFResult, aes(reorder(factor,MSE_mean), y=MSE))+labs(y = "MSE", x = "") }
        
        if(GoF == "MAE"){
          #p<-ggplot(RFResult, aes(reorder(factor,MAE_mean), y=MAE, fill = mode))+labs(y = "MAE", x = "") }
          p<-ggplot(RFResult, aes(reorder(factor,MAE_mean), y=MAE))+labs(y = "MAE", x = "") }
        
        if(GoF == "RMSE"){
          #p<-ggplot(RFResult, aes(reorder(factor,RMSE_mean), y=RMSE, fill = mode))+labs(y = "RMSE", x = "") }
          p<-ggplot(RFResult, aes(reorder(factor,RMSE_mean), y=RMSE))+labs(y = "RMSE", x = "") }
      
      }
      
      if(my_stat=="Median"){
        
        if(GoF == "R.squared"){
          
          #p<-ggplot(RFResult, aes(x = reorder(factor,R2_median), y=R.squared, fill = mode))+labs(y = expression('R'^2), x = "")+ylim(rangeR2)
          p<-ggplot(RFResult, aes(x = reorder(factor,R2_median), y=R.squared))+labs(y = expression('R'^2), x = "")+ylim(rangeR2)
          
          p1 = p
          }
        if(GoF == "MSE"){
          #p<-ggplot(RFResult, aes(reorder(factor,MSE_median), y=MSE, fill = mode))+labs(y = "MSE", x = "") }
          p<-ggplot(RFResult, aes(reorder(factor,MSE_median), y=MSE))+labs(y = "MSE", x = "") }
        
        if(GoF == "MAE"){
          #p<-ggplot(RFResult, aes(reorder(factor,MAE_median), y=MAE, fill = mode))+labs(y = "MAE", x = "") }
          p<-ggplot(RFResult, aes(reorder(factor,MAE_median), y=MAE))+labs(y = "MAE", x = "") }
    
        if(GoF == "RMSE"){
          
          #p<-ggplot(RFResult, aes(reorder(factor,RMSE_median), y=RMSE, fill = mode))+labs(y = "RMSE", x = "") }
          p<-ggplot(RFResult, aes(reorder(factor,RMSE_median), y=RMSE))+labs(y = "RMSE", x = "") }
  
        if(GoF == "NRMSE"){
          
          #p<-ggplot(RFResult,aes(reorder(factor,NRMSE_median), y=NRMSE, fill = mode))+labs(y = "NRMSE", x = "")+ylim(rangeNRMSE[1],rangeNRMSE[2])
          p<-ggplot(RFResult,aes(reorder(factor,NRMSE_median), y=NRMSE))+labs(y = "NRMSE", x = "")+ylim(rangeNRMSE[1],rangeNRMSE[2]) 
          
          p2 = p
           }
        if(GoF == "NMAE"){
      
          #p<-ggplot(RFResult, aes(reorder(factor,NMAE_median), y=NMAE, fill = mode))+labs(y = "NMAE", x = "")+ylim(rangeNMAE)
          p<-ggplot(RFResult, aes(reorder(factor,NMAE_median), y=NMAE))+labs(y = "NMAE", x = "")+ylim(rangeNMAE) 
        
          p3 = p
          
          }
       }
      
      if(my_stat=="byFactor"){
        
        RFResult$nb_test = factor(RFResult$nb_test, levels = c(paste0("RFE-",1:length(levels(RFResult$nb_test)) )))
        
        GoF = GoF_list[o]
        if(GoF == "R.squared"){
          #p<-ggplot(RFResult, aes(factor(nb_test), y=R.squared))+labs(y = "R.squared", x = "") }
        p<-ggplot(RFResult, aes(factor(nb_test), y=R.squared))+labs(y = "R.squared", x = "") }
      
        if(GoF == "MSE"){
          #p<-ggplot(RFResult, aes(factor(nb_test), y=MSE))+labs(y = "MSE", x = "") }
        p<-ggplot(RFResult, aes(factor(nb_test), y=MSE))+labs(y = "MSE", x = "") }
      
        if(GoF == "MAE"){
          #p<-ggplot(RFResult, aes(factor(nb_test), y=MAE))+labs(y = "MAE", x = "") }
        p<-ggplot(RFResult, aes(factor(nb_test), y=MAE))+labs(y = "MAE", x = "") }
      
        if(GoF == "RMSE"){
          #p<-ggplot(RFResult, aes(factor(nb_test), y=RMSE))+labs(y = "RMSE", x = "") }
        p<-ggplot(RFResult, aes(factor(nb_test), y=RMSE))+labs(y = "RMSE", x = "") }
      
      }
    }
  }

  
  # Plot - for vertical or horizontal boxplot
  # --------------------------------------------------------------------------------------
  # Vertical plot - usual one 
  if(orientation == "vertical"){
    p1<-p1 +
      geom_violin(adjust = 2, width=0.5)+
      geom_boxplot(width=0.1)+
      
      scale_fill_manual(values=c("#E69F00", "#999999")) +
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
    
    p2<-p2 +
      geom_violin(adjust = 2, width=0.5)+
      geom_boxplot(width=0.1)+
      
      scale_fill_manual(values=c("#E69F00", "#999999")) +
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
    
    p3<-p3 +
      geom_violin(adjust = 2, width=0.5)+
      geom_boxplot(width=0.1)+
      
      scale_fill_manual(values=c("#E69F00", "#999999")) +
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
    # ggsave(paste0("Boxplot_",output_file_name,"_orderBy",my_stat,"_",GoF,".tiff"),
    #        height=8, width=12,
    #        units="in", dpi=300, plot= p, bg = "transparent",
    #        path = paste(output_file))
  }    
  
  # Horizontal
  if(orientation == "horizontal"){
   
    p1 = p1 + geom_violin(adjust = 2, width=0.8)+ 
      geom_boxplot(width=0.25)+
      coord_flip() + 
      
      scale_fill_manual(values=c("#CCCCCC", "#0072B2")) + # light grey and blue 
      
      #theme_classic()+
      ggtitle(paste0(RFResult$city[1]))+
      theme(axis.title.y = element_text(size = rel(2.8), margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
            axis.text.y = element_text(size = rel(3.5), margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            #face = "bold", 
            axis.ticks.x = element_blank(),
            axis.title.x = element_text(size = rel(3.5), margin = ggplot2::margin(t = 15, r = 0, b = 0, l = 0)), 
            axis.text.x = element_text(size = rel(3.2), margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)),
            plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 0.5, l = 1.5, "cm"), 
            plot.title = element_text(size=17,hjust = 0.5,  margin = ggplot2::margin(t = 0, r = 0, b = 30, l = 0) ), 
            
            panel.background = element_rect(fill = "transparent") # bg of the panel
            , plot.background = element_rect(fill = "transparent") # bg of the plot
            #, panel.grid.major = element_blank() # get rid of major grid
            #, panel.grid.minor = element_blank() # get rid of minor grid
            , panel.grid.minor = element_line(color = "grey", linetype = 3, size = 1)
            , panel.grid.major.y = element_line(color = "grey", linetype = 3, size = 1)
            , panel.grid.major.x = element_line(color = "grey", linetype = 3, size = 1)
            
            # legend
            , legend.key = element_rect(fill = "transparent")  # Key background
            , legend.background = element_rect(fill= "transparent") # Legend background
            , legend.title = element_blank()
            #legend.position="none",
      )
    

    p2 = p2 + geom_violin(adjust = 2, width=0.8)+
      geom_boxplot(width=0.25)+
      coord_flip() +
      scale_fill_manual(values=c("#CCCCCC", "#0072B2")) + # light grey and blue 
      ggtitle(paste0(RFResult$city[1]))+
      theme(axis.title.y = element_text(size = rel(2.8), margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
            axis.text.y = element_text(size = rel(3.5), margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            #face = "bold", 
            axis.ticks.x = element_blank(),
            axis.title.x = element_text(size = rel(3.5), margin = ggplot2::margin(t = 15, r = 0, b = 0, l = 0)), 
            axis.text.x = element_text(size = rel(3.2), margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)),
            plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 0.5, l = 1.5, "cm"), 
            plot.title = element_text(size=17,hjust = 0.5,  margin = ggplot2::margin(t = 0, r = 0, b = 30, l = 0) ), 
            
            panel.background = element_rect(fill = "transparent") # bg of the panel
            , plot.background = element_rect(fill = "transparent") # bg of the plot
            #, panel.grid.major = element_blank() # get rid of major grid
            #, panel.grid.minor = element_blank() # get rid of minor grid
            , panel.grid.minor = element_line(color = "grey", linetype = 3, size = 1)
            , panel.grid.major.y = element_line(color = "grey", linetype = 3, size = 1)
            , panel.grid.major.x = element_line(color = "grey", linetype = 3, size = 1)
            
            # legend
            , legend.key = element_rect(fill = "transparent")  # Key background
            , legend.background = element_rect(fill= "transparent") # Legend background
            , legend.title = element_blank()
            #legend.position="none",
      )
    

    p3 = p3 +      
      geom_violin(adjust = 2, width=0.8)+
      geom_boxplot(width=0.25)+
      coord_flip()+
      scale_fill_manual(values=c("#CCCCCC", "#0072B2")) + # light grey and blue 
      ggtitle(paste0(RFResult$city[1]))+
      theme(axis.title.y = element_text(size = rel(2.8), margin = ggplot2::margin(t = 0, r = 20, b = 0, l = 0)), 
            axis.text.y = element_text(size = rel(3.5), margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
            #face = "bold", 
            axis.ticks.x = element_blank(),
            axis.title.x = element_text(size = rel(3.5), margin = ggplot2::margin(t = 15, r = 0, b = 0, l = 0)), 
            axis.text.x = element_text(size = rel(3.2), margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)),
            plot.margin = ggplot2::margin(t = 1, r = 0.5, b = 0.5, l = 1.5, "cm"), 
            plot.title = element_text(size=17,hjust = 0.5,  margin = ggplot2::margin(t = 0, r = 0, b = 30, l = 0) ), 
            
            panel.background = element_rect(fill = "transparent") # bg of the panel
            , plot.background = element_rect(fill = "transparent") # bg of the plot
            #, panel.grid.major = element_blank() # get rid of major grid
            #, panel.grid.minor = element_blank() # get rid of minor grid
            , panel.grid.minor = element_line(color = "grey", linetype = 3, size = 1)
            , panel.grid.major.y = element_line(color = "grey", linetype = 3, size = 1)
            , panel.grid.major.x = element_line(color = "grey", linetype = 3, size = 1)
            
            # legend
            , legend.key = element_rect(fill = "transparent")  # Key background
            , legend.background = element_rect(fill= "transparent") # Legend background
            , legend.title = element_blank()
            #legend.position="none",
      )
    
    
    # p
    # ggsave(paste0("Boxplot_",output_file_name,"_orderBy",my_stat,"_",GoF,".tiff"),
    #        height=14, width=10,
    #        units="in", dpi=300, plot= p, bg = "transparent",
    #        path = paste(output_file))
    # ggsave(paste0("Boxplot_",output_file_name,"_orderBy",my_stat,"_",GoF,".png"),
    #        height=12, width=12,
    #        units="in", dpi=350, plot= p, bg = "transparent",
    #        path = paste(output_file))
    
    
  }  
  
  figure <- ggpubr::ggarrange(p1, p2, p3, 
                              #labels = c("a) Dakar", "b) Dar es Salaam", "c) Kampala", "d) Ouagadougou"),
                              #align = c("hv"),
                              ncol = 3, nrow = 1) 

  ggsave(paste0("Boxplot_",output_file_name,"_orderBy",my_stat,"_",GoF,".tiff"),
         height=15, width=40,
         units="in", dpi=300, plot= figure, bg = "transparent",
         path = paste(output_file)) 
  
  
}
  













