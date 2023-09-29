
library(tidyverse)
library(tidybayes)
library(bayesplot)
library(dplyr)


library(ggplot2)
library(corrplot)
library(ggcorrplot)
library(paletteer)

# Simulation Design

#----

extract_smdf_id <- function(big_sim_result, id){
  
  smdf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  for(i in ((id-1)*nsims+1):(id*nsims)){
    #Extract row corresponding to nth simulation
    row_big_sim <- big_sim_result[i,]
    
    ptemp <- dim(row_big_sim$r2d2$sm)[1]
    for(j in 1:length(names_list)){
      rtheta <- as.numeric( do.call(rbind, row_big_sim[[j]]$rtheta))[-ptemp]
      temp <- row_big_sim[[j]]$sm %>%
        add_column(rtheta= rtheta,  .before= "mean") %>% 
        add_column(summary_name= names_list[j], .before=1) %>% 
        add_column(simnr=i, .before=1 ) %>%
        add_column( id=id, .before=1 )
      
      smdf <- rbind(smdf, temp)
    }
    
    
  }
  
  smdf
  
}



extract_preddf_id <- function(big_sim_result, id){
  
  preddf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  for(i in ((id-1)*nsims+1):(id*nsims)){
    #Extract row corresponding to nth simulation
    row_big_sim <- big_sim_result[i,]
    
    for(j in 1:length(names_list)){
      
      temp <- as_tibble(row_big_sim[[j]]$perf) %>%
        add_column(name= names_list[j], .before=1) %>% 
        add_column(simnr=i, .before=1 ) %>%
        add_column( id=id, .before=1 )
      
      preddf <- rbind(preddf, temp)
    }
    
    
  }
  
  preddf
  
}


#----- Analyze a specific condition

# Condition to analyze

exp_name="full_sim_comp_result_new_GGG1_117"

names_list = c("r2d2", 
               "r2d2grouped",
               "horseshoe", 
               "rhorseshoe",
               "gigg"  )
               

# ggplot2 settings

#my_palette <- "ggthemes::Classic_Green_Orange_6"

my_theme <- theme(plot.title = element_text(size=18),
                  axis.text=element_text(size=14),
                  axis.line.y = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.text.x = element_text(size = 14),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_text(size = 16),
                  axis.text.x = element_text(size = 12),
                  axis.text.y = element_text(size = 12))+
  theme_set(theme_bw())


nids <-  nrow(sim_conds)

big_sim_result <- readRDS(paste0("R/big_simulation/final_results/",exp_name,".RDS"))

#big_sim_result <- do.call(rbind, big_sim_result)

for(id in 1:nids){
  
  smdf <- extract_smdf_id(big_sim_result, id)
  
  
  #sim_result <- readRDS(paste0("sim_comparison/final_results/",exp_name,"_" ,id,".RDS"))
  dir.create(paste0("R/big_simulation/final_results/", stringi::stri_sub(exp_name,-6,-1), "_", id ))
  #dir.create(paste0("comparison/final_results/", exp_name,"/",id))
  results_folder <- paste0("R/big_simulation/final_results/", stringi::stri_sub(exp_name,-6,-1), "_", id )
  
  #--- Post Processing
  
  #Extract all summaries
  
  {
    
    #smdf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
    #for(i in 1:length(names_list)){
    # smdf <- rbind(smdf,modified_summary(summary_name=names_list[i], 
    #                                   summaries= sim_result))
    #}
    
    smdf$summary_name <- as.factor(smdf$summary_name)
    
    smdf$summary_name <- recode( smdf$summary_name , 
                                 r2d2= "R2",
                                 r2d2grouped= "GR2",
                                 gigg= "GIGG")
    
    
    smdf <- smdf  %>% 
      mutate( bias= mean-rtheta ) %>% 
      mutate( sqerror= (mean-rtheta)^2 ) 
    
    #--- Summaries across sims
    
    #Description <- c("bs non-zero", "bs zero", "bs",
    #                "phis", "sigma", "R2")
    #stats <- list(mean=mean, sd=sd)
    
    
    
    #bs (nonzero)
    
    # bsnz <- smdf %>% mutate_at( vars(variable) , factor)  %>% 
    #   filter( grepl('^b' ,variable) &!grepl("Intercept", variable)   ) %>%
    #   filter(rtheta!=0) %>% 
    #   group_by(summary_name) %>% 
    #   dplyr::select(bias, sqerror, coverage, width)  %>% 
    #   summarise(across(c(bias, sqerror, coverage, width), .fns= stats))  
    # 
    # #bs (zero)
    # bsz <- smdf %>% mutate_at( vars(variable) , factor)  %>% 
    #   filter( grepl('^b' ,variable) &!grepl("Intercept", variable)   ) %>%
    #   filter(rtheta==0) %>% 
    #   group_by(summary_name) %>% 
    #   dplyr::select(bias, sqerror, coverage, width)  %>% 
    #   summarise(across(c(bias, sqerror, coverage, width), .fns= stats))  
    # 
    # #bs (all)
    # bsall <-smdf %>% mutate_at( vars(variable) , factor)  %>% 
    #   filter( grepl('^b' ,variable) &!grepl("Intercept", variable)   ) %>%
    #   group_by(summary_name) %>% 
    #   dplyr::select(bias, sqerror, coverage, width)  %>% 
    #   summarise(across(c(bias, sqerror, coverage, width), .fns= stats))  
    # 
    # 
    # Type I and Type II error
    # Description <- c("Type I bs zero",
    #                  "Type II bs non-zero")
    # 
    # stats <- list(sum=sum)
    # 
    #----- RMSEs
    
    # Calculate RMSE per simulation! 
    # Group by summary name and simnr
    # Calculate their distribution
    
    # rmsez <- smdf %>%
    #   group_by(summary_name, simnr) %>% 
    #   filter( grepl('^b' ,variable) &!grepl("Intercept", variable)   ) %>%
    #   filter(rtheta==0) %>% 
    #   dplyr::select(sqerror)  %>% 
    #   summarise(across(c(sqerror), .fns= list(sum=sum) ), n=n()) %>% 
    #   mutate( rmsez = sqrt(sqerror_sum/n) )
    # 
    
    # rmsenz <- smdf %>%
    #   group_by(summary_name, simnr) %>% 
    #   filter( grepl('^b' ,variable) &!grepl("Intercept", variable)   ) %>%
    #   filter(rtheta!=0) %>% 
    #   dplyr::select(sqerror)  %>% 
    #   summarise(across(c(sqerror), .fns= list(sum=sum) ), n=n()) %>% 
    #   mutate( rmsenz = sqrt(sqerror_sum/n) )
    # 
    # 
    # rmse <- smdf %>%
    #   group_by(summary_name, simnr) %>% 
    #   filter( grepl('^b' ,variable) &!grepl("Intercept", variable)   ) %>%
    #   dplyr::select(sqerror) %>% 
    #   summarise(across(c(sqerror), .fns= list(sum=sum) ), n=n()) %>% 
    #   mutate( rmse = sqrt(sqerror_sum/n) )
    
  }
  
  
  #----- PLOTS
  
  #---RMSE plots
  
  #rmse all plot 
  
  # rmse %>% 
  #   ggplot(aes(x= summary_name, y= rmse, fill= summary_name))+
  #   geom_boxplot(alpha=0.8) +
  #   scale_color_paletteer_d(my_palette)+
  #   ggtitle( paste("RMSE comparison"))+
  #   labs(
  #     y = "RMSE",
  #     x = "fit name")+ 
  #   my_theme 
  
  
  #results_sub %>%
  #  ggplot(aes(type, {{ var }})) +
  
  # rmse %>% 
  #   ggplot(aes(x= summary_name, y= rmse, fill= summary_name))+  
  #   geom_violin() +
  #   geom_boxplot(width = 0.1) +
  #   stat_summary(fun=median, geom = "point", color = "black", size = 2)
  #   
  
  # geom_hline( aes(yintercept = mean_var,  col=summary_name),
  #             data = rmse %>%
  #               group_by( summary_name ) %>%
  #               summarise(mean_var = mean(( rmse)), .groups = "drop"),
  #             linetype = 2
  # )
  # 
  
  # fname= paste0(results_folder,"/","rmse_id",id,".pdf",sep="")  
  # ggsave(filename=fname, plot = last_plot())
  
  
  
  
  # rmse zeros
  # rmsez %>% 
  #   ggplot(aes(x= summary_name, y= rmsez, fill= summary_name))+  
  #   geom_violin() +
  #   geom_boxplot(width = 0.1) +
  #   stat_summary(fun=median, geom = "point", color = "black", size = 2)+
  #   scale_color_paletteer_d(my_palette)+
  #   ggtitle( paste("RMSE noise comparison"))+
  #   labs(
  #     y = "RMSE",
  #     x = "Fit name")+ 
  #   my_theme
  # 
  # fname= paste0(results_folder,"/","rmsez_id",id,".pdf",sep="")  
  # ggsave(filename=fname, plot = last_plot())
  
  
  # rmse signal
  # rmsenz %>% 
  #   ggplot(aes(x= summary_name, y= rmsenz, fill= summary_name))+  
  #   geom_violin() +
  #   geom_boxplot(width = 0.1) +
  #   stat_summary(fun=median, geom = "point", color = "black", size = 2)+
  #   scale_color_paletteer_d(my_palette)+
  #   ggtitle( paste("RMSE zeros comparison"))+
  #   labs(
  #     y = "RMSE",
  #     x = "Fit name")+ 
  #   my_theme+
  #   ggtitle( paste("RMSE signal comparison"))+
  #   labs(
  #     y = "RMSE",
  #     x = "fit name")+ 
  #   my_theme
  
  
  # fname= paste0(results_folder,"/","rmsenz_id",id,".pdf",sep="")  
  # ggsave(filename=fname, plot = last_plot())
  
  # r2d2 
  # add line for mean and boxplot
  # r folders script  overlay boxplot and violin plot
  # remove naive
  
  # starting from a non optimal baseline: careful here and be
  # consistent
  # add a double dipping with logit r2d2 kl
  # we do not see any benefits
  # we see its not beneficial and have 2 starting points
  # Shrinkage plots
  
  # r2d2
  
  
  #--- Plots
  # Shrinkage Plots
  
  new_names_list <- levels(as.factor(smdf$summary_name))
  
  for(i in 1:length(new_names_list) ){
    
    
    smdf %>%  mutate_at( vars(variable) , factor)  %>% 
      filter( summary_name == new_names_list[i]  )  %>%
      filter( grepl('^b' ,variable) &!grepl("b_Intercept", variable)   ) %>%
      ggplot(aes(x = rtheta, y = mean, color= as.factor(simnr) )) +
      geom_point(size= 0.3, alpha=0.9) +
      labs(
        x = "Real value of b",
        y = "Posterior mean") +
      ggtitle(new_names_list[i] )+
      my_theme+
      theme(legend.position="none")
    
    fname= paste0(results_folder,"/","post_shrinkage_",new_names_list[i],"_id",
                  id,".pdf",sep="")  
    ggsave(filename=fname, plot = last_plot())
    
  }
  
  
  
  #---
  
  
  
  # temp <- smdf %>%  mutate_at( vars(variable) , factor)  %>% 
  #   filter( summary_name == new_names_list[i]  )  %>%
  #   filter( grepl('^b' ,as.factor(variable)) &!grepl("Intercept", variable)   ) %>% 
  #   filter(rtheta!=0)
  
  
  # smdf %>%  mutate_at( vars(variable) , factor)  %>% 
  #   filter( summary_name == new_names_list[i]  )  %>%
  #   filter( grepl('^b' ,as.factor(variable)) &!grepl("Intercept", variable)   ) %>%
  #   ggplot(aes(x = rtheta, y = mean, color= as.factor(simnr) )) +
  #   geom_point(size= 1, alpha=1) +
  #   labs(
  #     x = "Real value of b",
  #     y = "Posterior mean") +
  #   ggtitle(new_names_list[i] )+
  #   my_theme
  # 
  # fname= paste0(results_folder,"/","post_shrinkage_",new_names_list[i],"_id",
  #               id,".pdf",sep="")  
  # ggsave(filename=fname, plot = last_plot())
  
  #----------------
  
  
  #---- Prediction summary
  
  
  
  
  #preddf <- as.data.frame(matrix(nrow =0  , ncol=0)) 
  # for(i in 1:length(names_list)){
  #   preddf <- rbind(preddf,modified_performance(names_list[i], sim_result)) 
  # }
  
  
  preddf <- extract_preddf_id(big_sim_result, id)
  
  preddf <- preddf %>% group_by(simnr) %>% 
    mutate(delta_rmse_b = rmse_b - min(rmse_b, na.rm = TRUE),
           delta_rmse_b0 = rmse_b0 - min(rmse_b0, na.rm = TRUE),
           delta_rmse_bpp = rmse_bpp - min(rmse_bpp, na.rm = TRUE),
           delta_lpd_train= lpd_train - max(lpd_train, na.rm = TRUE),
           delta_lpd_test= lpd_test - max(lpd_test, na.rm = TRUE),
           delta_elpd_loo= elpd_loo - max(elpd_loo, na.rm = TRUE),
           delta_elpd_loo_sd= elpd_loo_sd - min(elpd_loo_sd, na.rm = TRUE),
           delta_p_loo= p_loo - min(p_loo, na.rm = TRUE)) %>% 
    ungroup()
  
  
  preddf$name <- as.factor(preddf$name)
  preddf$name <- recode( preddf$name , 
                         r2d2= "R2",
                         r2d2grouped= "GR2",
                         gigg= "GIGG")
  
  
  
  #best within each data set instead of model
  
  # RMSE b
  
  preddf %>% dplyr::select(name, starts_with("rmse" ) & ends_with("_b")) %>% 
    pivot_longer(!name, names_to = "rmse", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ rmse, scales = "free", nrow=2) +
    ggtitle( paste("RMSE coefficients"))+
    labs(
      y = "RMSE",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","rmse_post_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  preddf %>% dplyr::select(name, starts_with("rmse" ) & ends_with("_b0")) %>% 
    pivot_longer(!name, names_to = "rmse", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ rmse, scales = "free", nrow=2) +
    ggtitle( paste("RMSE coefficients"))+
    labs(
      y = "RMSE",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","rmse_post0_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  
  preddf %>% dplyr::select(name, starts_with("rmse" ) & ends_with("_bpp")) %>% 
    pivot_longer(!name, names_to = "rmse", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ rmse, scales = "free", nrow=2) +
    ggtitle( paste("RMSE coefficients"))+
    labs(
      y = "RMSE",
      x = "Type of rmse")+ 
    my_theme
  
  
  fname= paste0(results_folder,"/","rmse_postpp_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  #Train rmse 
  
  preddf %>% dplyr::select(name, starts_with("rmse" ) & ends_with("_train")) %>% 
    pivot_longer(!name, names_to = "rmse", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ rmse, scales = "free", nrow=2) +
    ggtitle( paste("Comparison Train data"))+
    labs(
      y = "RMSE",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","TrainRMSE_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  #Test rmse 
  
  preddf %>% dplyr::select(name, starts_with("rmse" ) & ends_with("_test")) %>% 
    pivot_longer(!name, names_to = "rmse", values_to = "value") %>% 
    #filter( value < 5) %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+facet_wrap(~ rmse, scales = "free", nrow=2) +
    ggtitle( paste(" Comparison Test"))+
    labs(
      y = "RMSE",
      x = "Type of rmse")+ 
    my_theme
  
  fname= paste0(results_folder,"/","TestRMSE_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  # elpd
  
  preddf %>% dplyr::select(name, starts_with("lpd" ) & ends_with("_train")) %>% 
    pivot_longer(!name, names_to = "lpd", values_to = "value") %>% 
    filter (value > -10000) %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    ggtitle( paste(" elpd Train comparison"))+
    labs(
      y = "elpd",
      x = "name")+ 
    my_theme
  
  fname= paste0(results_folder,"/","lpdtrain_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  preddf %>% dplyr::select(name, "lpd_test" ) %>% 
    pivot_longer(!name, names_to = "lpd", values_to = "value") %>% 
    filter( value > -10000) %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    ggtitle( paste("elpd test comparison"))+
    labs(
      y = "elpd",
      x = "name")+
    my_theme
  
  fname= paste0(results_folder,"/","lpdtest_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  
  preddf %>% dplyr::select(name, p_loo) %>% 
    pivot_longer(!name, names_to = "p_loo", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    ggtitle( paste("p_loo"))+
    labs(
      y = "p_loo",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","p_loo_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  #--- delta rmse
  
  # delta rmse b
  preddf %>% dplyr::select(name, starts_with("delta_rmse" ) & ends_with("_b")) %>% 
    pivot_longer(!name, names_to = "delta_rmse", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ delta_rmse, scales = "free", nrow=2) +
    ggtitle( paste("delta RMSE coefficients"))+
    labs(
      y = "delta RMSE",
      x = "Model")+ 
    my_theme
  
  
  
  fname= paste0(results_folder,"/","delta_rmse_postb_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  # delta rmse b0
  preddf %>% dplyr::select(name, starts_with("delta_rmse" ) & ends_with("_b0")) %>% 
    pivot_longer(!name, names_to = "delta_rmse", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ delta_rmse, scales = "free", nrow=2) +
    ggtitle( paste("delta RMSE coefficients"))+
    labs(
      y = "delta RMSE",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","delta_rmse_postb0_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  # delta rmse bpp
  preddf %>% dplyr::select(name, starts_with("delta_rmse" ) & ends_with("_bpp")) %>% 
    pivot_longer(!name, names_to = "delta_rmse", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ delta_rmse, scales = "free", nrow=2) +
    ggtitle( paste("delta RMSE coefficients"))+
    labs(
      y = "delta RMSE",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","delta_rmse_postbpp_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  
  # delta elpd
  
  #delta lpd train
  preddf %>% dplyr::select(name, delta_lpd_train) %>% 
    pivot_longer(!name, names_to = "delta_lpd", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ delta_lpd, scales = "free", nrow=2) +
    ggtitle( paste("delta lpd train"))+
    labs(
      y = "delta lpd",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","delta_lpd_train_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  # delta lpd test
  preddf %>% dplyr::select(name, delta_lpd_test) %>% 
    pivot_longer(!name, names_to = "delta_lpd", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ delta_lpd, scales = "free", nrow=2) +
    ggtitle( paste("delta lpd test"))+
    labs(
      y = "delta lpd",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","delta_lpd_test_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
  #delta elpd_loo
  preddf %>% dplyr::select(name, delta_elpd_loo ) %>% 
    pivot_longer(!name, names_to = "delta_elpd", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ delta_elpd, scales = "free", nrow=2) +
    ggtitle( paste("delta elpd loo "))+
    labs(
      y = "delta elpd",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","delta_elpd_loo_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  #delta elpd_loo_sd
  preddf %>% dplyr::select(name, delta_elpd_loo_sd) %>% 
    pivot_longer(!name, names_to = "delta_elpd", values_to = "value") %>% 
    ggplot(aes(x= as.factor(name), y= value, fill=as.factor(name)))+
    geom_violin() +
    geom_boxplot(width = 0.1) +
    stat_summary(fun=median, geom = "point", color = "black", size = 2)+
    scale_color_paletteer_d(my_palette)+
    facet_wrap(~ delta_elpd, scales = "free", nrow=2) +
    ggtitle( paste("delta elpd loo "))+
    labs(
      y = "delta elpd",
      x = "Model")+ 
    my_theme
  
  fname= paste0(results_folder,"/","delta_elpd_loo_sd_id",
                id,".pdf",sep="")  
  ggsave(filename=fname, plot = last_plot())
  
  
}


# TODO:

# 
# 
# #--- Type I and Type II error
# 
# quantiles <- c("q1", "q2.5","q5", "q10", "q15", "q20" , "q25", "q33" , "q40" , "q50",       
#                "q60", "q67","q75", "q80","q85",  "q90", "q95", "q97.5","q99")
# 
# errordf <- matrix(0,0,0)
# for(i in 1:floor(length(quantiles)/2)){
#   
#   lower=quantiles[i]
#   upper=quantiles[length(quantiles)-i+1]
#   
#   smL <-  list(upper=upper,
#                lower=lower,
#                smdf= smdf %>% dplyr::select(!c(upper, lower)))
#   
#   errordf <- rbind(errordf,sim_error(smL))
#   
# }
# 
# alphas=c(1,2.5, 5, 10,15, 20, 25, 33, 40 )*2/100
# alphas=(1-alphas)
# 
# 
# nl <- nlevels(factor(smdf$summary_name))
# 
# # Type I error
# errordf %>% filter(Description == "Type I bs zero") %>% 
#   add_column(alphas=rep(alphas, each= nl )) %>% 
#   ggplot(aes(  x= alphas, y= error, color= factor(summary_name)))+
#   geom_line( linewidth=0.7)+
#   geom_point(aes(colour = factor(summary_name)), size = 1.5, alpha=1)+
#   scale_colour_manual(values= viridis::viridis(nl), name="p") +  
#   #facet_wrap(~ summary_name, scales = "free", ncol= nl) +
#   labs(
#     x = expression(paste("Credibility ", 1-alpha)),
#     y = "Type I") +
#   ggtitle("Type I error")+
#   my_theme
# 
# fname= paste0(results_folder,"/","Type_I_id",
#               id,".pdf",sep="")  
# ggsave(filename=fname, plot = last_plot())
# 
# # Type II error
# errordf %>% filter(Description == "Type II bs non-zero") %>% 
#   add_column(alphas=rep(alphas, each= nl )) %>% 
#   ggplot(aes(  x= sort(alphas), y= error, color= factor(summary_name)))+
#   geom_line( linewidth=0.7)+
#   geom_point(aes(colour = factor(summary_name)), size = 1.5, alpha=1)+
#   scale_colour_manual(values= viridis::viridis(nl), name="p") +  
#   labs(
#     x = expression(paste("Credibility ", 1-alpha)),
#     y = "Type II") +
#   ggtitle("Type II error")+
#   my_theme
# 
# 
# fname= paste0(results_folder,"/","Type_II_id",
#               id,".pdf",sep="")  
# ggsave(filename=fname, plot = last_plot())
# 
# 
# errordf %>% filter(Description == "Type II bs non-zero") %>% 
#   add_column(alphas=rep(alphas, each=nl)) %>% 
#   ggplot(aes(  x= sort(alphas), y= cov_error, color= factor(summary_name)))+
#   geom_line( linewidth=0.7)+
#   geom_point(aes(colour = factor(summary_name)), size = 1.5, alpha=1)+
#   scale_colour_manual(values= viridis::viridis(nl), name="p") +  
#   #facet_wrap(~ summary_name, scales = "free", ncol=nl) +
#   labs(
#     x = expression(paste("Credibility ", 1-alpha)),
#     y = "Proportion of coverage") +
#   ggtitle("Coverage ")+
#   my_theme 
# 
# fname= paste0(results_folder,"/","coverage_error_id",
#               id,".pdf",sep="")  
# ggsave(filename=fname, plot = last_plot())
# 
# 
# 
# # TODO ROC 
# 
# 
# # Some summaries
# bsnz #non zero
# bsz #zero
# bsall #all
# rmsez #zeros
# rmsenz #non zeroes
# rmse #all
# 
# 
# # TODO:
# # ROC curves
# 
#   
#  
#  
#  our r2 uses more information than the r2 in lm

# Conditions in which the grouping is beneficial. Their simulation conditions
# do not allow to find this situation 
# change importance of the grouping such that the grouped priors are better
# are groups really useful?

# If we had more efficient ways it would be doneso via the prior



