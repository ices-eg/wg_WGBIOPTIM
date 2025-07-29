#######################################################################################################
#######################################################################################################
#######################################################################################################
##
##
##   Biological sampling optimization (Script "SampleOptim")
##   Developed by: Patricia Goncalves (patricia@ipma.pt)
##   Last version development period: may 2023
##   Version: v4.1
##
##   Reference:
##   Gonçalves, Patrícia 2019. "SampleOptim" a data analysis R-tool to optimize fish sampling for
##   biological parameters as input on fish stock assessment.
##
##
##
#######################################################################################################
#######################################################################################################
######################################################Ó#################################################

# In this script we plot results
#' we use functions to 
#' - 
#' 
#'  


#============================================


##Packages:

library(FSA)
library(FSAdata)
library(nlstools)
library(reshape)
library(ggplot2)
library(ggthemes)
library(cvTools)
library(dplyr)
library("robustbase")
library(MASS)
library(psyphy)
library(boot)
library(RCurl)
library(dplyr)
library(tidyr)
library(raster)

library(here)

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
############# 3. Simulation results (Plots preparation)
#############   Aggregate all the data available from the different simulation runs
#############   For all the matrices generated from the code

output_dir <- here('output')
# 
output_sub_dir <- file.path(output_dir, "simulation_results")
dir.create(output_sub_dir, showWarnings = FALSE)


SEP <-","

source("00_helper_functions.R") # or maybe on the first script? Is this standalone?

#############################################
############   READ INPUT SETTINGS FROM FILE ###################################

# param defined in the 2_Simulations.R script


########
getAllData<-function(tabname){
  ress <- NULL
  for(n in otolitSet){
    tval <- read.table(paste0(output_dir, tabname, n, ".csv"), sep = SEP, header = TRUE)
    ress <- rbind(ress, tval)
  }
  return(ress)
}

# check that the files are there
list.files(here::here(output_dir))

#### Data of age, length by year (from the individuals selected in each simulation run)

data_selected_lt_age <- getAllData("dados_lt_age_")
summary(data_selected_lt_age)

############################################################################
####   von Bertalanffy growth model parameters by year and simulation run
############################################################################
#############################################################################

simulvbgm_param_year <- getAllData("results_simulvbgm_")
summary(simulvbgm_param_year)

####For all the years, to compare the VGBGM parameters between years
####Figure 7 - Summary of parameters by year for the full set of simulations (n=100), for j's (number of selected otoliths)

# the following figures are all the same except for y = K, t0, Linf
# -----------------------
#' Create a Boxplot
#'
#' This function generates a boxplot for a specified variable from the input data.
#' The plot is saved as a PNG file in the specified output directory.
#'
#' @param data A data frame containing the data to be plotted.
#' @param y_variable A string representing the name of the variable to be plotted on the y-axis.
#' @param output_dir A string specifying the directory where the output PNG file will be saved.
#' @param timeInterval A string representing the time interval for the plot.
#' @param y_limits Optional; a numeric vector of length 2 specifying the limits for the y-axis.
#'
#' @return NULL; the function saves a boxplot as a PNG file.

create_boxplot <- function(data, y_variable, output_dir, timeInterval, y_limits = NULL) {
  
  # Construct the file name for the output PNG file
  file_name <- file.path(paste0(output_dir, paste0("Fig7_", y_variable, "_VBGM_", timeInterval, ".png")))
  png(file = file_name)
  
  plot <- ggplot(data, aes(x = factor(type), y = .data[[y_variable]], fill = factor(type))) +
    geom_boxplot(outlier.shape = NA) +  # Add boxplot layer without outliers
    xlab("Number of otoliths selected by length class (cm)") +
    theme_classic() +
    facet_wrap(~year) + # Create separate panels for each year
    sample_theme() # Custom theme defined in 0_helper_functions.R
  
  # Check if y_limits is provided and apply it to the plot
  if (!is.null(y_limits)) {
    plot <- plot + ylim(y_limits)  # Set y-axis limits if specified
  }
  
  # Print the plot to the device
  print(plot)
  
  # Close the PNG device
  dev.off()
}

# Create the boxplots
create_boxplot(simulvbgm_param_year, "K", output_dir, timeInterval)
create_boxplot(simulvbgm_param_year, "t0", output_dir, timeInterval)
create_boxplot(simulvbgm_param_year, "Linf", output_dir, timeInterval, c(20, 100))


####################################################################################
#####  Mean length at age data original and by simulation run
###################################################################################
####################################################################################

table_mla_sims <- getAllData("table_original_simulsub_mla_")
summary(table_mla_sims)

###############################################################################################################
###### Figure 8 - compare mean length at age (oringinal versus simulation) by year and simulation run
#######

years <- unique(table_mla_sims$year)

for(i in years){
  fig8_mla <- ggplot(data = subset(table_mla_sims, year==i), aes(x=factor(age), y=m_lt,fill=factor(type))) +
    geom_boxplot(outlier.shape = NA)+
    xlab("age") +
    ylab("mean length at age (cm)") + 
    scale_color_brewer(palette = "Paired") + 
    facet_wrap(~year) +
    ggtitle(i) +
    sample_theme()
  ggsave(fig8_mla, file = paste0(output_sub_dir,"Fig8_mla_",i,timeInterval,".png"),width=14, height = 10, units="cm")
}



###########################################################################################
### Standard deviation from length at age by simualtion run
#####################################################################
#############################################################################################


table_sdla_sims <- getAllData("table_original_simulsub_sd_")
summary(table_sdla_sims)

################################################################################################
###### Figure 9 - compare the sd length at age (oringinal versus simulation) by year and simulation run
#######
#table_sdla_sims<-table_sdla_sims[!is.na(table_sdla_sims$sd_lt),]##remove NAs on data[table_sdla_sims$year==years[[nb]],]

years <- unique(table_sdla_sims$year)

for(i in years){
  fig9_sdla<- ggplot(data=subset(table_sdla_sims, year==i),aes(x=factor(age), y=sd_lt,fill=factor(type))) +
    geom_boxplot(outlier.shape =NA) + 
    xlab("age") +
    ylab("sd length at age (cm)") +
    scale_color_brewer(palette = "Paired") + 
    facet_wrap(~year) +
    ggtitle(i)+
    sample_theme()
  ggsave(fig9_sdla, file=paste0(output_sub_dir,"Fig9_sdla_",i,timeInterval,".png"),width=14, height = 10, units="cm")
}


##################################################################################################
###
### Data combine (table_mla_sims, table_sdla_sims)
##################################################################################################

data_mla_sd <- merge(table_mla_sims, table_sdla_sims)
data_mla_sd <- data_mla_sd[complete.cases(data_mla_sd$sd_lt),]
data_mla_sd$coefv <- data_mla_sd$m_lt/data_mla_sd$sd_lt


years <- unique(data_mla_sd$year)

for(i in years){
  fig9_cvla <- ggplot(data=subset(data_mla_sd, year==i),aes(x=factor(age), y=coefv,fill=factor(type))) +
    geom_boxplot(outlier.shape =NA) + 
    xlab("age") + 
    ylab("CV (lenght)") + 
    scale_color_brewer(palette = "Paired") +
    facet_wrap(~year) + 
    ggtitle(i) +
    sample_theme()
  ggsave(fig9_cvla, file=paste0(output_sub_dir,"Fig9a_cvla_",i,timeInterval,".png"),width=14, height = 10, units="cm")
}


####################################################################################################
### Stats (mape, rmspe, mspe) from each simualations run and year
######################################################################
######################################################################################################

stats_simul <- getAllData("table_res_stat_")
summary(stats_simul)

###############################################################################################################
###### Figure 10 - compare the stats (mape, mspe, rtmspe) by year and by simulation type (number of otoliths/length class)
#######

# the following Fig10 figures are all the same except for y = mspe, mape or rtmspe

#---------------------------------------------
#' Plot Metrics Function
#'
#' This function generates a plot for specified metrics from simulation statistics.
#' 
#' @param stats_simul A data frame containing simulation statistics.
#' @param output_dir A string representing the directory where the output file will be saved.
#' @param timeInterval A string indicating the time interval for the plot.
#' @param metric A string specifying the metric to be plotted. Options are 'mspe', 'mape', or 'rtmspe'.
#' 
#' @return A plot saved as a PNG file in the specified output directory.

# Function to plot figures 10 above for MSPE, MAPE, and RTMSPE
plot_metrics <- function(stats_simul, output_dir, timeInterval, metric) {
  
  # Construct the file name based on the metric
  file_name <- file.path(output_dir, paste0(paste0("Fig10_", metric, "_", timeInterval, ".png")))
  # Open a PNG device
  png(file = file_name)
  
  # Create the plot
  plot <- ggplot(stats_simul, aes_string(x = "factor(type)", y = metric, group = 1)) +
    geom_step() +
    xlab("Number of Otoliths Selected by Length Class (cm)") +
    facet_wrap(~year) +
    sample_theme() # in 00_helper_functions
  
  # Print the plot to the device
  print(plot)
  
  # Close the PNG device
  dev.off()
}#



#'----------------------
# Use the function
plot_metrics(stats_simul, output_dir, timeInterval, "mspe")
plot_metrics(stats_simul, output_dir, timeInterval, "mape")
plot_metrics(stats_simul, output_dir, timeInterval, "rtmspe")



#--------------
###### Normalized data stats

##Mean by stats by year (any later use??)
mean_mspe <- stats_simul %>% 
  group_by(year) %>% 
  summarize(Mean = mean(mspe, na.rm=TRUE))
mean_mape <- stats_simul %>% 
  group_by(year) %>% 
  summarize(Mean = mean(mape, na.rm=TRUE))
mean_rtmspe <- stats_simul %>% 
  group_by(year) %>% 
  summarize(Mean = mean(rtmspe, na.rm=TRUE))



### Normalization of data 
nor_stats <- stats_simul %>%
  group_by(year) %>%
  mutate(
    ### Normalization by stats and year (max)
    norm_mspe = mspe / max(mspe),
    norm_mape = mape / max(mape),
    norm_rtmspe = rtmspe / max(rtmspe), 
    ### Normalization by stats and year (mean)
    norm_mn_mspe = mspe/mean(mspe),
    norm_mn_mape = mape/mean(mape),
    norm_mn_rtmspe = rtmspe/mean(rtmspe))


# Call the plotting function for each normalized metric
plot_metrics(nor_stats, output_dir, timeInterval, "norm_mspe")
plot_metrics(nor_stats, output_dir, timeInterval, "norm_mape")
plot_metrics(nor_stats, output_dir, timeInterval, "norm_rtmspe")

plot_metrics(nor_stats, output_dir, timeInterval, "norm_mn_mspe")
plot_metrics(nor_stats, output_dir, timeInterval, "norm_mn_mape")
plot_metrics(nor_stats, output_dir, timeInterval, "norm_mn_rtmspe")


###################################################################################################################
####### Comparison of von Bertalanffy growth model parameters estimated by year
###### and for each j (number of otoliths by length class) the result of the 100 simulations aggregated
###################################################################################################################

# Define a function to generate plots for Fig11
plot_fig11 <- function(stats_simul, output_dir, timeInterval, metric) {
  
  # Construct the file name based on the metric
  file_name <- file.path(output_dir, paste0("Fig11_", metric, "_VBGM_", timeInterval, ".png"))
  
  # Open a PNG device
  png(file = file_name)
  
  # Create the plot
  plot <- ggplot(stats_simul, aes_string(x = "factor(type)", y = metric)) +
    geom_point(size = 2, colour = "green") +
    xlab("Number of Otoliths Selected by Length Class (cm)") +
    facet_wrap(~year) +
    sample_theme()  # Assuming sample_theme() is defined elsewhere
  
  # Print the plot to the device
  print(plot)
  
  # Close the PNG device
  dev.off()
}


plot_fig11(stats_simul, output_dir, timeInterval, "k")
plot_fig11(stats_simul, output_dir, timeInterval, "t0")
plot_fig11(stats_simul, output_dir, timeInterval, "Linf")


################################################################################################################
###    Based on the VBGM by sim and year, the length for each selected fish was predicted based on age data
###########################################################################################################
#################################################################################################################

lpredict_vbsim<-getAllData("vb_predict_melt_")
summary(lpredict_vbsim)

###############################################################################################################
###### Figure 12 - compare the stats (mape, mspe, rtmspe) by year and by simulation type (number of otoliths/length class)
#######

# To revise. Not sure about the best way to show this. 

Fig12_predictLt <- file.path(paste(output_dir,"Fig12_predictLt_", timeInterval, ".png", sep = ""))
png(file=Fig12_predictLt)
fig12_predictLt <- ggplot(data=lpredict_vbsim, aes(x=ID_ind, y=pred_lt,colour=type)) +
  geom_line() +
  xlab("ID_ind")+
  ylab("Predicted length") +
  facet_wrap(~year) +
  sample_theme()
print(fig12_predictLt)
dev.off()


Fig12a_predictLt <- file.path(paste(output_dir,"Fig12a_predictLt_", timeInterval, ".png", sep = ""))
png(file=Fig12a_predictLt)
fig12a_predictLt <- ggplot(data=lpredict_vbsim, aes(x=ID_ind, 
                                                  y=pred_lt,
                                                  colour=factor(year))) +
  geom_point() + 
  xlab("ID_ind") +
  ylab("Predicted length")+
  facet_wrap(~factor(type))+
  sample_theme()
print(fig12a_predictLt)
dev.off()



  ###############################################################################################################################
  ###############################################################################################################################
  ### Data from biological samples
  ###############################################################################################################################
  ###############################################################################################################################

  dados_bio_simul <- getAllData("dados_bio_")
  summary(dados_bio_simul)
  # Clean the dataset by removing NAs
  dados_bio_simul <- dados_bio_simul[complete.cases(dados_bio_simul$wt),]
  
  # plot function
  figwg_plot <- function(data, year, x_var, x_label, y_var, y_label, facet_by_quarter = FALSE) {
    
    # Ensure the specified x_var is treated as a factor
    data[[x_var]] <- as.factor(data[[x_var]])
    
    # Create the plot using aes_string for dynamic x-axis
    plot <- ggplot(data = subset(data, year == year),
                   aes_string(x = x_var, y = y_var, fill = "factor(type)")) +
      geom_boxplot(outlier.shape = NA) +
      xlab(x_label) +  # Use the provided x-axis label
      ylab(paste0(y_label, " at ",  x_label)) +
      scale_color_brewer(palette = "Paired") +
      ggtitle(paste("Year:", year)) +
      sample_theme()  # sample_theme() is defined in 00_helper_functions.R
    
    # Conditional faceting based on user input
    if (facet_by_quarter) {
      plot <- plot + facet_grid(quarter ~ year)
    } else {
      plot <- plot + facet_wrap(~year)
    }
    
    return(plot)
  } # figwgplot
  
  
  
  ###### MEAN WEIGHT at age

  table_meanweight <- group_by(dados_bio_simul, age, year, ID_sim, type) %>% 
    summarize(m_wg = mean(wt))
  
  years <- unique(table_meanweight$year)
  
  ###### MEAN WEIGHT at length
  table_meanweight_Lt <- group_by(dados_bio_simul, Lt, year, ID_sim,type) %>% 
    summarize(m_wg = mean(wt))

  
  ###### SD WEIGHT at age
  table_sdweight <- dados_bio_simul  %>% 
    group_by(age, year, ID_sim,type) %>%
    summarise(sdwt = sd(wt))
  summary(table_sdweight)
  

  ###### CV WEIGHT at age
  table_cvweight<-group_by(dados_bio_simul, age, year, ID_sim,type) %>% 
    summarize(CVwg = cv(wt))
  
  #---------------------------
  # Loop for yearly plots
  for (i in years) {
    Figweight <- figwg_plot(table_meanweight, i, x_var = "age", x_label = "Age", y_var =  "m_wg", y_label = "Mean Weight", facet_by_quarter = FALSE)
    ggsave(Figweight, file = paste0(output_sub_dir, "Figweight", i, timeInterval, ".png"), width = 14, height = 10, units = "cm")
  
    Figweight_lt<- figwg_plot(table_meanweight_Lt, i, x_var = "Lt", x_label = "Length", y_var =  "m_wg", y_label = "Mean Weight", facet_by_quarter = FALSE)
    ggsave(Figweight_lt, file=paste0(output_sub_dir,"Figweight_lt",i,timeInterval,".png"),width=14, height = 10, units="cm")
  
    Figweightsd<- figwg_plot(table_sdweight, i, x_var = "age", x_label = "Age", y_var =  "sdwt", y_label = "SD Weight", facet_by_quarter = FALSE) 
    ggsave(Figweightsd, file=paste0(output_sub_dir,"Figweightsd",i,timeInterval,".png"),width=14, height = 10, units="cm")
    
   
    FigweightCV<- figwg_plot(table_cvweight, i, x_var = "age", x_label = "Age", y_var =  "CVwg", y_label = "CV Weight", facet_by_quarter = FALSE) 
    ggsave(FigweightCV, file=paste0(output_sub_dir,"FigweightCV",i,timeInterval,".png"),width=14, height = 10, units="cm")
    
    
     
    } #end loop over years
  
  # Check if timeInterval is "Q" for quarterly data
  if (timeInterval == "Q") {
    table_meanweight <- group_by(dados_bio_simul, age, year, quarter, ID_sim, type) %>%
      summarize(m_wg = mean(wt))
   
    table_meanweight_Lt_quarter <- group_by(dados_bio_simul, Lt, year, quarter, ID_sim,type) %>% 
      summarize(m_wg = mean(wt))
    write.csv(table_meanweight_Lt_quarter,"output/simulation_results/table_meanweight_Lt_quarter.csv")
 
    table_sdweight_quarter<-dados_bio_simul  %>% group_by(age, year, quarter, ID_sim,type) %>% summarise(sdwt = sd(wt))
    write.csv(table_sdweight_quarter,"output/simulation_results/table_sdweight_quarter.csv")
 
    table_cvweight_quarter<-group_by(dados_bio_simul, age, year, quarter, ID_sim,type) %>% summarize(CVwg = cv(wt))
    write.csv(table_cvweight_quarter,"output/simulation_results/table_cvweight_quarter.csv")
    
     
    # Loop for quarterly plots
    for (i in years) {
      Figweight_quarter <- figwg_plot(table_meanweight, i, x_var = "age", x_label = "Age", y_var =  "m_wg", y_label = "Mean Weight", facet_by_quarter = TRUE)
      ggsave(Figweight_quarter, file = paste0(output_sub_dir, "Figweight_quarter", i, timeInterval, ".png"), width = 14, height = 10, units = "cm")
    
      Figweight_ltquarter <- figwg_plot(table_meanweight_Lt_quarter, i, x_var = "Lt", x_label = "Length", y_var =  "m_wg", y_label = "Mean Weight", facet_by_quarter = TRUE) 
      ggsave(Figweight_ltquarter, file=paste0(output_sub_dir,"Figweight_ltquarter",i,timeInterval,".png"),width=14, height = 10, units="cm")
      
      Figweightsd_quarter<- figwg_plot(table_sdweight_quarter, i, x_var = "age", x_label = "Age", y_var =  "sdwt", y_label = "SD Weight", facet_by_quarter = TRUE)  
      ggsave(Figweightsd_quarter, file=paste0(output_sub_dir,"Figweightsd_quarter",i,timeInterval,".png"),width=14, height = 10, units="cm")
                                   
      
      FigweightCV_quarter<- figwg_plot(table_cvweight_quarter, i, x_var = "age", x_label = "Age", y_var =  "CVwg", y_label = "CV Weight", facet_by_quarter = TRUE)  
      ggsave(FigweightCV_quarter, file=paste0(output_sub_dir,"FigweightCV_quarter",i,timeInterval,".png"),width=14, height = 10, units="cm")
      
      
      } #loop over years
  } # if timeInterval Q
#'------------------  



  ####################################################################################################################
  ####################################################################################################################
  ##### Data from the maturity ogive model adjustment
  ###################################################################################################################
  ####################################################################################################################
  if(param.age_only==FALSE){

  dados_mo_simul<-getAllData("table_res_mo_")
  summary(dados_mo_simul)


  dados_mo_re <- dados_mo_simul %>% gather(type_mat, L_mean, L25:L75)

  years<- unique(dados_mo_re$year) ##list of years on the samples data

  for(i in years)
  {
    Figmo <- file.path(paste(output_sub_dir,"FigmoLt_", i, ".png", sep = ""))
    png(file=Figmo)
    plotmo<-
      ggplot(dados_mo_re[dados_mo_re$year==i,], aes(x=factor(type),y=L_mean, fill=type_mat))+xlab("Number of otoliths selected by length class (cm)")+ ylab("Length")+
      geom_boxplot()+ theme_classic() +
      theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
            axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
    print(plotmo)
    dev.off()
  }

  ###Sem outliers



  for(i in years)
  {
    Figmo_sna <- file.path(paste(output_sub_dir,"Figmo_snaLt_", i, ".png", sep = ""))
    png(file=Figmo_sna)
    plotmo_sna<-
      ggplot(dados_mo_re[dados_mo_re$year==i,], aes(x=factor(type),y=L_mean, fill=type_mat))+xlab("Number of otoliths selected by length class (cm)")+ ylab("Length")+
      geom_boxplot(outlier.shape = NA)+ theme_classic() +
      theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
            axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
    print(plotmo_sna)
    dev.off()
  }

}

##########################################################################################################
###########################################################################################################
#############################     END CODE ;)    ###########################################################
###########################################################################################################
###########################################################################################################

