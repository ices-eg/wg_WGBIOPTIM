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
#######################################################################################################

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


###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
############# 3. Simulation results (Plots preparation)
#############   Agregate all the data available from the different simulation runs
#############   For all the matrixes generated from the code

#setwd("C:\\Users\\patricia\\Documents\\WHB\\PNAB2018\\optimizacao_amostragem2018\\excercicio_teste_whb")

setwd("C:\\Users\\patricia\\Documents\\WHB\\PNAB2021\\Optimizacao_amostragem_PNAB\\04.Testes_simulaPG\\WHB") ##Set directory
output_dir<-"output/"

output_sub_dir<-paste(output_dir,"simulation_results/", sep="/")
dir.create(output_sub_dir, showWarnings = FALSE)


SEP<-","


getSep<-function(fname){
  linha<-readLines(con = fname, n = 1L, ok = TRUE, warn = TRUE,encoding = "unknown", skipNul = TRUE)
  a<-unlist(strsplit(linha, ","))
  b<-unlist(strsplit(linha, ";"))
  if(length(b)>length(a)) return(";")
  else return(",")
}

#############################################
############   READ INPUT SETTINGS FROM FILE ###################################
#File Header
#Variable name;Mandatory;Variable.value;Definition
#
NAME=0
MANDATORY=1
VALUE=2
DEFINITION=3
###### Input settings file
fname<-"input_params.csv"
sep_is<-getSep(fname)
sett<- read.table(fname,sep=sep_is, header=F,row.names = c(), stringsAsFactors=FALSE)
stt<- data.frame(t(sett), stringsAsFactors=FALSE)
rownames(stt)<- c()
colnames(stt) <- as.character(unlist(stt[1,]))
stt <- stt[-1, ]

param.age_only=stt$AGE_ONLY[[VALUE]]
param.timeStrata=stt$TIME_STRATA[[VALUE]]
param.min_otol=as.numeric(stt$MIN_OTOL[[VALUE]])
param.max_otol=as.numeric(stt$MAX_OTOL[[VALUE]])
param.interval_otol=as.numeric(stt$interval_OTOL[[VALUE]])
param.extra_otol=as.array(stt$EXTRA_OTOL[[VALUE]])

param.age_only
param.timeStrata
param.min_otol
param.max_otol
param.interval_otol
param.extra_otol

timeInterval<- param.timeStrata # "T" ##Time interval (T="quarter"), for the otoliths selection by length class

extra<-strsplit(param.extra_otol, " ")
otolitSet<-c(seq(param.min_otol,param.max_otol, by=param.interval_otol), unlist(extra))

########
getAllData<-function(tabname){
  ress<-NULL
  for(n in otolitSet){
    tval<-read.table(paste(output_dir,tabname,n,".csv",sep=""),sep=SEP, header=T)
    ress<-rbind(ress, tval)
  }
  return(ress)
}


#### Data of age, length by year (from the individuals selected in each simulation run)

data_selected_lt_age<-getAllData("dados_lt_age_")
summary(data_selected_lt_age)

############################################################################
####   von Bertalanffy growth model parameters by year and simulation run
############################################################################
#############################################################################

simulvbgm_param_year<-getAllData("results_simulvbgm_")
summary(simulvbgm_param_year)

####For all the years, to compare the VGBGM parameters between years
####Figure 7 - Summary of parameters by year for the full set of simulations (n=100), for j's (number of selected otoliths)
Fig7_K_VBGM <- file.path(paste(output_dir,"Fig7_K_VBGM_", timeInterval, ".png", sep = ""))
png(file=Fig7_K_VBGM)
fig7_K<-ggplot(simulvbgm_param_year, aes(x=factor(type), y=K, fill=factor(type))) +
  geom_boxplot(outlier.shape=NA)+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig7_K)
dev.off()

Fig7_t0_VBGM <- file.path(paste(output_dir,"Fig7_t0_VBGM_", timeInterval, ".png", sep = ""))
png(file=Fig7_t0_VBGM)
fig7_t0<-ggplot(simulvbgm_param_year, aes(x=factor(type), y=t0,  fill=factor(type))) +
  geom_boxplot(outlier.shape=NA)+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig7_t0)
dev.off()

Fig7_Linf_VBGM <- file.path(paste(output_dir,"Fig7_Linf_VBGM_", timeInterval, ".png", sep = ""))
png(file=Fig7_Linf_VBGM)
fig7_Linf<-ggplot(simulvbgm_param_year, aes(x=factor(type), y=Linf, fill=factor(type))) +
  geom_boxplot(outlier.shape=NA)+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+ ylim(c(20,100))+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig7_Linf)
dev.off()



####################################################################################
#####  Mean length at age data original and by simulation run
###################################################################################
####################################################################################

table_mla_sims<-getAllData("table_original_simulsub_mla_")
summary(table_mla_sims)

###############################################################################################################
###### Figure 8 - compare mean length at age (oringinal versus simulation) by year and simulation run
#######
#timeInterval<- "T"
years<-unique(table_mla_sims$year)

for(i in years)
{
  fig8_mla<- ggplot(data=subset(table_mla_sims, year==i),aes(x=factor(age), y=m_lt,fill=factor(type))) +
    geom_boxplot(outlier.shape = NA)+xlab("age")+ylab("mean length at age (cm)")+scale_color_brewer(palette = "Paired")+ theme_classic()+
    facet_wrap(~year)+ggtitle(i)+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  ggsave(fig8_mla, file=paste0(output_sub_dir,"Fig8_mla_",i,timeInterval,".png"),width=14, height = 10, units="cm")
}



###########################################################################################
### Standard deviation from length at age by simualtion run
#####################################################################
#############################################################################################


table_sdla_sims<-getAllData("table_original_simulsub_sd_")
summary(table_sdla_sims)

################################################################################################
###### Figure 9 - compare the sd length at age (oringinal versus simulation) by year and simulation run
#######
#table_sdla_sims<-table_sdla_sims[!is.na(table_sdla_sims$sd_lt),]##remove NAs on data[table_sdla_sims$year==years[[nb]],]

years<-unique(table_sdla_sims$year)

for(i in years)
{
  fig9_sdla<- ggplot(data=subset(table_sdla_sims, year==i),aes(x=factor(age), y=sd_lt,fill=factor(type))) +
    geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("sd length at age (cm)")+scale_color_brewer(palette = "Paired")+ theme_classic()+
    facet_wrap(~year)+ggtitle(i)+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  ggsave(fig9_sdla, file=paste0(output_sub_dir,"Fig9_sdla_",i,timeInterval,".png"),width=14, height = 10, units="cm")
}


##################################################################################################
###
### Data combine (table_mla_sims, table_sdla_sims)
##################################################################################################

data_mla_sd<- cbind(table_mla_sims, table_sdla_sims)
data_mla_sd<- data_mla_sd[,c(1,2,3,9,10,11,12)]
data_mla_sd<-data_mla_sd[complete.cases(data_mla_sd$sd_lt),]
data_mla_sd$coefv<-data_mla_sd$m_lt/data_mla_sd$sd_lt


years<-unique(data_mla_sd$year)

for(i in years)
{
  fig9_cvla<- ggplot(data=subset(data_mla_sd, year==i),aes(x=factor(age), y=coefv,fill=factor(type))) +
    geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("CV (lenght)")+scale_color_brewer(palette = "Paired")+ theme_classic()+
    facet_wrap(~year)+ggtitle(i)+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  ggsave(fig9_cvla, file=paste0(output_sub_dir,"Fig9a_cvla_",i,timeInterval,".png"),width=14, height = 10, units="cm")
}


####################################################################################################
### Stats (mape, rmspe, mspe) from each simualations run and year
######################################################################
######################################################################################################

stats_simul<-getAllData("table_res_stat_")
summary(stats_simul)

###############################################################################################################
###### Figure 10 - compare the stats (mape, mspe, rtmspe) by year and by simulation type (number of otoliths/length class)
#######

Fig10_mspe <- file.path(paste(output_dir,"Fig10_mspe_", timeInterval, ".png", sep = ""))
png(file=Fig10_mspe)
fig10_mspe<-ggplot(stats_simul, aes(x=factor(type), y=mspe, group=1)) + geom_step()+
  #geom_point(size=2)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig10_mspe)
dev.off()

Fig10_mape <- file.path(paste(output_dir,"Fig10_mape_", timeInterval, ".png", sep = ""))
png(file=Fig10_mape)
fig10_mape<-ggplot(stats_simul, aes(x=factor(type), y=mape, group=1)) + geom_step()+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig10_mape)
dev.off()

Fig10_rtmspe <- file.path(paste(output_dir,"Fig10_rtmspe_", timeInterval, ".png", sep = ""))
png(file=Fig10_rtmspe)
fig10_rtmspe<-ggplot(stats_simul, aes(x=factor(type), y=rtmspe, group=1)) + geom_step()+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig10_rtmspe)
dev.off()


###### Normalized data stats


##Mean by stats by year
mean_mspe <- stats_simul %>% group_by(year) %>% summarize(Mean = mean(mspe, na.rm=TRUE))
mean_mape <- stats_simul %>% group_by(year) %>% summarize(Mean = mean(mape, na.rm=TRUE))
mean_rtmspe <- stats_simul %>% group_by(year) %>% summarize(Mean = mean(rtmspe, na.rm=TRUE))

### Normalization of data by stats and year (max)
nor_mspe<- stats_simul %>% group_by(year) %>% mutate(norm_mspe = mspe/max(mspe))
nor_mape<- stats_simul %>% group_by(year) %>% mutate(norm_mape = mape/max(mape))
nor_rtmspe<- stats_simul %>% group_by(year) %>% mutate(norm_rtmspe = rtmspe/max(rtmspe))


Fig10a_norm_mspe <- file.path(paste(output_dir,"Fig10_norm_mspe_", timeInterval, ".png", sep = ""))
png(file=Fig10a_norm_mspe)
fig10a_mspe<-ggplot(nor_mspe, aes(x=factor(type), y=norm_mspe, group=1)) + geom_step()+
  #geom_point(size=2)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig10a_mspe)
dev.off()

Fig10a_norm_mape <- file.path(paste(output_dir,"Fig10_norm_mape_", timeInterval, ".png", sep = ""))
png(file=Fig10a_norm_mape)
fig10a_mape<-ggplot(nor_mape, aes(x=factor(type), y=norm_mape, group=1)) + geom_step(stat="summary", fun.y = mean)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig10a_mape)
dev.off()

Fig10a_norm_rtmspe <- file.path(paste(output_dir,"Fig10_norm_rtmspe_", timeInterval, ".png", sep = ""))
### tinha isto png(file=Fig10a_norm_rtimeIntervalspe)
png(file=Fig10a_norm_rtmspe)
fig10a_rtmspe<-ggplot(nor_rtmspe, aes(x=factor(type), y=norm_rtmspe, group=1)) + geom_step(stat="summary", fun.y = mean)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig10a_rtmspe)
dev.off()

### Normalization of data by stats and year (mean)
nor_mn_mspe<- stats_simul %>% group_by(year) %>% mutate(norm_mn_mspe = mspe/mean(mspe))
nor_mn_mape<- stats_simul %>% group_by(year) %>% mutate(norm_mn_mape = mape/mean(mape))
nor_mn_rtmspe<- stats_simul %>% group_by(year) %>% mutate(norm_mn_rtmspe = rtmspe/mean(rtmspe))


Fig10a_norm_mn_mspe <- file.path(paste(output_dir,"Fig10_norm_mn_mspe_", timeInterval, ".png", sep = ""))
png(file=Fig10a_norm_mn_mspe)
fig10a_mn_mspe<-ggplot(nor_mn_mspe, aes(x=factor(type), y=norm_mn_mspe, group=1)) + geom_step()+
  #geom_point(size=2)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig10a_mn_mspe)
dev.off()

Fig10a_norm_mn_mape <- file.path(paste(output_dir,"Fig10_norm_mn_mape_", timeInterval, ".png", sep = ""))
png(file=Fig10a_norm_mn_mape)
fig10a_mn_mape<-ggplot(nor_mn_mape, aes(x=factor(type), y=norm_mn_mape, group=1)) + geom_step(stat="summary", fun.y = mean)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig10a_mn_mape)
dev.off()

Fig10a_norm_mn_rtmspe <- file.path(paste(output_dir,"Fig10_norm_mn_rtmspe_", timeInterval, ".png", sep = ""))
png(file=Fig10a_norm_mn_rtmspe)
fig10a_mn_rtmspe<-ggplot(nor_mn_rtmspe, aes(x=factor(type), y=norm_mn_rtmspe, group=1)) + geom_step(stat="summary", fun.y = mean)+
  xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig10a_mn_rtmspe)
dev.off()



###################################################################################################################
####### Comparison of von Bertalanffy growth model parameters estimated by year
###### and for each j (number of otoliths by length class) the result of the 100 simulations aggregated
###################################################################################################################

Fig11_K_VBGM <- file.path(paste(output_dir,"Fig11_K_VBGM_", timeInterval, ".png", sep = ""))
png(file=Fig11_K_VBGM)
fig11_K<-ggplot(stats_simul, aes(x=factor(type), y=k)) +
  geom_point(size=2,colour="green")+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig11_K)
dev.off()

Fig11_t0_VBGM <- file.path(paste(output_dir,"Fig11_t0_VBGM_", timeInterval, ".png", sep = ""))
png(file=Fig11_t0_VBGM)
fig11_t0<-ggplot(stats_simul, aes(x=factor(type), y=t0)) +
  geom_point(size=2,colour="green")+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig11_t0)
dev.off()

Fig11_Linf_VBGM <- file.path(paste(output_dir,"Fig11_Linf_VBGM_", timeInterval, ".png", sep = ""))
png(file=Fig11_Linf_VBGM)
fig11_Linf<-ggplot(stats_simul, aes(x=factor(type), y=Linf)) +
  geom_point(size=2,colour="green")+xlab("number of otoliths selected by length class (cm)")+theme_classic()+
  facet_wrap(~year)+#ylim(c(20,100))+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig11_Linf)
dev.off()



################################################################################################################
###    Based on the VBGM by sim and year, the length for each selected fish was predicted based on age data
###########################################################################################################
#################################################################################################################

lpredict_vbsim<-getAllData("vb_predict_melt_")
summary(lpredict_vbsim)

###############################################################################################################
###### Figure 12 - compare the stats (mape, mspe, rtmspe) by year and by simulation type (number of otoliths/length class)
#######
Fig12_predictLt <- file.path(paste(output_dir,"Fig12_predictLt_", timeInterval, ".png", sep = ""))
png(file=Fig12_predictLt)
fig12_predictLt<-ggplot(data=lpredict_vbsim, aes(x=ID_ind, y=pred_lt,colour=type)) +
  geom_line()+xlab("ID_ind")+ylab("Predicted length")+theme_classic()+
  facet_wrap(~year)+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig12_predictLt)
dev.off()


Fig12a_predictLt <- file.path(paste(output_dir,"Fig12a_predictLt_", timeInterval, ".png", sep = ""))
png(file=Fig12a_predictLt)
fig12a_predictLt<-ggplot(data=lpredict_vbsim, aes(x=ID_ind, y=pred_lt,colour=factor(year))) +
  geom_point()+xlab("ID_ind")+ylab("Predicted length")+theme_classic()+
  facet_wrap(~factor(type))+
  theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
        axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
print(fig12a_predictLt)
dev.off()



  ###############################################################################################################################
  ###############################################################################################################################
  ### Data from biological samples
  ###############################################################################################################################
  ###############################################################################################################################

  dados_bio_simul<-getAllData("dados_bio_")
  summary(dados_bio_simul)

 ###### MEAN WEIGHT at age
  dados_bio_simul<- dados_bio_simul[complete.cases(dados_bio_simul$wt),] ##to remove the NAs

  table_meanweight<-group_by(dados_bio_simul, age, year, ID_sim,type) %>% summarize(m_wg = mean(wt))
  years<-unique(table_meanweight$year)

  for(i in years)
  {
    Figweight<- ggplot(data=subset(table_meanweight, year==i),aes(x=factor(age), y=m_wg,fill=factor(type))) +
      geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("mean weight at age")+scale_color_brewer(palette = "Paired")+ theme_classic()+
      facet_wrap(~year)+ggtitle(i)+
      theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
            axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
    ggsave(Figweight, file=paste0(output_sub_dir,"Figweight",i,timeInterval,".png"),width=14, height = 10, units="cm")
  }


  if(timeInterval=="T"){

    table_meanweight<-group_by(dados_bio_simul, age, year, quarter, ID_sim,type) %>% summarize(m_wg = mean(wt))

  for(i in years)
  {
    Figweight_quarter<- ggplot(data=subset(table_meanweight, year==i),aes(x=factor(age), y=m_wg,fill=factor(type))) +
      geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("mean weight at age")+scale_color_brewer(palette = "Paired")+ theme_classic()+
      facet_grid(quarter~year)+ggtitle(i)+
      theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
            axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
    ggsave(Figweight_quarter, file=paste0(output_sub_dir,"Figweight_quarter",i,timeInterval,".png"),width=14, height = 10, units="cm")
  }

}

  ###### MEAN WEIGHT at length
  table_meanweight_Lt<-group_by(dados_bio_simul, Lt, year, ID_sim,type) %>% summarize(m_wg = mean(wt))
  years<-unique(table_meanweight_Lt$year)

  for(i in years)
  {
    Figweight_lt<- ggplot(data=subset(table_meanweight_Lt, year==i),aes(x=as.factor(Lt), y=m_wg,fill=factor(type))) +
      geom_boxplot(outlier.shape =NA)+xlab("length")+ylab("mean weight at length")+scale_color_brewer(palette = "Paired")+ theme_classic()+
      facet_wrap(~year)+ggtitle(i)+
      theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
            axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
    ggsave(Figweight_lt, file=paste0(output_sub_dir,"Figweight_lt",i,timeInterval,".png"),width=14, height = 10, units="cm")
  }

  if(timeInterval=="T"){

    table_meanweight_Lt_quarter<-group_by(dados_bio_simul, Lt, year, quarter, ID_sim,type) %>% summarize(m_wg = mean(wt))
    write.csv(table_meanweight_Lt_quarter,"output/simulation_results/table_meanweight_Lt_quarter.csv")

    for(i in years)
    {
      Figweight_ltquarter<- ggplot(data=subset(table_meanweight_Lt_quarter, year==i),aes(x=as.factor(Lt), y=m_wg,fill=factor(type))) +
        geom_boxplot(outlier.shape =NA)+xlab("length")+ylab("mean weight at length")+scale_color_brewer(palette = "Paired")+ theme_classic()+
        facet_grid(quarter~year)+ggtitle(i)+
        theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
              axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
      ggsave(Figweight_ltquarter, file=paste0(output_sub_dir,"Figweight_ltquarter",i,timeInterval,".png"),width=14, height = 10, units="cm")
    }
  }

  ###### SD WEIGHT at age
  table_sdweight<-dados_bio_simul  %>% group_by(age, year, ID_sim,type) %>% summarise(sdwt = sd(wt))
  summary(table_sdweight)
  years<-unique(table_sdweight$year)

  for(i in years)
  {
    Figweightsd<- ggplot(data=subset(table_sdweight, year==i),aes(x=factor(age), y=sdwt,fill=factor(type))) +
      geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("SD weight at age")+scale_color_brewer(palette = "Paired")+ theme_classic()+
      facet_wrap(~year)+ggtitle(i)+
      theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
            axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
    ggsave(Figweightsd, file=paste0(output_sub_dir,"Figweightsd",i,timeInterval,".png"),width=14, height = 10, units="cm")
  }

if(timeInterval=="T"){

  table_sdweight_quarter<-dados_bio_simul  %>% group_by(age, year, quarter, ID_sim,type) %>% summarise(sdwt = sd(wt))
  write.csv(table_sdweight_quarter,"output/simulation_results/table_sdweight_quarter.csv")
  for(i in years)
  {
    Figweightsd_quarter<- ggplot(data=subset(table_sdweight_quarter, year==i),aes(x=factor(age), y=sdwt,fill=factor(type))) +
      geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("SD weight at age")+scale_color_brewer(palette = "Paired")+ theme_classic()+
      facet_grid(quarter~year)+ggtitle(i)+
      theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
            axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
    ggsave(Figweightsd_quarter, file=paste0(output_sub_dir,"Figweightsd_quarter",i,timeInterval,".png"),width=14, height = 10, units="cm")
  }
}


  ###### CV WEIGHT at age
  table_cvweight<-group_by(dados_bio_simul, age, year, ID_sim,type) %>% summarize(CVwg = cv(wt))
  years<-unique(table_cvweight$year)

  for(i in years)
  {
    FigweightCV<- ggplot(data=subset(table_cvweight, year==i),aes(x=factor(age), y=CVwg,fill=factor(type))) +
      geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("CV weight at age")+scale_color_brewer(palette = "Paired")+ theme_classic()+
      facet_wrap(~year)+ggtitle(i)+
      theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
            axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
    ggsave(FigweightCV, file=paste0(output_sub_dir,"FigweightCV",i,timeInterval,".png"),width=14, height = 10, units="cm")
  }


  if(timeInterval=="T"){
    table_cvweight_quarter<-group_by(dados_bio_simul, age, year, quarter, ID_sim,type) %>% summarize(CVwg = cv(wt))
    write.csv(table_cvweight_quarter,"output/simulation_results/table_cvweight_quarter.csv")

    for(i in years)
    {
      FigweightCV_quarter<- ggplot(data=subset(table_cvweight_quarter, year==i),aes(x=factor(age), y=CVwg,fill=factor(type))) +
        geom_boxplot(outlier.shape =NA)+xlab("age")+ylab("CV weight at age")+scale_color_brewer(palette = "Paired")+ theme_classic()+
        facet_grid(quarter~year)+ggtitle(i)+
        theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
              axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
      ggsave(FigweightCV_quarter, file=paste0(output_sub_dir,"FigweightCV_quarter",i,timeInterval,".png"),width=14, height = 10, units="cm")
    }
  }

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

