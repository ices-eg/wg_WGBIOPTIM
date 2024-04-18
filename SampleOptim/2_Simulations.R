#######################################################################################################
#######################################################################################################
#######################################################################################################
##
##
##   Biological sampling optimization (Script "SampleOptim")
##   Developed by: Patricia Goncalves (patricia@ipma.pt)
##   Last version development period: may 2023
##   Version: v4.1
##   Updated: 18 april 2024
##
##   Reference:
##   Gonçalves, Patrícia 2021. "SampleOptim" a data analysis R-tool to optimize fish sampling for
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



########################################################################################################
#### Files path and function source:
#setwd("~/...") ##Set directory

setwd("C:\\Users\\gonpa\\Documents\\WGBIOPTIM\\04.SampleOptim\\SampleOptim") ##Set directory
output_dir<-"output/"
output_dir<-"output/"
#dir.create(output_dir, showWarnings = FALSE) #this folder should have been created on R-script 1.

###Biological sample data (Applied to a period of years)
data<-read.table("SampleOptim_input_fileformat_example.csv",sep=";", header=T)

data_samplebio<- data %>%
  rename(
    ID_BIO_PEIXES= ID_BIO_FISH, ### this column needs to be created corresponds to the number of the line
    DATA = date,
    MONTH = Month,
    ANO = Year,
    COD_FAO = Species,
    PORTO_NOME = Port,
    C_CLASSE = LengthClass,
    P_INDIVIDUAL = Weight,
    SEXO = Sex,
    EST_MATURACAO = MaturityStage,
    IDADE = Age)


SEP<-","
VERBOSE<-F


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

param.species<-stt$species[[VALUE]]
param.area<-stt$AREA[[VALUE]]
param.age_only<-as.logical(stt$AGE_ONLY[[VALUE]])
param.port<-as.logical(stt$PORT[[VALUE]])
param.distUniPort<-as.logical(stt$distUniPorto[[VALUE]])
param.timeStrata<-stt$TIME_STRATA[[VALUE]]
param.sex_ratio<-as.numeric(stt$SEX_RATIO[[VALUE]])
param.min_lc<-as.numeric(stt$MIN_LC[[VALUE]])
param.max_lc<-as.numeric(stt$MAX_LC[[VALUE]])
param.interval_lc<-as.numeric(stt$interval_LC[[VALUE]])
param.min_age<-as.numeric(stt$MIN_age[[VALUE]])
param.max_age<-as.numeric(stt$MAX_age[[VALUE]])
param.min_otol<-as.numeric(stt$MIN_OTOL[[VALUE]])
param.max_otol<-as.numeric(stt$MAX_OTOL[[VALUE]])
param.interval_otol<-as.numeric(stt$interval_OTOL[[VALUE]])
param.extra_otol<-as.array(stt$EXTRA_OTOL[[VALUE]])
param.Linf<-as.numeric(stt$Linf[[VALUE]])
param.K<-as.numeric(stt$K[[VALUE]])
param.t0<-as.numeric(stt$t0[[VALUE]])
param.year_start<-as.numeric(stt$year_start[[VALUE]])
param.year_end<-as.numeric(stt$year_end[[VALUE]])
param.stage_mature<-as.numeric(stt$stage_mature[[VALUE]])
param.numSubSamples<-as.numeric(stt$n[[VALUE]])

## Check the input values for the simulations
param.species
param.area
param.age_only
param.port
param.distUniPort
param.timeStrata
param.sex_ratio
param.min_lc
param.max_lc
param.interval_lc
param.min_age
param.max_age
param.min_otol
param.max_otol
param.interval_otol
param.extra_otol
param.Linf
param.K
param.t0
param.year_start
param.year_end
param.age_only
param.numSubSamples

if(
  (param.min_lc > param.max_lc) |
  (param.min_age > param.max_age) |
  (param.min_otol > param.max_otol) |
  (param.year_start > param.year_end) |
  (param.numSubSamples <=0)
){
  print("ERROR in input parameter table")
}

extra<-strsplit(param.extra_otol, " ")
otolitSet<-as.numeric( c(seq(param.min_otol,param.max_otol, by=param.interval_otol), unlist(extra)) )

class_length_set<- seq(param.min_lc,param.max_lc,param.interval_lc); #set length class range and the intervals
Linf<-param.Linf #50
K<-param.K #0.1
t0<-param.t0 #-3
timeInterval<- param.timeStrata # "T" ##Time interval (T="quarter"), for the otoliths selection by length class
sexRatio<-param.sex_ratio #0.5
numSubSamples<-param.numSubSamples #n
newdata<-seq(param.min_age,param.max_age,0.1) ###set age distribution vector for predictions
newdata<-data.frame(newdata)
colnames(newdata)<- "IDADE"
year_init<- param.year_start
year_last<- param.year_end
yearSet<- c(year_init:year_last)

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

#### 2.1 Data preparation for run the simulations:
###REMOVE NAs in age data (Note: I only use the individuals where age has atributted)
### Age classes
data_samplebio<-data_samplebio[c("ID_BIO_PEIXES", "MONTH", "ANO", "PORTO_NOME", "IDADE",
                      "COD_FAO", "EST_MATURACAO", "C_CLASSE", "SEXO", "P_INDIVIDUAL")]

data_samplebio<-data_samplebio[!is.na(data_samplebio$IDADE),]
data_samplebio<-data_samplebio[!is.na(data_samplebio$MONTH),]

if(dim(data_samplebio[is.na(data_samplebio$IDADE),])[1]>0) message("ERRO-> IDADE COM NA")
if(dim(data_samplebio[is.na(data_samplebio$MONTH),])[1]>0) message("ERRO-> MONTH COM NA")
if(dim(data_samplebio[is.na(data_samplebio$ID_BIO_PEIXES),])[1]>0) message("ERRO-> ID_BIO_PEIXES COM NA")
if(dim(data_samplebio[is.na(data_samplebio$C_CLASSE),])[1]>0) message("ERRO-> C_CLASSE COM NA")

summary(data_samplebio)

valid_years<-unique(data_samplebio$ANO)
anos<-intersect(yearSet, valid_years)



#class_length_set<- seq(param.min_lc,param.max_lc,param.interval_lc); #set length class range and the intervals

source("sample_selection_function.R") ##Function for randomly select the samples by length class

set.seed(2019) #any value

########################################################################################################
########################################################################################################
########################################################################################################
####2.2  Set function for von Bertalanffy growth model parameters (Linf, K and t0)
########################################################################################################

vonberPorAno<-function(ano){
  print(paste("Vonber: ",ano,sep=""))
  system.time(
    vvv <- lapply(1:numSubSamples, function(i,ano) {
      print(paste("V:",ano,":",numOtolitsPerClass,":",i,sep=""))
      dataV<- amostraTemporal(data_samplebio[data_samplebio$ANO==ano,], numOtolitsPerClass, class_length_set, sexRatio,timeInterval, porto=FALSE, distUniPorto=FALSE, verbose=VERBOSE)
      #print(dim(dataV))
      ##Retira a amostra definida
      fitTypical <-NULL
      salta<-FALSE  # só um <
      tryCatch(
        {
          fitTypical <- nls(vbTypical, data=dataV, start=svTypical, control)
        },
        error = function(e) { print(e); salta<<-TRUE }, ## caso n?o tenha convergido dois <<
        warning = function(w) { print(w); salta<<-TRUE },
        message = function(m) { print(m); salta<<-TRUE },
        finally = function(m) { print("vai saltar") }
      )
      if(salta | is.null(fitTypical)){
        return(list(ANO=ano,data=dataV,coef=NA,Linf=NA,K=NA,t0=NA,predict_values=NA,n=i,j=numOtolitsPerClass))
      }
      coef<-summary(fitTypical)$coefficients
      Linf<-summary(fitTypical)$coefficients[[1]]
      K<-summary(fitTypical)$coefficients[[2]]
      t0<-summary(fitTypical)$coefficients[[3]]
      predict_values<- predict(fitTypical,newdata=newdata)
      return(list(ANO=ano,data=dataV,coef=coef,Linf=Linf,K=K,t0=t0,predict_values=predict_values,n=i,j=numOtolitsPerClass))
    }, ano)
  )[[3]] #
  return(vvv)
}


#########################################################################################################
#########################################################################################################
#########################################################################################################
####            **  START SIMULATIONS  **    ###########################################################
#########################################################################################################
### Note: repeat the code between lines 100 and 388 for each "numOtolitsPerClass" (number of otoliths by length class)
#########################################################################################################
#########################################################################################################
### Quarter/Sample
## SR=1:1
# numOtolitsPerClass=1  (Conditions:porto=FALSE, distUniPorto=FALSE)
# numOtolitsPerClass - is the number of otoliths selected by length class
#########################################################################################################
#########################################################################################################
##
#numOtolitsPerClass<-3
#otolitSet<-c(1:10,15,20)
for(numOtolitsPerClass in otolitSet){
  ##### Initial values set and fixed by species
  ###Initial values for parameters to all the simulations
  
  svTypical <- list(Linf=Linf,K=K,t0=t0) ##Initial parameters values for the growth curve
  vbTypical <- C_CLASSE~Linf*(1-exp(-K*(IDADE-t0))) ##von Bertallanfy growth model
 #control<- nls.control(maxiter=10000, minFactor = 1/4096, printEval=F)
  control<- nls.control(maxiter=10000)
 
  vonber<-sapply(anos, vonberPorAno)


  ######################################
  ######################################
  #extracting the variables from each of the 100 samples (subsamples)
  #
  vb_ano_x100<- sapply(vonber, function(x) x$ANO)
  vb_ano<- sapply(vonber, function(x) x$data$ANO)
  trocas<-typeof(vb_ano)!="list"

  vb_idade<- sapply(vonber, function(x) x$data$IDADE)
  #if(typeof(vb_idade)!="list") vb_idade<-as.list(vb_idade)

  vb_ct<- sapply(vonber, function(x) x$data$C_CLASSE)
  #if(typeof(vb_ct)!="list") vb_ct<-as.list(vb_ct)

  vb_linf<- as.list(sapply(vonber, function(x) x$Linf))
  vb_k<- as.list(sapply(vonber, function(x) x$K))
  vb_t0<- as.list(sapply(vonber, function(x) x$t0))
  vb_n<- as.list(sapply(vonber, function(x) x$n))
  vb_sex<- sapply(vonber, function(x) x$data$SEXO)
  vb_matur<- sapply(vonber, function(x) x$data$EST_MATURACAO)
  vb_weight<- sapply(vonber, function(x) x$data$P_INDIVIDUAL)
  vb_month<- sapply(vonber, function(x) x$data$MONTH)

  ###construir a matrix com os dados das vari?veis de cada 1 das 100 amostras(sub-amostras)
  # ##Dados de ct, idade
  vb_ano_melt<- melt(vb_ano)
  vb_idade_melt <- melt(vb_idade)
  vb_ct_melt<- melt(vb_ct)
  vb_month_melt<- melt(vb_month)

  if(trocas){

    lixo_a3<-vb_ano_melt
    lixo_i3<-vb_idade_melt
    lixo_c3<-vb_ct_melt
    lixo_c4<-vb_month_melt

    lixo_aa3<-cbind(lixo_a3$value, lixo_a3$X2)
    lixo_ii3<-cbind(lixo_i3$value, lixo_i3$X2)
    lixo_cc3<-cbind(lixo_c3$value, lixo_c3$X2)
    lixo_cc4<-cbind(lixo_c4$value, lixo_c4$X2)
    vb_ano_melt<-as.data.frame(lixo_aa3)
    vb_idade_melt<-as.data.frame(lixo_ii3)
    vb_ct_melt<-as.data.frame(lixo_cc3)
    vb_month_melt<-as.data.frame(lixo_cc4)

  }

  dados_lt_age<-cbind(vb_ct_melt,vb_idade_melt,vb_ano_melt,vb_month_melt)
  dados_lt_age<-dados_lt_age[,c(1,3,5,7,8)]
  colnames(dados_lt_age)<-c("Lt","age","year","month","ID_sim")
  dados_lt_age$quarter<-factor(NA,levels=c("1","2","3","4"))
  dados_lt_age[dados_lt_age$month<=3,"quarter"]<-"1"
  dados_lt_age[dados_lt_age$month>3 & dados_lt_age$month<=6,"quarter"]<-"2"
  dados_lt_age[dados_lt_age$month>6 & dados_lt_age$month<=9,"quarter"]<-"3"
  dados_lt_age[dados_lt_age$month>9, "quarter"]<-"4"
  dados_lt_age$type<-numOtolitsPerClass


  write.table(dados_lt_age,paste(output_dir, "dados_lt_age_",numOtolitsPerClass,".csv",sep=""),sep=SEP, row.names=FALSE)

  ### Predict
  total<-array()
  for(x in vonber){
    total<-rbind(total,cbind(c(1:length(x$predict_values)), x$n, x$predict_values, x$j, x$ANO))
  }
  total<-total[-1,]

  colnames(total)<-c("ID_ind","ID_sim","pred_lt","type","year")

  vb_predict_melt<-total
  write.table(vb_predict_melt,paste(output_dir,"vb_predict_melt_",numOtolitsPerClass,".csv",sep=""),sep=SEP, row.names=FALSE)

## Bio data
  vb_wt_melt<-melt(vb_weight) #fish total weight
  vb_wt_melt<-vb_wt_melt["value"]
  vb_month_melt<-vb_month_melt[,1:2]
  vb_ct_melt<-vb_ct_melt[1]
  vb_idade_melt<-vb_idade_melt[1]
  vb_ano_melt<-vb_ano_melt[1]
  dados_bio<-cbind(vb_ct_melt,vb_idade_melt,vb_ano_melt,vb_wt_melt,vb_month_melt)
  colnames(dados_bio)<-c("Lt","age","year","wt","month","ID_sim")
  dados_bio$quarter<-factor(NA,levels=c("1","2","3","4"))
  dados_bio[dados_bio$month<=3,"quarter"]<-"1"
  dados_bio[dados_bio$month>3 & dados_bio$month<=6,"quarter"]<-"2"
  dados_bio[dados_bio$month>6 & dados_bio$month<=9,"quarter"]<-"3"
  dados_bio[dados_bio$month>9, "quarter"]<-"4"
  dados_bio$type<-numOtolitsPerClass
  write.table(dados_bio,paste(output_dir, "dados_bio_",numOtolitsPerClass,".csv",sep=""),sep=SEP, row.names=FALSE)


#####
  if(param.age_only==FALSE){
    ##Data of length, age, sex, maturity stage, weight, month and year
    vb_ct_melt<-vb_ct_melt[1] ## length
    vb_idade_melt<-vb_idade_melt[1] ##age
    vb_ano_melt<-vb_ano_melt[1]  ##year
    vb_month_melt<-vb_month_melt[,1:2]
    vb_sex_melt<-melt(vb_sex) #sex
    vb_sex_melt<-vb_sex_melt["value"]
    vb_mat_melt<-melt(vb_matur) ##maturity stage
    vb_mat_melt<-vb_mat_melt["value"]
    dados_bio<-cbind(vb_ct_melt,vb_idade_melt,vb_ano_melt,vb_sex_melt,vb_mat_melt,vb_wt_melt,vb_month_melt)
    colnames(dados_bio)<-c("Lt","age","year","sex","mat_stg","wt","month","ID_sim")

    dados_bio$mat_stg<-as.numeric(dados_bio$mat_stg)
    dados_bio$maturity<-ifelse(dados_bio$mat_stg<param.stage_mature | is.na(dados_bio$mat_stg) , 0, ifelse(dados_bio$mat_stg>=param.stage_mature,1,NA))

    #dados_bio$quarter<-ifelse(dados_bio$month<6, 1, ifelse(dados_bio$mat_stg>1,1,NA))
    dados_bio$quarter<-factor(NA,levels=c("1","2","3","4"))
    dados_bio[dados_bio$month<=3,"quarter"]<-"1"
    dados_bio[dados_bio$month>3 & dados_bio$month<=6,"quarter"]<-"2"
    dados_bio[dados_bio$month>6 & dados_bio$month<=9,"quarter"]<-"3"
    dados_bio[dados_bio$month>9, "quarter"]<-"4"
    dados_bio$type<-numOtolitsPerClass
    write.table(dados_bio,paste(output_dir, "dados_bio_",numOtolitsPerClass,".csv",sep=""),sep=SEP, row.names=FALSE)


    ###############################################################################
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### MATURITY OGIVE
    ### DETERMINE: L25, L50, L75
    ## Confindence intervals by year
    ##
    ##
    ## NOTA: SUBSET 1? QUARTER (SPAWNING SEASON)
    ##############################################################################
    ##############################################################################

    years<-unique(dados_bio$year)

    table_mature<- function (data=dados_bio){
      sim<-unique(data$ID_sim)
      results<-matrix(nrow=length(sim),ncol=6)
      for(nb in 1: length(sim))
      {
        glm1 <- glm(factor(maturity)~Lt,family=binomial,data=data[data$ID_sim==nb,])
        Lmat <- signif(dose.p (glm1, p = c(0.25, 0.50, 0.75)), digits = 3)
        results[nb,1]<-unique(data$year[data$ID_sim==nb]) #year
        results[nb,2]<-as.numeric(Lmat[[1]]) #L25
        results[nb,3]<-as.numeric(Lmat[[2]]) #L50
        results[nb,4]<-as.numeric(Lmat[[3]]) #L75
        results[nb,5]<-nb ##ID_simul
        results[nb,6]<-unique(data$type) ##numOtolitsPerClass
      }
      colnames(results)<-c("year","L25","L50","L75","ID_sim","type")
      results
    }

    #
    #table_mo<-table_mature(data=dados_bio[dados_bio$quarter=="1",]) ##Only quarter 1

    table_mo<-table_mature(data=dados_bio) ##Data from the whole year
    write.table(table_mo,paste(output_dir, "table_res_mo_",numOtolitsPerClass,".csv",sep=""),sep=SEP, row.names=FALSE)
  }

  ####Figure 3 - compare length and age distributions by year (by simulations)
  ##Note: the length distribution did not change by simulations
  ##(because the selection is based in the number of otoliths by length class)
  ################ Data from simulations - figures (length and age distribution)

  Fig3_length <- file.path(paste(output_dir, "Fig3_length_", numOtolitsPerClass,"_",timeInterval, ".png", sep = ""))
  png(file=Fig3_length)
  Fig3_length<- ggplot(dados_lt_age, aes(x=Lt, colour=factor(ID_sim))) +
                geom_density(show.legend = FALSE)+facet_wrap(~ year, ncol=2)+theme_classic()
  print(Fig3_length)
  dev.off()


  Fig3_age <- file.path(paste(output_dir, "Fig3_age_", numOtolitsPerClass,"_",timeInterval, ".png", sep = ""))
  png(file=Fig3_age)
  Fig3_age<- ggplot(dados_lt_age, aes(x=age, colour=factor(ID_sim))) +
             geom_density(show.legend = FALSE)+facet_wrap(~ year, ncol=2)+theme_classic()
  print(Fig3_age)
  dev.off()

  ####Determine mean length at age - original data and by simulation (for each year)
  table_original<-group_by(data_samplebio, IDADE, ANO) %>% summarize(m_lt = mean(C_CLASSE))
  table_original$data<-"original"
  table_original$ID_sim<-0
  table_original$type<-numOtolitsPerClass
  colnames(table_original)<- c("age","year","m_lt","data","ID_sim","type")


  table_simul<-group_by(dados_lt_age, age, year, ID_sim) %>% summarize(m_lt = mean(Lt))
  table_simul$data<-"simulations"
  table_simul$type<-numOtolitsPerClass
  table_simul<-table_simul[,c(1,2,4,5,3,6)]###organize columns
  ######

  ## Combine data original and simulations in one table (data frame)
  table_original_simul<- merge(table_original,table_simul,all=TRUE)
  #year_simul<-c("2003","2004","2005","2008","2009","2011","2014")
  year_simul<-c(year_init:year_last)

  table_original_simulsub<-table_original_simul[table_original_simul$year %in% year_simul,]
  write.table(table_original_simulsub,paste(output_dir, "table_original_simulsub_mla_",numOtolitsPerClass,".csv",sep=""),sep=SEP, row.names=FALSE)

  #############
  ###Compare mean length at age from distributions of original data with the data from simulations
  Fig4_length <- file.path(paste(output_dir, "Fig4_length_", numOtolitsPerClass,"_",timeInterval,".png", sep = ""))
  png(file=Fig4_length)
  Fig4_length<- ggplot(table_original_simulsub, aes(x=factor(age), y=m_lt, fill=factor(type))) +
    geom_bar(stat="identity", position=position_dodge())+facet_wrap(~ year, ncol=2)+xlab("Age")+
    ylab("Mean length (cm)")+
    theme_classic()
  print(Fig4_length)
  dev.off()


  Fig4_age <- file.path(paste(output_dir, "Fig4_age_", numOtolitsPerClass,"_",timeInterval,".png", sep = ""))
  png(file=Fig4_age)
  Fig4_age<- ggplot(table_original_simulsub, aes(x=factor(age),y=m_lt,colour=factor(type))) +
    geom_boxplot()+
    facet_wrap(~ year, ncol=2)+xlab("Age")+
    ylab("Mean length (cm)")+
    theme_classic()
  print(Fig4_age)
  dev.off()


  ####Determine sd (length) at age - original data and by simulation (for each year)
  table_originalsd<-group_by(data_samplebio, IDADE, ANO) %>% summarize(sd_lt = sd(C_CLASSE))
  table_originalsd$data<-"original"
  table_originalsd$ID_sim<-0
  table_originalsd$type<-numOtolitsPerClass
  colnames(table_originalsd)<- c("age","year","sd_lt","data","ID_sim","type")

  table_simulsd<-group_by(dados_lt_age, age, year, ID_sim) %>% summarize(sd_lt = sd(Lt))
  table_simulsd$data<-"simulations"
  table_simulsd$type<-numOtolitsPerClass
  table_simulsd<-table_simulsd[,c(1,2,4,5,3,6)]###organize columns

  ######
  ## Combine data original and simulations in one table (data frame)
  table_original_simulsd<- merge(table_originalsd,table_simulsd,all=TRUE)
  #year_simul<-c("2003","2004","2005","2008","2009","2011","2012")
  table_original_simulsubsd<-table_original_simulsd[table_original_simulsd$year %in% year_simul,]
  write.table(table_original_simulsubsd,paste(output_dir, "table_original_simulsub_sd_",numOtolitsPerClass,".csv",sep=""),sep=SEP, row.names=FALSE)

  #############
  ###Compare standard deviation of length at age from distributions of original data with the data from simulations
  Fig5_length <- file.path(paste(output_dir, "Fig5_length_", numOtolitsPerClass,"_",timeInterval, ".png", sep = ""))
  png(file=Fig5_length)
  Fig5_length<- ggplot(table_original_simulsubsd, aes(x=factor(age), y=sd_lt, fill=factor(type))) +
    geom_bar(stat="identity", position=position_dodge())+facet_wrap(~ year, ncol=2)+xlab("Age")+
    ylab("Standard deviation of length (cm)")+
    theme_classic()
  print(Fig5_length)
  dev.off()


  Fig5_age <- file.path(paste(output_dir, "Fig5_age_", numOtolitsPerClass,"_",timeInterval, ".png", sep = ""))
  png(file=Fig5_age)
  Fig5_age<- ggplot(table_original_simulsubsd, aes(x=factor(age),y=sd_lt,colour=factor(type))) +
    geom_boxplot()+
    facet_wrap(~ year, ncol=2)+xlab("Age")+
    ylab("Standard deviation of length (cm)")+
    theme_classic()
  print(Fig5_age)
  dev.off()


  ################################################################
  #################################################################
  ### Growth parameters from the von Bertallanfy model by year
  vb_anounique_melt<-melt(vb_ano_x100)
  colnames(vb_anounique_melt)<- c("year")
  vb_linf_melt <- melt(vb_linf)
  colnames(vb_linf_melt)<- c("Linf","ID_sim")
  vb_k_melt<- melt(vb_k)
  colnames(vb_k_melt)<- c("k","ID_sim")
  vb_t0_melt<- melt(vb_t0)
  colnames(vb_t0_melt)<- c("t0","ID_sim")
  trimestre_simulvb<- cbind(vb_linf_melt, vb_k_melt, vb_t0_melt,vb_anounique_melt)
  trimestre_simulvb$type<-numOtolitsPerClass
  colnames(trimestre_simulvb)<- c("Linf","ID_sim","K","ID_sim","t0","ID_sim","year","type")
  trimestre_simulvb<- trimestre_simulvb[,c(1,2,3,5,7,8)]
  write.table(trimestre_simulvb,paste(output_dir, "results_simulvbgm_",numOtolitsPerClass,".csv",sep=""),sep=SEP, row.names=FALSE)


  ####For all the years, to compare the VGBGM parameters between years
  ####Figure 6 - Summary of parameters by year for the full set of simulations (n=100), by numOtolitsPerClass (number of selected otoliths)
  Fig6_K_VBGM <- file.path(paste(output_dir, "Fig6_K_VBGM_", numOtolitsPerClass,"_",timeInterval, ".png", sep = ""))
  png(file=Fig6_K_VBGM)
  fig6_K<-ggplot(trimestre_simulvb, aes(x=factor(year), y=K)) +
    geom_boxplot()+xlab("year")+theme_classic()+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  print(fig6_K)
  dev.off()


  Fig6_t0_VBGM <- file.path(paste(output_dir, "Fig6_t0_VBGM_", numOtolitsPerClass,"_",timeInterval, ".png", sep = ""))
  png(file=Fig6_t0_VBGM)
  fig6_t0<-ggplot(trimestre_simulvb, aes(x=factor(year), y=t0)) +
    geom_boxplot()+xlab("year")+theme_classic()+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  print(fig6_t0)
  dev.off()


  Fig6_Linf_VBGM <- file.path(paste(output_dir, "Fig6_Linf_VBGM_", numOtolitsPerClass,"_",timeInterval, ".png", sep = ""))
  png(file=Fig6_Linf_VBGM)
  fig6_Linf<-ggplot(trimestre_simulvb, aes(x=factor(year), y=Linf)) +
    geom_boxplot()+xlab("year")+theme_classic()+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  print(fig6_Linf)
  dev.off()

  #### Von bertallanfy growth model
  svTypical <- list(Linf=Linf,K=K,t0=t0) ##Initial growth parameters
  vbTypical <- Lt~Linf*(1-exp(-K*(age-t0))) ##von Bertallanfy growth model
  control<- nls.control(maxiter=10000)

  ###############################################################################
  ###############################################################################
  ### Determining mean square error between sim and obs (cost=mspe, mape, rtmspe)
  ### mspe(sim, obs, na.rm=TRUE)
  ##############################################################################
  ##############################################################################
  table_results<- function (data=dados_lt_age){
    years<-unique(dados_lt_age$year)
    results<-matrix(nrow=length(years),ncol=8)
    for(nb in 1: length(years)) {
      #fitTypical <- nls(vbTypical,data=data[data$year==years[nb],],start=svTypical,control)

      salta<-FALSE  # só um <
      tryCatch(
        {
          fitTypical <- nls(vbTypical,data=data[data$year==years[nb],],start=svTypical,control)
        },
        error = function(e) { print(e); salta<<-TRUE }, ## caso n?o tenha convergido dois <<
        warning = function(w) { print(w); salta<<-TRUE },
        message = function(m) { print(m); salta<<-TRUE },
        finally = function(m) { print("vai saltar") }
      )
      if(salta){
        results[nb,1]<-as.numeric(years[nb])
        results[nb,2]<-NA
        results[nb,3]<-NA
        results[nb,4]<-NA
        results[nb,5]<-NA
        results[nb,6]<-NA
        results[nb,7]<-NA
        results[nb,8]<-unique(data$type) ##numOtolitsPerClass
        next
      }


      mspe_fit<-cvFit(fitTypical,Lt ~ Linf * (1 - exp(-K * (age - t0))),data=data[data$year==years[nb],],y=data$age[data$year==years[nb]], cost=mspe, k=10)
      mape_fit<-cvFit(fitTypical,Lt~Linf*(1-exp(-K*(age-t0))),data=data[data$year==years[nb],],y=data$age[data$year==years[nb]], cost=mape, k=10)
      rtmspe_fit<-cvFit(fitTypical,Lt~Linf*(1-exp(-K*(age-t0))),data=data[data$year==years[nb],],y=data$age[data$year==years[nb]], cost=rtmspe, k=10)
      results[nb,1]<-as.numeric(years[nb])
      results[nb,2]<-as.numeric(coef(fitTypical)[1])
      results[nb,3]<-as.numeric(coef(fitTypical)[2])
      results[nb,4]<-as.numeric(coef(fitTypical)[3])
      results[nb,5]<-as.numeric(mspe_fit$cv)
      results[nb,6]<-as.numeric(mape_fit$cv)
      results[nb,7]<-as.numeric(rtmspe_fit$cv)
      results[nb,8]<-unique(data$type) ##numOtolitsPerClass
    }
    colnames(results)<-c("year","Linf","k","t0","mspe","mape","rtmspe","type")
    results
  }


  table_res<-table_results(data=dados_lt_age)
  write.table(table_res,paste(output_dir, "table_res_stat_",numOtolitsPerClass,".csv",sep=""),sep=SEP, row.names=FALSE)
}
######################################## END Simulations ##################################################
###########################################################################################################
