#######################################################################################################
#######################################################################################################
#######################################################################################################
##
##
##   Biological sampling optimization (Script "SampleOptim")
##   Developed by: Patricia Goncalves (patricia@ipma.pt)
##   Last version development period: June 2021
##   Version: v3.1
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



########################################################################################################
#### Files path and function source:
#setwd("C:\\Users\\patricia\\Documents\\WHB\\PNAB2021\\Optimizacao_amostragem_PNAB\\04.Testes_simulaPG\\WHB\\fev2019_testes") ##Set directory

setwd("C:\\Users\\gonpa\\Documents\\WGBIOPTIM\\04.SampleOptim\\SampleOptim") ##Set directory
output_dir<-"output/"
dir.create(output_dir, showWarnings = FALSE)

###Biological sample data (Applied to a period of years)
data<-read.table("SampleOptim_input_fileformat_example.csv",sep=";", header=T)

 data_samplebio<- data %>%
  rename(
     ID_BIO_PEIXES = ID_BIO_FISH, ### this column needs to be created corresponds to the number of the line
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
DATA_COL_NAME<-"data_fin"


#########################################################################################################
#### 1. Data Preliminary analysis: (Exploratory analysis)
#### Summary

summary(data_samplebio)
#names(data_samplebio)

table(data_samplebio$C_CLASSE, data_samplebio$MONTH) ##summary of the number of individuals by length class and month
plot(data_samplebio$C_CLASSE,data_samplebio$MONTH, xlab="Length class", ylab="Month of sampling")

### Number of samples by Port, year and month
nsamples_ano_mes<- data_samplebio %>% group_by(PORTO_NOME, MONTH, ANO) %>% count(paste(DATA_COL_NAME))
#write.table(nsamples_ano_mes, "numbersamples_summary_WHB.csv",sep=SEP)


#### Length classes of the samples by Port, year and month
lengthclass_samples_ano_mes<- data_samplebio %>% group_by(PORTO_NOME, MONTH, ANO) %>% count(C_CLASSE)
#write.table(lengthclass_samples_ano_mes, "numberlengthclasses_samples_summary_WHB.csv",sep=SEP)

########Figure a - length distribution samples by Port by year and month
Port<-unique(data_samplebio$PORTO_NOME) ## list of Ports names
year<- sort(unique(lengthclass_samples_ano_mes$ANO)) ##list of years on the samples data

for(bb in 1:length(year)){
  plota<- ggplot(lengthclass_samples_ano_mes[lengthclass_samples_ano_mes$ANO==year[bb],], aes(x=C_CLASSE,y=n, colour=PORTO_NOME))+xlab("length")+ ylab("number of individuals")+
    geom_line()+ theme_classic() + facet_wrap(~factor(MONTH))+
    theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
          axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
  dev.copy(png, paste(output_dir, year[bb],"_length_distribution_samples_Port_year",".png",sep=""))
  print(plota)
  dev.off()
}


#### Age of the samples by Port, year and month
age_samples_ano_mes<- data_samplebio %>% group_by(PORTO_NOME, MONTH, ANO) %>% count(IDADE)

########Figure b - age distribution samples by Port by year and month
#Port<-unique(data_samplebio$PORTO_NOME) ## list of Ports names
year<- sort(unique(age_samples_ano_mes$ANO)) ##list of years on the samples data

for(bb in 1:length(year)){
  plotb<- ggplot(age_samples_ano_mes[age_samples_ano_mes$ANO==year[bb],], aes(x=IDADE,y=n, colour=PORTO_NOME))+xlab("age")+ ylab("number of individuals")+ geom_line()+ theme_classic() + facet_wrap(~factor(MONTH))+
   theme(axis.title.y = element_text(size = 14),axis.title.x=element_text(size=14),
         axis.line = element_line(size = 0.5),axis.text = element_text(size = 10))
   dev.copy(png, paste(output_dir, year[bb],"_age_distribution_samples_Port_year",".png",sep=""))
   print(plotb)
   dev.off()
}



## Figura 1 - Length distribution by year
years<- year[year %in% unique(data_samplebio$ANO)]
#lines_plot<-round(length(years)/2)
#par(mfrow=c(1,1))
for(nb in 1: length(years)){
  fig1<- hist(data_samplebio$C_CLASSE[data_samplebio$ANO==years[nb]],xlab="length", ylab="number of individuals", main=years[nb])
  dev.copy(png, paste(output_dir, years[nb],"_length_distribution",".png",sep=""))
  dev.off()
}


###Figure 2 - Age distribution by year

for(nb in 1: length(years)){
  fig2<- hist(data_samplebio$IDADE[data_samplebio$ANO==years[nb]],xlab="age", ylab="number of individuals", main=years[nb])
  dev.copy(png, paste(output_dir, years[nb],"_age_distribution",".png",sep=""))
  dev.off()
}

#########################################################################################################
#########################################################################################################
#########################################################################################################
