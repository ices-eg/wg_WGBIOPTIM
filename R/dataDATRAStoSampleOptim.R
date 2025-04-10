isurvey<-"BTS"
iyears<-2021
iquarters<-1:4
icountry<-"BE"

dataDATRAStoSampleOptim("NS-IBTS",2021,1:4,"NL")

#' dataDATRAStoSampleOptim
#'
#' @param isurvey contains survey name (E.g. "BTS","BITS","NS-IBTS")  to get a list of available surveys: icesDatras::getSurveyList()
#' @param iyears contains years of interest (E.g. 2020 or 2020:2021)
#' @param iquarters contains quarters of interest (E.g. 3 or 3:4)
#' @param icountry contains country of interest, using ISO 3166  (E.g. GB-ENG, BE)
#'
#' @return list(data_samplebio = data_samplebio) SampleOptim input table
#' @export
#'
#' @examples
#'dataDATRAStoSampleOptim("BTS",2021,1:4,"BE")
#'

dataDATRAStoSampleOptim <- function(isurvey,iyears,iquarters,icountry)

{
  
  
# Load packages
require(dplyr)

require(icesDatras)
    if (!require(icesDatras)) {
    # Handle the case where icesDatras is not installed
    library(devtools)
    install.packages("icesDatras")
    library(icesDatras)
    }

require(icesVocab)
    if (!require(icesVocab)) {
    # Handle the case where icesVocab is not installed
    library(devtools)
    install.packages("icesVocab")
    library(icesVocab)
    }
  
  
# Retrieve the data from DATRAS
print("downloading data from https://datras.ices.dk/WebServices/, this might take awhile")
ca_data <- getDATRAS(record = "CA", survey = isurvey, years = iyears, quarters = iquarters)%>%
filter(Country == icountry)
if (nrow(ca_data) == 0) {
  stop("No CA data, probably incorrect country code or no data for a given survey x year x quarter x country combination")
}
print("Downloaded ca data from DATRAS")
hh_data <- getDATRAS(record = "HH", survey = isurvey, years = iyears, quarters = iquarters)%>%
  filter(Country == icountry)
print("Downloaded hh data from DATRAS")





#  replace AphiaID with FAO codes
FAOcodes <- getCodeList("SpecASFIS")
FAOcodes <- FAOcodes %>%
  rename_at("Description",~"Scientific") %>%
  rename_at("Key",~"FAO")
FAOcodes<-FAOcodes[,c(2,3)]

Aphiacodes <- getCodeList("SpecWoRMS")
Aphiacodes <- Aphiacodes %>%
  mutate("SpecCode" = as.integer(Key))%>%
  rename_at("Description",~"Scientific")
Aphiacodes<-Aphiacodes[,c(3,7)]
Aphiacodes <-Aphiacodes %>% left_join(FAOcodes,by= "Scientific",relationship = "many-to-many")

ca_data<-ca_data  %>% left_join(Aphiacodes,by= "SpecCode",relationship = "many-to-many")

# Get CA data into format for SampleOptim 
data_samplebio <- data.frame(id_bio_fish = as.integer(c(1:nrow(ca_data))),  stringsAsFactors = FALSE)
data_samplebio$date<-as.factor(ca_data$DateofCalculation)
# get month from HH, copy to CA
hh_data$Trip_number_StNo <-  paste(hh_data$Country,hh_data$Survey,hh_data$Year,hh_data$StNo,sep="_")
ca_data$Trip_number_StNo <-  paste(ca_data$Country,ca_data$Survey,ca_data$Year,ca_data$StNo,sep="_")
ca_data$month <- as.integer(hh_data$Month[match(ca_data$Trip_number_StNo,hh_data$Trip_number_StNo)])
data_samplebio$month<-ca_data$month
data_samplebio$year<-as.integer(ca_data$Year)
data_samplebio$fao_code<-ca_data$FAO
data_samplebio$port_name<-''
data_samplebio$length_class<-as.integer(ca_data$LngtClass)  
data_samplebio$individual_wg<-as.numeric(ca_data$IndWgt) 
ca_data <- ca_data %>%
  mutate(Sex = if_else(is.na(Sex) | Sex == "U", "I", Sex))
data_samplebio$sex<-as.factor(ca_data$Sex)
data_samplebio$maturity_stage<-ca_data$Maturity
data_samplebio$age<-as.integer(ca_data$Age) 

return(list(data_samplebio=data_samplebio))

}
