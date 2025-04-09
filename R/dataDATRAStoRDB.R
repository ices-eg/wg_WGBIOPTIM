#' dataRDBEStoRDB
#'
#' @param isurvey contains survey name (E.g. "BTS","BITS","NS-IBTS")  to get a list of available surveys: icesDatras::getSurveyList()
#' @param iyears contains years of interest (E.g. 2020 or 2020:2021)
#' @param iquarters contains quarters of interest (E.g. 3 or 3:4)
#' @param icountry contains country of interest, using ISO 3166  (E.g. GB-ENG, BE)
#'
#' @return list(HL = HL, HH = HH, SL = SL, CA = CA) RDB tables
#' @export
#'
#' @examples
#'dataDATRAStoRDB("BTS",2021,3,"BE")
#'

dataDATRAStoRDB <- function(isurvey,iyears,iquarters,icountry)

{
  
  
# Load packages
require(icesDatras)
require(dplyr)
require(vmstools)
if (!require(vmstools)) {
    # Handle the case where ggplot2 is not installed
  library(devtools)
  install_github("nielshintzen/vmstools/vmstools/")
  library(vmstools)
  }
  
  
  
# Retrieve the data from DATRAS
print("downloading data from https://datras.ices.dk/WebServices/, this might take awhile")
hh_data <- getDATRAS(record = "HH", survey = isurvey, years = iyears, quarters = iquarters)%>%
filter(Country == icountry)
print("Downloaded hh data from DATRAS")
hl_data <- getDATRAS(record = "HL", survey = isurvey, years = iyears, quarters = iquarters)%>%
filter(Country == icountry)
print("Downloaded hl data from DATRAS")
ca_data <- getDATRAS(record = "CA", survey = isurvey, years = iyears, quarters = iquarters)%>%
filter(Country == icountry)
print("Downloaded ca data from DATRAS")


########################## HH ###################
HH <- data.frame(Record_type = rep("HH", nrow(hh_data)),  stringsAsFactors = FALSE)
HH$Record_type <- "HH"
HH$Sampling_type <- "S"
HH$Landing_country <- hh_data$Country
HH$Vessel_flag_country <- hh_data$Country
HH$Year <- hh_data$Year
HH$Project <- hh_data$Survey
HH$Trip_number <-  paste(hh_data$Country,hh_data$Survey,hh_data$Year,sep="_")
HH$Station_number <- hh_data$StNo
HH$Fishing_validity <- hh_data$HaulVal
HH$Aggregation_level <- "H"
HH$Catch_registration <- "All"
HH$Species_registration <- "All"
HH$Date <- as.Date(ISOdate(hh_data$Year, hh_data$Month, hh_data$Day))
HH$Time <- sprintf("%02d:%02d", as.integer(substr(hh_data$TimeShot , 1, 2)), as.integer(substr(hh_data$TimeShot, 3, 4)))
HH$Fishing_duration <- hh_data$HaulDur
HH$Pos_Start_Lat_dec <- hh_data$ShootLat
HH$Pos_Start_Lon_dec <- hh_data$ShootLong 
HH$Pos_Stop_Lat_dec <- hh_data$HaulLat
HH$Pos_Stop_Lon_dec <- hh_data$HaulLong
# get ICES areas based on haul coordinates, using R package VMStools
hh_data$SI_LONG<-hh_data$HaulLong
hh_data$SI_LATI<-hh_data$HaulLat
data(ICESareas)
areas<-ICESarea(hh_data,ICESareas,fast=T)
hh_data$Area<-ICESareas$Area_Full[areas]
HH$Area <- hh_data$Area
HH$Statistical_rectangle <- hh_data$StatRec
HH$Sub_polygon <- ''
HH$Main_fishing_depth <- hh_data$MinTrawlDepth
HH$Main_water_depth <- hh_data$DepthStratum
HH$FAC_National <- hh_data$Gear
HH$FAC_EC_lvl5 <- hh_data$Gear
HH$FAC_EC_lvl6 <- hh_data$Gear
HH$Gear_type <- hh_data$Gear
HH$Mesh_size <- hh_data$CodendMesh
HH$Selection_device <- 0
HH$Mesh_size_selection_device <- ''
HH[is.na(HH)] <- ''
HH <- distinct(HH)


########################## HL ###################

HL <- data.frame(Record_type = rep("HL", nrow(hl_data)),  stringsAsFactors = FALSE)
HL$Sampling_type <- "S"
HL$Landing_country <- hl_data$Country
HL$Vessel_flag_country <- hl_data$Country
HL$Year <- hl_data$Year
HL$Project <- hl_data$Survey
HL$Trip_number <- paste(hl_data$Country,hl_data$Survey,hl_data$Year,sep="_")
HL$Station_number <- hl_data$StNo
HL$Species <- hl_data$Valid_Aphia
HL$Catch_category <- ''
HL$Landing_category <-''
HL$Comm_size_cat_scale <- ''
HL$Comm_size_cat <- ''
HL$Subsampling_category <- ''								
HL$Sex <- hl_data$Sex
HL$Individual_sex <- hl_data$Sex
HL$Length_class <- hl_data$LngtClass
HL$Number_at_length <- hl_data$HLNoAtLngt 
HL[is.na(HL)] <- ''
HL <- aggregate(Number_at_length ~ ., data=HL, function(x) sum(x,na.rm=TRUE))



########################## SL ###################

SL <- data.frame(Record_type = rep("SL", nrow(hl_data)),  stringsAsFactors = FALSE)
SL$Sampling_type <- "S"
SL$Landing_country <- hl_data$Country
SL$Vessel_flag_country <- hl_data$Country
SL$Year <- hl_data$Year
SL$Project <- hl_data$Survey
SL$Trip_number <- paste(hl_data$Country,hl_data$Survey,hl_data$Year,sep="_")
SL$Station_number <- hl_data$StNo
SL$Species <- hl_data$Valid_Aphia
SL$Catch_category <- ''
SL$Landing_category <- ''
SL$Comm_size_cat_scale <- ''
SL$Comm_size_cat <- ''
SL$Subsampling_category <- ''								
SL$Sex <- hl_data$Sex
SL$Weight <- hl_data$CatCatchWgt
SL$Subsample_weight <- hl_data$SubWgt
SL$Length_code <- "unknown"
SL[is.na(SL)] <- ''
SL <- distinct(SL)


########################## CA ###################

CA <- data.frame(Record_type = rep("CA", nrow(ca_data)),  stringsAsFactors = FALSE)
CA$Record_type <- "CA"
CA$Sampling_type <- "S"
CA$Landing_country <- ca_data$Country
CA$Vessel_flag_country <- ca_data$Country
CA$Year <- ca_data$Year
CA$Project <- ca_data$Survey
CA$Trip_number <- ca_data$FTunitName
CA$Trip_number <- paste(ca_data$Country,ca_data$Survey,ca_data$Year,sep="_")
CA$Station_number <- ca_data$StNo
CA$Quarter <- ca_data$Quarter
hh_data$Trip_number_StNo <-  paste(hh_data$Country,hh_data$Survey,hh_data$Year,hh_data$StNo,sep="_")
ca_data$Trip_number_StNo <-  paste(ca_data$Country,ca_data$Survey,ca_data$Year,ca_data$StNo,sep="_")
ca_data$Month <- hh_data$Month[match(ca_data$Trip_number_StNo,hh_data$Trip_number_StNo)]
CA$Month<-ca_data$Month
CA$Species <- ca_data$Valid_Aphia
CA$Sex <- ca_data$Sex
CA$Catch_category <- ''
CA$Landing_category <-'' 
CA$Comm_size_cat_scale <-'' 
CA$Comm_size_cat <- ''
CA$Stock <- ''
ca_data$Area<-hh_data$Area[match(ca_data$Trip_number_StNo,hh_data$Trip_number_StNo)]
CA$Area <- ca_data$Area
ca_data$Rect<-hh_data$StatRec[match(ca_data$Trip_number_StNo,hh_data$Trip_number_StNo)]
CA$Statistical_rectangle <- ca_data$Rect  # rectangle info should also be available in the Ca, but here we use the data stored in the HH to ensure consistency
CA$Sub_polygon <- ''
CA$Length_class <- ca_data$LngtClass
CA$Age <- ca_data$Age
CA$Single_fish_number <- ca_data$FishID
CA$Length_code <- "unknown"
CA$Aging_method <- ca_data$AgePrepMet
CA$Age_plus_group <- ca_data$PlusGr 
CA$Otolith_weight <- ''
CA$Otolith_side <- ''
CA$Weight <- ca_data$IndWgt
CA$Maturity_staging_method <- ''
CA$Maturity_scale <- ca_data$MaturityScale
CA$Maturity_stage <- ca_data$Maturity
CA[is.na(CA)] <- ''


return(list(HL = HL, HH = HH, SL = SL, CA = CA))

}
