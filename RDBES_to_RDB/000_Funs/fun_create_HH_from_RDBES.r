#' Create a HH table from RDBES data
#'
#' ATT needs developments for different lower hierarchies  
#' (as of 20230418 only hierarchy = 1 and hierarchy = 5 specified)
#'
#' ATT requires pre-loading of create_RDBES_Id_table
#'
#' ATT some specifications are SLU specific (may need generalization, change to RDB colnames, etc.)
#'
#' @param x RDBESdataObject
#' @param hierarchy is RDBES hierarchy present in the RDBESdataObject (as integer)
#'
#' @return a data.frame that approximates an RDB HH table
#' 
#' @export
#'
#' @examples
#' \dontrun{
#'
#' myH1RawObject <-
#'   createRDBESDataObject(rdbesExtractPath = "tests\\testthat\\h1_v_1_19_13")
#'   
#' create_HH_from_RDBES(x = myH1RawObject, hierarchy = 1)
#' }



create_HH_from_RDBES<-function(x, hierarchy){

if(!exists("createTableOfRDBESIds")) stop ("you need to load createTableOfRDBESIds ahead of runing this function")
if(!hierarchy %in% c(1,5)) stop ("hierarchy not defined")

id_table <- createTableOfRDBESIds(x, hierarchy)

# adds the remaining info [hierarchy specific]
if(hierarchy == 1){
	HH_base<-data.table(FTid=x$FO$FTid, FOid=x$FO$FOid)
	HH_base$Record_type <- "HH"
	HH_base$Sampling_type <- "S"
	HH_base$Landing_country <- substring(x$FT$FTarvLoc[match(HH_base$FTid,x$FT$FTid)],1,2) # approximated from FTarvLoc (no LE table in hierarchy)
	HH_base$Vessel_flag_country <- x$VD$VDctry[match(x$FT$FTencrVessCode[match(HH_base$FTid,x$FT$FTid)],x$VD$VDencrVessCode)]
	HH_base$Year <- x$DE$DEyear[match(id_table$DEid[match(HH_base$FOid, id_table$FOid)],x$DE$DEid)]
	HH_base$Project <- x$DE$DEstratumName[match(id_table$DEid[match(HH_base$FOid, id_table$FOid)],x$DE$DEid)]
	HH_base$Trip_number <- do.call("rbind", strsplit(x$FT$FTunitName[match(HH_base$FTid,x$FT$FTid)]," "))[,1]
	HH_base$Station_number <- do.call("rbind", strsplit(x$FO$FOunitName[match(HH_base$FOid,x$FO$FOid)]," "))[,1]
	HH_base$Fishing_validity <- x$FO$FOval
	HH_base$Aggregation_level <- x$FO$FOaggLev
	HH_base$Catch_registration <- x$FO$FOcatReg
	HH_base$Species_registration <- "All" # sweden specific
	HH_base$Date <- x$FO$FOstartDate
	HH_base$Time <- x$FO$FOstartTime
	HH_base$Fishing_duration <- x$FO$FOdur
	HH_base$Pos_Start_Lat_dec <- x$FO$FOstartLat
	HH_base$Pos_Start_Lon_dec <- x$FO$FOstartLon
	HH_base$Pos_Stop_Lat_dec <- x$FO$FOstopLat
	HH_base$Pos_Stop_Lon_dec <- x$FO$FOstopLon
	HH_base$Area <- x$FO$FOarea
	HH_base$Statistical_rectangle <- x$FO$FOstatRect
	HH_base$Sub_polygon <- ""
	HH_base$Main_fishing_depth <- ""
	HH_base$Main_water_depth <- ""
	HH_base$FAC_National <- x$FO$FOnatFishAct
	HH_base$FAC_EC_lvl5 <- x$FO$FOmetier5
	HH_base$FAC_EC_lvl6 <- x$FO$FOmetier6
	HH_base$Gear_type <- x$FO$FOgear
	HH_base$Mesh_size <- x$FO$FOmeshSize
	HH_base$Selection_device <- x$FO$FOselDev
	HH_base$Mesh_size_selection_device <- x$FO$FOselDevMeshSize
}

if(hierarchy == 5){
	HH_base<-data.table(LEid=x$LE$LEid)
	HH_base$Record_type<-"TR"
	HH_base$Sampling_type<-"M"
	HH_base$Landing_country <- x$LE$LEctry[match(HH_base$LEid,x$LE$LEid)]
	HH_base$Vessel_flag_country<-x$VD$VDctry[match(x$LE$LEencrVessCode[match(HH_base$LEid,x$LE$LEid)],x$VD$VDencrVessCode)]
	HH_base$Year<-x$DE$DEyear[match(id_table$DEid[match(HH_base$LEid, id_table$LEid)],x$DE$DEid)]
	HH_base$Project<-x$DE$DEstratumName[match(id_table$DEid[match(HH_base$LEid, id_table$LEid)],x$DE$DEid)]
	HH_base$Trip_number<-do.call("rbind", strsplit(x$LE$LEunitName[match(HH_base$LEid,x$LE$LEid)]," "))[,1]
	HH_base$Station_number <- 999
	HH_base$Fishing_validity <- "V"
	HH_base$Aggregation_level <- "T"
	HH_base$Catch_registration <- x$LE$LEcatReg
	HH_base$Species_registration <- "Par" # sweden specific
	HH_base$Date <- x$LE$LEdate
	HH_base$Time <- x$LE$LEdate
	HH_base$Fishing_duration <- ""
	HH_base$Pos_Start_Lat_dec <- ""
	HH_base$Pos_Start_Lon_dec <- ""
	HH_base$Pos_Stop_Lat_dec <- ""
	HH_base$Pos_Stop_Lon_dec <- ""
	HH_base$Area <- x$LE$LEarea
	HH_base$Statistical_rectangle <- x$LE$LEstatRect
	HH_base$Sub_polygon <- ""
	HH_base$Main_fishing_depth <- ""
	HH_base$Main_water_depth <- ""
	HH_base$FAC_National <- x$LE$LEnatFishAct
	HH_base$FAC_EC_lvl5 <- x$LE$LEmetier5
	HH_base$FAC_EC_lvl6 <- x$LE$LEmetier6
	HH_base$Gear_type <- x$LE$LEgear
	HH_base$Mesh_size <- x$LE$LEmeshSize
	HH_base$Selection_device <- x$LE$LEselDev
	HH_base$Mesh_size_selection_device <- x$LE$LEselDevMeshSize

}

# FD2 specific format
HH_base$Landing_country<-ifelse(HH_base$Landing_country=="SE","SWE",HH_base$Landing_country)
HH_base$Vessel_flag_country<-ifelse(HH_base$Vessel_flag_country=="SE","SWE",HH_base$Vessel_flag_country)

	target_cols<-c('Record_type','Sampling_type','Landing_country','Vessel_flag_country','Year','Project','Trip_number','Station_number','Fishing_validity','Aggregation_level','Catch_registration','Species_registration','Date','Time','Fishing_duration','Pos_Start_Lat_dec','Pos_Start_Lon_dec','Pos_Stop_Lat_dec','Pos_Stop_Lon_dec','Area','Statistical_rectangle','Sub_polygon','Main_fishing_depth','Main_water_depth','FAC_National','FAC_EC_lvl5','FAC_EC_lvl6','Gear_type','Mesh_size','Selection_device','Mesh_size_selection_device')

# head(unique(data.frame(CA_base)[,target_cols]))
# target_cols[!target_cols %in% 	colnames(HH_base)]

unique(data.frame(HH_base)[,target_cols])
}
# e.g.,
# create_HH_from_RDBES(x = RDBESprepObj, hierarchy = 1)