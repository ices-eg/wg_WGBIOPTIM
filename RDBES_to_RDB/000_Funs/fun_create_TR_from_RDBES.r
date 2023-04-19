#' Create a RDB TR table from RDBES data
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
#' @return a data.frame that approximates an RDB TR table
#' 
#' @export
#'
#' @examples
#' \dontrun{
#'
#' myH1RawObject <-
#'   createRDBESDataObject(rdbesExtractPath = "tests\\testthat\\h1_v_1_19_13")
#'   
#' create_TR_from_RDBES(x = myH1RawObject, hierarchy = 1)
#' }


create_TR_from_RDBES<-function(x, hierarchy){

# note: some approximations - see notes in individual lines

# x is RDBESobj
# hierarchy is hierarchy (integer)
# outputs an approximate TR table that, e.g., can be compared with FD2 outputs

if(!exists("createTableOfRDBESIds")) stop ("you need to load 'createTableOfRDBESIds' ahead of runing this function")
if(!hierarchy %in% c(1,5)) stop ("hierarchy not defined")

id_table <- createTableOfRDBESIds(x, hierarchy)

# adds the remaining info [hierarchy specific]
if(hierarchy == 1){
	TR_base<-data.table(FTid=x$FT$FTid)
	TR_base$Record_type<-"TR"
	TR_base$Sampling_type<-"S"
	TR_base$Landing_country <- substring(x$FT$FTarvLoc,1,2) # approximated from FTarvLoc (no LE table in hierarchy)
	TR_base$Vessel_flag_country<-x$VD$VDctry[match(x$FT$FTencrVessCode[match(TR_base$FTid,x$FT$FTid)],x$VD$VDencrVessCode)]
	TR_base$Year<-x$DE$DEyear[match(id_table$DEid[match(TR_base$FTid, id_table$FTid)],x$DE$DEid)]
	TR_base$Project<-x$DE$DEstratumName[match(id_table$DEid[match(TR_base$FTid, id_table$FTid)],x$DE$DEid)]
	TR_base$Trip_number<-do.call("rbind", strsplit(x$FT$FTunitName[match(TR_base$FTid,x$FT$FTid)]," "))[,1]
	TR_base$Vessel_length<-x$VD$VDlen[match(x$FT$FTencrVessCode[match(TR_base$FTid,x$FT$FTid)],x$VD$VDencrVessCode)]
	TR_base$Vessel_power <- x$VD$VDpwr[match(x$FT$FTencrVessCode[match(TR_base$FTid,x$FT$FTid)],x$VD$VDencrVessCode)]
	TR_base$Vessel_size <- x$VD$VDton[match(x$FT$FTencrVessCode[match(TR_base$FTid,x$FT$FTid)],x$VD$VDencrVessCode)]
	TR_base$Vessel_type<-"" # not available but not important
	TR_base$Harbour <- x$FT$FTarvLoc # approximate
	TR_base$No_SetsHauls_on_trip <- x$FT$FTfoNum
	TR_base$Days_at_sea <- as.integer(as.Date(x$FT$FTarvDate)- as.Date(x$FT$FTdepDate)+1)
	TR_base$Vessel_identifier <- x$FT$FTencrVessCode # approximate
	TR_base$Sampling_country <- "SWE"
	TR_base$Sampling_method <- "#"
}

if(hierarchy == 5){
	TR_base<-data.table(LEid=x$LE$LEid)
	TR_base$Record_type<-"TR"
	TR_base$Sampling_type<-"M"
	TR_base$Landing_country <- x$LE$LEctry
	TR_base$Vessel_flag_country<-x$VD$VDctry[match(x$LE$LEencrVessCode[match(TR_base$LEid,x$LE$LEid)],x$VD$VDencrVessCode)]
	TR_base$Year<-x$DE$DEyear[match(id_table$DEid[match(TR_base$LEid, id_table$LEid)],x$DE$DEid)]
	TR_base$Project<-x$DE$DEstratumName[match(id_table$DEid[match(TR_base$LEid, id_table$LEid)],x$DE$DEid)]
	TR_base$Trip_number<-do.call("rbind", strsplit(x$LE$LEunitName[match(TR_base$LEid,x$LE$LEid)]," "))[,1]
	TR_base$Vessel_length<-x$VD$VDlen[match(x$LE$LEencrVessCode[match(TR_base$LEid,x$LE$LEid)],x$VD$VDencrVessCode)]
	TR_base$Vessel_power <- x$VD$VDpwr[match(x$LE$LEencrVessCode[match(TR_base$LEid,x$LE$LEid)],x$VD$VDencrVessCode)]
	TR_base$Vessel_size <- x$VD$VDton[match(x$LE$LEencrVessCode[match(TR_base$LEid,x$LE$LEid)],x$VD$VDencrVessCode)]
	TR_base$Vessel_type<-"" # not available but not important
	TR_base$Harbour <- x$LE$LElocode # approximate
	TR_base$No_SetsHauls_on_trip <- "" # ATT available only if FT is included in H5 [FTnumberOfHaulsOrSets]
	TR_base$Days_at_sea <- "" # ATT available only if FT is included in H5 [FTarrivalDate-FTdepDate+1]
	TR_base$Vessel_identifier <- x$LE$LEencrVessCode # approximate
	TR_base$Sampling_country <- "SWE"
	TR_base$Sampling_method <- "#"

}

# FD2 specific format
TR_base$Landing_country<-ifelse(TR_base$Landing_country=="SE","SWE",TR_base$Landing_country)
TR_base$Vessel_flag_country<-ifelse(TR_base$Vessel_flag_country=="SE","SWE",TR_base$Vessel_flag_country)


	target_cols<-c('Record_type','Sampling_type','Landing_country','Vessel_flag_country','Year','Project','Trip_number','Vessel_length','Vessel_power','Vessel_size','Vessel_type',
	'Harbour','No_SetsHauls_on_trip','Days_at_sea','Vessel_identifier','Sampling_country','Sampling_method')

# head(unique(data.frame(CA_base)[,target_cols]))
# target_cols[!target_cols %in% 	colnames(TR_base)]

unique(data.frame(TR_base)[,target_cols])
}

# e.g.,
# create_TR_from_RDBES(x = RDBESprepObj, hierarchy = 1)