#' Create a HL table from RDBES data
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
#' @return a data.frame that approximates an RDB HL table
#' 
#' @export
#'
#' @examples
#' \dontrun{
#'
#' myH1RawObject <-
#'   createRDBESDataObject(rdbesExtractPath = "tests\\testthat\\h1_v_1_19_13")
#'   
#' create_HL_from_RDBES(x = myH1RawObject, hierarchy = 1)
#' }

create_HL_from_RDBES<-function(x, hierarchy, dir_codelists){


if(!exists("codeLists")){ print("read RDBES codelists")
		filename<-"RDBES_codeLists.Rdata"
		pathfile<-paste0(dir_codelists, filename,sep="")
		load(pathfile)} 
if(!exists("createTableOfRDBESIds")) stop ("you need to load 'createTableOfRDBESIds' ahead of runing this function")
if(!hierarchy %in% c(1,5)) stop ("hierarchy not defined")

id_table <- createTableOfRDBESIds(x, hierarchy)


# adds the remaining info [hierarchy specific]
if(hierarchy == 1){
	HL_base <- x$FM[FMtypeMeas %in% c("LengthTotal","LengthCarapace"),list(FMid,SAid,FMclassMeas,FMnumAtUnit,FMconFacAssess)]
	HL_base$FTid<-id_table$FTid[match(HL_base$SAid, id_table$SAid)]
	HL_base$FOid<-id_table$FOid[match(HL_base$SAid, id_table$SAid)]

	HL_base$Record_type <- "HL"
	HL_base$Sampling_type <- "S"
	HL_base$Landing_country <- "SWE" # approximate (no LE table in hierarchy)
	HL_base$Vessel_flag_country <- x$VD$VDctry[match(x$FT$FTencrVessCode[match(HL_base$FTid,x$FT$FTid)],x$VD$VDencrVessCode)]
	HL_base$Year <- x$DE$DEyear[match(id_table$DEid[match(HL_base$FOid, id_table$FOid)],x$DE$DEid)]
	HL_base$Project <- x$DE$DEstratumName[match(id_table$DEid[match(HL_base$FOid, id_table$FOid)],x$DE$DEid)]
	HL_base$Trip_number <- do.call("rbind", strsplit(x$FT$FTunitName[match(HL_base$FTid,x$FT$FTid)]," "))[,1]
	HL_base$Station_number <- do.call("rbind", strsplit(x$FO$FOunitName[match(HL_base$FOid,x$FO$FOid)]," "))[,1]
	HL_base$Species <- codeLists[["SpecWoRMS"]]$Description[match(x$SA$SAspeCode[match(HL_base$SAid,x$SA$SAid)],codeLists[["SpecWoRMS"]]$Key)]
	HL_base$Catch_category <- toupper(x$SA$SAcatchCat[match(HL_base$SAid,x$SA$SAid)])
	HL_base$Landing_category <- toupper(x$SA$SAlandCat[match(HL_base$SAid,x$SA$SAid)])
	HL_base$Comm_size_cat_scale <- x$SA$SAcommCatScl[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Comm_size_cat <- x$SA$SAcommCat[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Subsampling_category <- x$SA$SAstratumName[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Sex <- x$SA$SAsex[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Individual_sex <- x$SA$SAsex[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Length_class <- HL_base$FMclassMeas*HL_base$FMconFacAssess
	HL_base$Number_at_length <- HL_base$FMnumAtUnit

}

 if(hierarchy == 5){
 	HL_base <- x$BV[BVtypeMeas %in% c("LengthTotal","LengthCarapace"),list(FMnumAtUnit=.N),list(SAid,FMclassMeas=as.numeric(BVvalueMeas),BVconFacAssess)]
	HL_base$LEid<-id_table$LEid[match(HL_base$SAid, id_table$SAid)]
 
	HL_base$Record_type <- "HL"
	HL_base$Sampling_type <- "M"
	HL_base$Landing_country <- x$LE$LEctry[match(HL_base$LEid,x$LE$LEid)]
	HL_base$Vessel_flag_country <- x$VD$VDctry[match(x$LE$LEencrVessCode[match(HL_base$LEid,x$LE$LEid)],x$VD$VDencrVessCode)]
	HL_base$Year <- x$DE$DEyear[match(id_table$DEid[match(HL_base$LEid, id_table$LEid)],x$DE$DEid)]
	HL_base$Project <- x$DE$DEstratumName[match(id_table$DEid[match(HL_base$LEid, id_table$LEid)],x$DE$DEid)]
	HL_base$Trip_number <- do.call("rbind", strsplit(x$LE$LEunitName[match(HL_base$LEid,x$LE$LEid)]," "))[,1]
	HL_base$Station_number <- 999
	HL_base$Species <- codeLists[["SpecWoRMS"]]$Description[match(x$SA$SAspeCode[match(HL_base$SAid,x$SA$SAid)],codeLists[["SpecWoRMS"]]$Key)]
	HL_base$Catch_category <- toupper(x$SA$SAcatchCat[match(HL_base$SAid,x$SA$SAid)])
	HL_base$Landing_category <- toupper(x$SA$SAlandCat[match(HL_base$SAid,x$SA$SAid)])
	HL_base$Comm_size_cat_scale <- x$SA$SAcommCatScl[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Comm_size_cat <- x$SA$SAcommCat[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Subsampling_category <- x$SA$SAstratumName[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Sex <- x$SA$SAsex[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Individual_sex <- x$SA$SAsex[match(HL_base$SAid,x$SA$SAid)]
	HL_base$Length_class <- HL_base$FMclassMeas*HL_base$BVconFacAssess
	HL_base$Number_at_length <- HL_base$FMnumAtUnit
 }

# FD2 specific format
HL_base$Landing_country<-ifelse(HL_base$Landing_country=="SE","SWE",HL_base$Landing_country)
HL_base$Vessel_flag_country<-ifelse(HL_base$Vessel_flag_country=="SE","SWE",HL_base$Vessel_flag_country)


	target_cols<-c('Record_type','Sampling_type','Landing_country','Vessel_flag_country','Year','Project','Trip_number','Station_number','Species','Catch_category','Landing_category','Comm_size_cat_scale','Comm_size_cat','Subsampling_category','Sex','Individual_sex','Length_class','Number_at_length')

# head(unique(data.frame(CA_base)[,target_cols]))
# target_cols[!target_cols %in% 	colnames(HL_base)]

unique(data.frame(HL_base)[,target_cols])
}
# e.g.,
# create_HL_from_RDBES(x = RDBESprepObj, hierarchy = 1, dir_codelists = "C:/WORK/Design&Analysis/RDBESaux/002_Outputs/v1.19.12/")