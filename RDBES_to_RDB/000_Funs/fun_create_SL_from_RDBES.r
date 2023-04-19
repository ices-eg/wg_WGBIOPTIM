#' Create a RDB SL table from RDBES data
#'
#' ATT needs developments for different other upper hierarchies  [only 1 and 5 for now] 
#' ATT needs developments for different lower hierarchies  [only C for now] 
#'
#' ATT requires pre-loading of create_RDBES_Id_table
#'
#' ATT some specifications are SLU specific (may need generalization, change to RDB colnames, etc.)
#'
#' @param x RDBESdataObject
#' @param hierarchy is RDBES hierarchy present in the RDBESdataObject (as integer)
#'
#' @return a data.frame that approximates an RDB SL table
#' 
#' @export
#'
#' @examples
#' \dontrun{
#'
#' myH1RawObject <-
#'   createRDBESDataObject(rdbesExtractPath = "tests\\testthat\\h1_v_1_19_13")
#'   
#' create_SL_from_RDBES(x = myH1RawObject, hierarchy = 1)
#' }


create_SL_from_RDBES<-function(x, hierarchy, dir_codelists){

if(!exists("codeLists")){ print("read RDBES codelists")
		filename<-"RDBES_codeLists.Rdata"
		pathfile<-paste0(dir_codelists, filename,sep="")
		load(pathfile)} 
if(!exists("createTableOfRDBESIds")) stop ("you need to load 'createTableOfRDBESIds' ahead of runing this function")
if(!hierarchy %in% c(1,5)) stop ("hierarchy not defined")

id_table <- createTableOfRDBESIds(x, hierarchy)


# removes the parent rows where they exist
table(x$SA$SSid,!is.na(x$SA$SAparSequNum))
t1<-table(paste(x$SA$SSid,x$SA$SAspeCode, x$SA$SAcatchCat),is.na(x$SA$SAparSequNum))
tmp<-rownames(t1)[apply(t1>0,1,sum)==2]
out<-x$SA[!(paste(x$SA$SSid,x$SA$SAspeCode, x$SA$SAcatchCat, x$SA$SAparSequNum) %in% paste(tmp,"NA")),]


# removes out-of-frame
out<-out[SAunitName>0,]

SL_base<-out[,list(SAid,SAspeCode,SAcatchCat,SAlandCat,SAcommCatScl,SAcommCat,SAstratumName,SAsex, SAtotalWtLive,SAsampWtLive,SAlowHierarchy)]

	if(any(SL_base$SAlowHierarchy %in% c("A","B"))) SL_base$FMaccuracy[SL_base$SAlowHierarchy %in% c("A","B")]<-x$FM$FMaccuracy[match(SL_base$SAid[SL_base$SAlowHierarchy %in% c("A","B")],x$FM$SAid)]
	if(any(SL_base$SAlowHierarchy =="C")) SL_base$FMaccuracy[SL_base$SAlowHierarchy=="C"]<-x$BV$BVaccuracy[x$BV$BVtypeMeas %in% c("LengthTotal","LengthCarapace")][match(SL_base$SAid[SL_base$SAlowHierarchy=="C"],x$BV$SAid[x$BV$BVtypeMeas %in% c("LengthTotal","LengthCarapace")])]
	SL_base$FMaccuracy[SL_base$SAlowHierarchy=="D"]<-""

# adds the remaining info [hierarchy specific]
if(hierarchy == 1){
	SL_base$FTid<-id_table$FTid[match(SL_base$SAid, id_table$SAid)]
	SL_base$FOid<-id_table$FOid[match(SL_base$SAid, id_table$SAid)]	
	
	SL_base$Record_type <- "SL"
	SL_base$Sampling_type <- "S"
	SL_base$Landing_country <- x$FT$FTarvLoc[match(SL_base$FTid,x$FT$FTid)]
	SL_base$Vessel_flag_country <- x$VD$VDctry[match(x$FT$FTencrVessCode[match(SL_base$FTid,x$FT$FTid)],x$VD$VDencrVessCode)]
	SL_base$Year <- x$DE$DEyear[match(id_table$DEid[match(SL_base$FOid, id_table$FOid)],x$DE$DEid)]
	SL_base$Project <- x$DE$DEstratumName[match(id_table$DEid[match(SL_base$FOid, id_table$FOid)],x$DE$DEid)]
	SL_base$Trip_number <- do.call("rbind", strsplit(x$FT$FTunitName[match(SL_base$FTid,x$FT$FTid)]," "))[,1]
	SL_base$Station_number <- do.call("rbind", strsplit(x$FO$FOunitName[match(SL_base$FOid,x$FO$FOid)]," "))[,1]
	SL_base$Species <- codeLists[["SpecWoRMS"]]$Description[match(SL_base$SAspeCode,codeLists[["SpecWoRMS"]]$Key)]
	SL_base$Catch_category <- toupper(SL_base$SAcatchCat)
	SL_base$Landing_category <- toupper(SL_base$SAlandCat)
	SL_base$Comm_size_cat_scale <- SL_base$SAcommCatScl
	SL_base$Comm_size_cat <- SL_base$SAcommCat
	SL_base$Subsampling_category <- SL_base$SAstratumName
	SL_base$Sex <- SL_base$SAsex
	SL_base$Weight <- SL_base$SAtotalWtLive
	SL_base$Subsample_weight <- SL_base$SAsampWtLive
	SL_base$Length_code <- SL_base$FMaccuracy
	
}

if(hierarchy == 5){
	SL_base$LEid<-id_table$LEid[match(SL_base$SAid, id_table$SAid)]

	SL_base$Record_type <- "SL"
	SL_base$Sampling_type <- "M"
	SL_base$Landing_country <- x$LE$LEctry[match(SL_base$LEid,x$LE$LEid)]
	SL_base$Vessel_flag_country <- x$VD$VDctry[match(x$LE$LEencrVessCode[match(SL_base$LEid,x$LE$LEid)],x$VD$VDencrVessCode)]
	SL_base$Year <- x$DE$DEyear[match(id_table$DEid[match(SL_base$LEid, id_table$LEid)],x$DE$DEid)]
	SL_base$Project <- x$DE$DEstratumName[match(id_table$DEid[match(SL_base$LEid, id_table$LEid)],x$DE$DEid)]
	SL_base$Trip_number <- do.call("rbind", strsplit(x$LE$LEunitName[match(SL_base$LEid,x$LE$LEid)]," "))[,1]
	SL_base$Station_number <- 999
	SL_base$Species <- codeLists[["SpecWoRMS"]]$Description[match(SL_base$SAspeCode,codeLists[["SpecWoRMS"]]$Key)]
	SL_base$Catch_category <- toupper(SL_base$SAcatchCat)
	SL_base$Landing_category <- toupper(SL_base$SAlandCat)
	SL_base$Comm_size_cat_scale <- SL_base$SAcommCatScl
	SL_base$Comm_size_cat <- SL_base$SAcommCat
	SL_base$Subsampling_category <- SL_base$SAstratumName
	SL_base$Sex <- SL_base$SAsex
	SL_base$Weight <- SL_base$SAtotalWtLive
	SL_base$Subsample_weight <- SL_base$SAsampWtLive
	SL_base$Length_code <- SL_base$FMaccuracy

}

# FD2 specific format
SL_base$Landing_country<-ifelse(SL_base$Landing_country=="SE","SWE",SL_base$Landing_country)
SL_base$Vessel_flag_country<-ifelse(SL_base$Vessel_flag_country=="SE","SWE",SL_base$Vessel_flag_country)


	target_cols<-c('Record_type','Sampling_type','Landing_country','Vessel_flag_country','Year','Project','Trip_number','Station_number','Species','Catch_category','Landing_category','Comm_size_cat_scale','Comm_size_cat','Subsampling_category','Sex','Weight','Subsample_weight','Length_code')

# head(unique(data.frame(CA_base)[,target_cols]))
# target_cols[!target_cols %in% 	colnames(SL_base)]

unique(data.frame(SL_base)[,target_cols])
}
# e.g.,
# create_SL_from_RDBES(x = RDBESprepObj, hierarchy = 1, dir_codelists = "C:/WORK/Design&Analysis/RDBESaux/002_Outputs/v1.19.12/")