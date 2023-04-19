#' Create a CA table from RDBES data
#'
#' ATT needs developments for different hierarchies 
#' (as of 20230418 only hierarchy = 1 and hierarchy = 5 specified)
#'
#' ATT requires pre-loading of createTableOfRDBESIds
#'
#' ATT some specifications are SLU specific (may need generalization, change to RDB colnames, etc.)
#'
#' @param x RDBESdataObject
#' @param hierarchy is RDBES hierarchy present in the RDBESdataObject (as integer)
#'
#' @return a data.frame that approximates an RDB CA table
#' 
#' @export
#'
#' @examples
#' \dontrun{
#'
#' myH1RawObject <-
#'   createRDBESDataObject(rdbesExtractPath = "tests\\testthat\\h1_v_1_19_13")
#'   
#' create_CA_from_RDBES(x = myH1RawObject, hierarchy = 1)
#' }


create_CA_from_RDBES<-function(x, hierarchy){

#if(!exists("createTableOfIds")) stop ("you need to load fun_createTableOfIds ahead of runing this function")
if(!exists("createTableOfRDBESIds")) stop ("you need to load 'createTableOfRDBESIds' ahead of runing this function")
if(!hierarchy %in% c(1,5)) stop ("hierarchy not defined")

id_table <- createTableOfRDBESIds(x, hierarchy)

# retrieves the basic BV info [common to all hierarchies]
CA_base <- x$BV[,list(Length_class=BVvalueMeas[BVtypeMeas=="LengthTotal"], Length_code = BVaccuracy[BVtypeMeas=="LengthTotal"], Sex=BVvalueMeas[BVtypeMeas=="Sex"],
			Individual_weight = as.numeric(BVvalueMeas[BVtypeMeas=="WeightMeasured"])*	BVconFacAssess[BVtypeMeas=="WeightMeasured"],
			Age = BVvalueMeas[BVtypeMeas=="Age"], AgeQuality = BVcertaintyQuali[BVtypeMeas=="Age"], Specimen_No=BVfishId),.(BVfishId)]	

# adds the remaining info [hierarchy specific]
if(hierarchy == 1){
	CA_base$Record_type<-"CA"
	CA_base$Sampling_type<-"S"
	CA_base$Landing_country <- "SWE" # could be made general
	CA_base$Vessel_flag_country<-"SWE" # could be made general
	CA_base$Year<-x$DE$DEyear[match(id_table$DEid[match(CA_base$BVfishId, id_table$BVfishId)],x$DE$DEid)]
	CA_base$Project<-x$DE$DEstratumName[match(id_table$DEid[match(CA_base$BVfishId, id_table$BVfishId)],x$DE$DEid)]
	CA_base$Trip_number<-do.call("rbind", strsplit(x$FT$FTunitName[match(id_table$FTid[match(CA_base$BVfishId, id_table$BVfishId)],x$FT$FTid)]," "))[,1] 
	CA_base$Station_number<-do.call("rbind", strsplit(x$FO$FOunitName[match(id_table$FOid[match(CA_base$BVfishId, id_table$BVfishId)],x$FO$FOid)]," "))[,1]
	# check with FD2 - FTdepDate? or FTarvDate
	CA_base$Quarter<-quarter(x$FT$FTdepDate[match(id_table$FTid[match(CA_base$BVfishId, id_table$BVfishId)],x$FT$FTid)]) 
	CA_base$Month<-month(x$FT$FTdepDate[match(id_table$FTid[match(CA_base$BVfishId, id_table$BVfishId)],x$FT$FTid)]) 
	
	CA_base$Species<-x$SA$SAspeCode[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)] 
	CA_base$Catch_category<-toupper(x$SA$SAcatchCat[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)]) 
	CA_base$Landing_category<-toupper(x$SA$SAlandCat[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)]) 
	CA_base$Comm_size_cat_scale<-x$SA$SAcommCatScl[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)] 
	CA_base$Comm_size_cat<-x$SA$SAcommCat[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)]  
	CA_base$Stock<-"" 
	CA_base$Area<-x$FO$FOarea[match(id_table$FOid[match(CA_base$BVfishId, id_table$BVfishId)],x$FO$FOid)]
	CA_base$Statistical_rectangle<-x$FO$FOstatRect[match(id_table$FOid[match(CA_base$BVfishId, id_table$BVfishId)],x$FO$FOid)] 
	CA_base$Subpolygon<-"" 
	CA_base$Aging_method<-"" 
	CA_base$Age_plus_group<-"-" 
	CA_base$Otolith_weight<-""
	CA_base$Otolith_side<-""
	CA_base$Maturity_staging<-""
	CA_base$Maturity_scale<-""
	CA_base$Maturity_stage<-""
}

if(hierarchy == 5){
	CA_base$Record_type<-"CA"
	CA_base$Sampling_type<-"M"
	CA_base$Landing_country <- x$LE$LEctry[match(id_table$LEid[match(CA_base$BVfishId, id_table$BVfishId)],x$LE$LEid)]
	aux<-x$LE$LEencrVessCode[match(id_table$LEid[match(CA_base$BVfishId, id_table$BVfishId)],x$LE$LEid)]
	CA_base$Vessel_flag_country <- x$VD$VDctry[match(aux,x$VD$VDencrVessCode)]
	#CA_base$Vessel_flag_country<-"SWE"
	CA_base$Year<-x$DE$DEyear[match(id_table$DEid[match(CA_base$BVfishId, id_table$BVfishId)],x$DE$DEid)]
	CA_base$Project<-x$DE$DEstratumName[match(id_table$DEid[match(CA_base$BVfishId, id_table$BVfishId)],x$DE$DEid)]
	CA_base$Trip_number<-do.call("rbind", strsplit(x$LE$LEunitName[match(id_table$LEid[match(CA_base$BVfishId, id_table$BVfishId)],x$LE$LEid)]," "))[,1] # TRIP_ID
	CA_base$Station_number<-999 
	CA_base$Quarter<-quarter(x$LE$LEdate[match(id_table$LEid[match(CA_base$BVfishId, id_table$BVfishId)],x$LE$LEid)]) 
	CA_base$Month<-month(x$LE$LEdate[match(id_table$LEid[match(CA_base$BVfishId, id_table$BVfishId)],x$LE$LEid)]) 
	CA_base$Species<-x$SA$SAspeCode[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)] 
	CA_base$Catch_category<-toupper(x$SA$SAcatchCat[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)]) 
	CA_base$Landing_category<-toupper(x$SA$SAlandCat[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)]) 
	CA_base$Comm_size_cat_scale<-x$SA$SAcommCatScl[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)] 
	CA_base$Comm_size_cat<-x$SA$SAcommCat[match(id_table$SAid[match(CA_base$BVfishId, id_table$BVfishId)],x$SA$SAid)]  
	CA_base$Stock<-"" 
	CA_base$Area<-x$LE$LEarea[match(id_table$LEid[match(CA_base$BVfishId, id_table$BVfishId)],x$LE$LEid)]
	CA_base$Statistical_rectangle<-x$LE$LEstatRect[match(id_table$LEid[match(CA_base$BVfishId, id_table$BVfishId)],x$LE$LEid)] 
	CA_base$Subpolygon<-"" 
	CA_base$Aging_method<-"" 
	CA_base$Age_plus_group<-"-" 
	CA_base$Otolith_weight<-""
	CA_base$Otolith_side<-""
	CA_base$Maturity_staging<-""
	CA_base$Maturity_scale<-""
	CA_base$Maturity_stage<-""
}

	target_cols<-c("Record_type","Sampling_type","Landing_country","Vessel_flag_country","Year","Project","Trip_number","Station_number",
			"Quarter","Month","Species","Sex","Catch_category","Landing_category","Comm_size_cat_scale","Comm_size_cat","Stock","Area",
			"Statistical_rectangle","Subpolygon","Length_class","Age","Specimen_No","Length_code","Aging_method","Age_plus_group","Otolith_weight","Otolith_side","Individual_weight","Maturity_staging","Maturity_scale","Maturity_stage")

# head(unique(data.frame(CA_base)[,target_cols]))
# target_cols[!target_cols %in% 	colnames(CA_base)]

unique(data.frame(CA_base)[,target_cols])
}
# e.g.,
# create_CA_from_RDBES(x = RDBESprepObj, hierarchy = 1)