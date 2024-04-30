##### Input: RDBES tables and Hierarchy
##### Output: RDB tables - CA, HL, HH and SL
##### At present stage - only for H1 (and H2?), but may be extended


doRDBEStoRDB <- 
function(DEtable=NULL, SDtable=NULL,
VStable=NULL, FTtable=NULL,
FOtable=NULL, SAtable=NULL,
LEtable=NULL, OStable=NULL,
SStable=NULL, TEtable=NULL,
FMtable=NULL, BVtable=NULL, 
VDtable=NULL, SLtable=NULL, 
selected_upperHierarchy=c(1:16))

{

	SAtable_for_HL <- subset(SAtable, SAlowerHierarchy %in% c("A","B"))
	SAtable_for_CA <- subset(SAtable, SAlowerHierarchy %in% c("A","B","C"))

	if (selected_upperHierarchy %in% c(1))
	{
	FTtable <- 	FTtable[,!(colnames(FTtable) %in% c("FOid","OSid","TEid"))]
	SStable <- 	SStable[,!(colnames(SStable) %in% c("OSid","TEid","LEid"))]

	DEtable <- subset(DEtable, DEhierarchy==1) 	


	#### for HL, SL, HH

	FMSAtable <- merge(FMtable, SAtable_for_HL, by=c("SAid"),all.x=TRUE)
	FMSASStable <- merge(FMSAtable, SStable, by=c("SSid"),all.x=TRUE)
	FMSASSFOtable <- merge(FMSASStable, FOtable, by=c("FOid","FTid"),all.x=TRUE)
	FMSASSFOFTtable <- merge(FMSASSFOtable, FTtable, by=c("FTid","SDid"),all.x=TRUE)
	FMSASSFOFTVDtable <- merge(FMSASSFOFTtable, VDtable, by=c("VDid"),all.x=TRUE)
	FMSASSFOFTVDVStable <- merge(FMSASSFOFTVDtable, VStable, by=c("VSid","VDid","SDid"),all.x=TRUE)
	FMSASSFOFTVDVSSDtable <- merge(FMSASSFOFTVDVStable, SDtable, by=c("SDid"),all.x=TRUE)
	FMSASSFOFTVDVSSDDEtable <- merge(FMSASSFOFTVDVSSDtable, DEtable, by=c("DEid"),all.x=TRUE)

	FMSASSFOFTVDVSSDDEtable <- subset(FMSASSFOFTVDVSSDDEtable, DEhierarchy==selected_upperHierarchy)  #### each hierarchy sequentially?

	### determine Landing_country from harbour name - how? ###

	library(devtools)
	#install_github("ices-tools-prod/icesVocab")
	library(icesVocab)
	
	L_C <- getCodeList("Harbour_LOCODE")

	FMSASSFOFTVDVSSDDE_Location_table <- merge(FMSASSFOFTVDVSSDDEtable, 
	L_C[,c("Key","Description","LongDescription")], by.x=c("FTarrivalLocation"),by.y=c("Key"),all.x=TRUE)

	#### for CA

	BVFMSASSFOFTVDVSSDDE_Location_table <- merge(BVtable, FMSASSFOFTVDVSSDDE_Location_table, by=c("SAid","FMid"),all.x=TRUE)

	BVFMSASSFOFTVDVSSDDE_Location_table_age <- subset(BVFMSASSFOFTVDVSSDDE_Location_table, substr(BVtypeMeasured,1,3)=="Age")

	BVFMSASSFOFTVDVSSDDE_Location_table_weight <- subset(BVFMSASSFOFTVDVSSDDE_Location_table, substr(BVtypeMeasured,1,3)=="Wei")[,c("SAid","FMid",
	"BVnationalUniqueFishId","BVunitName","FMclassMeasured","BVtypeMeasured", "BVvalueMeasured", "BVvalueUnitOrScale", "BVaccuracy","BVmethod")]
	
	names(BVFMSASSFOFTVDVSSDDE_Location_table_weight)[(ncol(BVFMSASSFOFTVDVSSDDE_Location_table_weight)-4):ncol(BVFMSASSFOFTVDVSSDDE_Location_table_weight)] <- 
	paste0(names(BVFMSASSFOFTVDVSSDDE_Location_table_weight)[(ncol(BVFMSASSFOFTVDVSSDDE_Location_table_weight)-4):ncol(BVFMSASSFOFTVDVSSDDE_Location_table_weight)],
	"_weight")

	BVFMSASSFOFTVDVSSDDE_Location_table_maturity <- subset(BVFMSASSFOFTVDVSSDDE_Location_table, substr(BVtypeMeasured,1,3)=="Mat")[,c("SAid","FMid",
	"BVnationalUniqueFishId","BVunitName","FMclassMeasured","BVtypeMeasured", "BVvalueMeasured", "BVvalueUnitOrScale", "BVaccuracy","BVmethod")]

	names(BVFMSASSFOFTVDVSSDDE_Location_table_maturity)[(ncol(BVFMSASSFOFTVDVSSDDE_Location_table_maturity)-4):ncol(BVFMSASSFOFTVDVSSDDE_Location_table_maturity)] <- 
	paste0(names(BVFMSASSFOFTVDVSSDDE_Location_table_maturity)[(ncol(BVFMSASSFOFTVDVSSDDE_Location_table_maturity)-4):ncol(BVFMSASSFOFTVDVSSDDE_Location_table_maturity)],
	"_maturity")

	BVFMSASSFOFTVDVSSDDE_Location_table_sex <- subset(BVFMSASSFOFTVDVSSDDE_Location_table, substr(BVtypeMeasured,1,3)=="Sex")[,c("SAid","FMid",
	"BVnationalUniqueFishId","BVunitName","FMclassMeasured","BVtypeMeasured", "BVvalueMeasured", "BVvalueUnitOrScale", "BVaccuracy","BVmethod")]

	names(BVFMSASSFOFTVDVSSDDE_Location_table_sex)[(ncol(BVFMSASSFOFTVDVSSDDE_Location_table_sex)-4):ncol(BVFMSASSFOFTVDVSSDDE_Location_table_sex)] <- 
	paste0(names(BVFMSASSFOFTVDVSSDDE_Location_table_sex)[(ncol(BVFMSASSFOFTVDVSSDDE_Location_table_sex)-4):ncol(BVFMSASSFOFTVDVSSDDE_Location_table_sex)],
	"_sex")

		if (nrow(BVFMSASSFOFTVDVSSDDE_Location_table_weight)>0) BVFMSASSFOFTVDVSSDDE_Location_table_age <- merge(BVFMSASSFOFTVDVSSDDE_Location_table_age,
		BVFMSASSFOFTVDVSSDDE_Location_table_weight, by=c("SAid","FMid",
		"BVnationalUniqueFishId","BVunitName","FMclassMeasured"), all.x=TRUE)

		if (nrow(BVFMSASSFOFTVDVSSDDE_Location_table_maturity)>0) BVFMSASSFOFTVDVSSDDE_Location_table_age <- merge(BVFMSASSFOFTVDVSSDDE_Location_table_age,
		BVFMSASSFOFTVDVSSDDE_Location_table_maturity, by=c("SAid","FMid",
		"BVnationalUniqueFishId","BVunitName","FMclassMeasured"), all.x=TRUE)

		if (nrow(BVFMSASSFOFTVDVSSDDE_Location_table_sex)>0) BVFMSASSFOFTVDVSSDDE_Location_table_age <- merge(BVFMSASSFOFTVDVSSDDE_Location_table_age,
		BVFMSASSFOFTVDVSSDDE_Location_table_sex, by=c("SAid","FMid",
		"BVnationalUniqueFishId","BVunitName","FMclassMeasured"), all.x=TRUE)

	########################## HL ###################

	HL <- data.frame(matrix(nrow = nrow(FMSASSFOFTVDVSSDDE_Location_table), ncol=length(HL_rdb.names)), stringsAsFactors = FALSE)
	colnames(HL) <- HL_rdb.names

	HL$Record_type <- "HL"
 	HL$Sampling_type <- FMSASSFOFTVDVSSDDE_Location_table$FTsamplingType
 	HL$Landing_country <- substr(FMSASSFOFTVDVSSDDE_Location_table$FTarrivalLocation,1,2)  ## Is it a correct way to extract a Landing Country?
	HL$Vessel_flag_country <- FMSASSFOFTVDVSSDDE_Location_table$VDflagCountry
	HL$Year <- FMSASSFOFTVDVSSDDE_Location_table$DEyear
	HL$Project <- FMSASSFOFTVDVSSDDE_Location_table$DEsamplingScheme
	HL$Trip_number <- FMSASSFOFTVDVSSDDE_Location_table$FTunitName
	HL$Station_number <- FMSASSFOFTVDVSSDDE_Location_table$FOunitName
	HL$Species <- FMSASSFOFTVDVSSDDE_Location_table$SAspeciesCode
	HL$Catch_category <- FMSASSFOFTVDVSSDDE_Location_table$SAcatchCategory
	HL$Landing_category <- FMSASSFOFTVDVSSDDE_Location_table$SAlandingCategory
	HL$Comm_size_cat_scale <- FMSASSFOFTVDVSSDDE_Location_table$SAcommSizeCatScale
	HL$Comm_size_cat <- FMSASSFOFTVDVSSDDE_Location_table$SAcommSizeCat
	HL$Subsampling_category <- ''										## I didn't find this field in the RDBES
	HL$Sex <- FMSASSFOFTVDVSSDDE_Location_table$SAsex
	HL$Individual_sex <- ''
	HL$Length_class <- FMSASSFOFTVDVSSDDE_Location_table$FMclassMeasured
	HL$Number_at_length <- FMSASSFOFTVDVSSDDE_Location_table$FMnumberAtUnit
	#HL$Hierarchy <- FMSASSFOFTVDVSSDDE_Location_table$DEhierarchy				## Optional

	HL[is.na(HL)] <- ''
	HL <- aggregate(Number_at_length ~ ., data=HL, function(x) sum(x,na.rm=TRUE))


	########################## HH ###################
	
	HH <- data.frame(matrix(nrow = nrow(FMSASSFOFTVDVSSDDE_Location_table), ncol=length(HH_rdb.names)), stringsAsFactors = FALSE)
	colnames(HH) <- HH_rdb.names

	HH$Record_type <- "HH"
 	HH$Sampling_type <- FMSASSFOFTVDVSSDDE_Location_table$FTsamplingType
 	HH$Landing_country <- substr(FMSASSFOFTVDVSSDDE_Location_table$FTarrivalLocation,1,2)  ## Is it a correct way to extract a Landing Country?
	HH$Vessel_flag_country <- FMSASSFOFTVDVSSDDE_Location_table$VDflagCountry
	HH$Year <- FMSASSFOFTVDVSSDDE_Location_table$DEyear
	HH$Project <- FMSASSFOFTVDVSSDDE_Location_table$DEsamplingScheme
	HH$Trip_number <- FMSASSFOFTVDVSSDDE_Location_table$FTunitName
	HH$Station_number <- FMSASSFOFTVDVSSDDE_Location_table$FOunitName
	HH$Fishing_validity <- FMSASSFOFTVDVSSDDE_Location_table$FOvalidity
	HH$Aggregation_level <- FMSASSFOFTVDVSSDDE_Location_table$FOaggregationLevel
	HH$Catch_registration <- FMSASSFOFTVDVSSDDE_Location_table$FOcatchReg 
	HH$Species_registration <- ''										## Didn't find what is that in RDBES       
	HH$Date <- FMSASSFOFTVDVSSDDE_Location_table$FOstartDate  
	HH$Time <- FMSASSFOFTVDVSSDDE_Location_table$FOstartTime
	HH$Fishing_duration <- FMSASSFOFTVDVSSDDE_Location_table$FOduration
	HH$Pos_Start_Lat_dec <- FMSASSFOFTVDVSSDDE_Location_table$FOstartLat
	HH$Pos_Start_Lon_dec <- FMSASSFOFTVDVSSDDE_Location_table$FOstartLon
	HH$Pos_Stop_Lat_dec <- FMSASSFOFTVDVSSDDE_Location_table$FOstopLat
	HH$Pos_Stop_Lon_dec <- FMSASSFOFTVDVSSDDE_Location_table$FOstopLon
	HH$Area <- FMSASSFOFTVDVSSDDE_Location_table$FOarea
	HH$Statistical_rectangle <- FMSASSFOFTVDVSSDDE_Location_table$FOrectangle
	HH$Sub_polygon <- FMSASSFOFTVDVSSDDE_Location_table$FOjurisdictionArea
	HH$Main_fishing_depth <- FMSASSFOFTVDVSSDDE_Location_table$FOfishingDepth
	HH$Main_water_depth <- FMSASSFOFTVDVSSDDE_Location_table$FOwaterDepth
	HH$FAC_National <- FMSASSFOFTVDVSSDDE_Location_table$FOnationalFishingActivity
	HH$FAC_EC_lvl5 <- FMSASSFOFTVDVSSDDE_Location_table$FOmetier5
	HH$FAC_EC_lvl6 <- FMSASSFOFTVDVSSDDE_Location_table$FOmetier6
	HH$Gear_type <- FMSASSFOFTVDVSSDDE_Location_table$FOgear
	HH$Mesh_size <- FMSASSFOFTVDVSSDDE_Location_table$FOmeshSize
	HH$Selection_device <- FMSASSFOFTVDVSSDDE_Location_table$FOselectionDevice
	HH$Mesh_size_selection_device <- FMSASSFOFTVDVSSDDE_Location_table$FOselectionDeviceMeshSize
	#HH$Hierarchy <- FMSASSFOFTVDVSSDDE_Location_table$DEhierarchy 				## Optional

	HH[is.na(HH)] <- ''
	HH <- distinct(HH)

	########################## SL ###################
	
	SL <- data.frame(matrix(nrow = nrow(FMSASSFOFTVDVSSDDE_Location_table), ncol=length(SL_rdb.names)), stringsAsFactors = FALSE)
	colnames(SL) <- SL_rdb.names

	SL$Record_type <- "SL"
 	SL$Sampling_type <- FMSASSFOFTVDVSSDDE_Location_table$FTsamplingType
 	SL$Landing_country <- substr(FMSASSFOFTVDVSSDDE_Location_table$FTarrivalLocation,1,2)  ## Is it a correct way to extract a Landing Country?
	SL$Vessel_flag_country <- FMSASSFOFTVDVSSDDE_Location_table$VDflagCountry
	SL$Year <- FMSASSFOFTVDVSSDDE_Location_table$DEyear
	SL$Project <- FMSASSFOFTVDVSSDDE_Location_table$DEsamplingScheme
	SL$Trip_number <- FMSASSFOFTVDVSSDDE_Location_table$FTunitName
	SL$Station_number <- FMSASSFOFTVDVSSDDE_Location_table$FOunitName
	SL$Species <- FMSASSFOFTVDVSSDDE_Location_table$SAspeciesCode
	SL$Catch_category <- FMSASSFOFTVDVSSDDE_Location_table$SAcatchCategory
	SL$Landing_category <- FMSASSFOFTVDVSSDDE_Location_table$SAlandingCategory
	SL$Comm_size_cat_scale <- FMSASSFOFTVDVSSDDE_Location_table$SAcommSizeCatScale
	SL$Comm_size_cat <- FMSASSFOFTVDVSSDDE_Location_table$SAcommSizeCat
	SL$Subsampling_category <- ''										## I didn't find this field in the RDBES
	SL$Sex <- FMSASSFOFTVDVSSDDE_Location_table$SAsex
	SL$Weight <- FMSASSFOFTVDVSSDDE_Location_table$SAtotalWeightLive 
	SL$Subsample_weight <- FMSASSFOFTVDVSSDDE_Location_table$SAsampleWeightLive
	SL$Length_code <- FMSASSFOFTVDVSSDDE_Location_table$FMaccuracy
	#SL$Hierarchy <- FMSASSFOFTVDVSSDDE_Location_table$DEhierarchy 				## Optional
	
	SL[is.na(SL)] <- ''
	SL <- distinct(SL)


	########################## CA ###################
	
	CA <- data.frame(matrix(nrow = nrow(BVFMSASSFOFTVDVSSDDE_Location_table_age), ncol=length(CA_rdb.names)), stringsAsFactors = FALSE)
	colnames(CA) <- CA_rdb.names

	CA$Record_type <- "CA"
 	CA$Sampling_type <- BVFMSASSFOFTVDVSSDDE_Location_table_age$FTsamplingType
 	CA$Landing_country <- substr(BVFMSASSFOFTVDVSSDDE_Location_table_age$FTarrivalLocation,1,2)  ## Is it a correct way to extract a Landing Country?
	CA$Vessel_flag_country <- BVFMSASSFOFTVDVSSDDE_Location_table_age$VDflagCountry
	CA$Year <- BVFMSASSFOFTVDVSSDDE_Location_table_age$DEyear
	CA$Project <- BVFMSASSFOFTVDVSSDDE_Location_table_age$DEsamplingScheme
	CA$Trip_number <- BVFMSASSFOFTVDVSSDDE_Location_table_age$FTunitName
	CA$Station_number <- BVFMSASSFOFTVDVSSDDE_Location_table_age$FOunitName
	CA$Quarter <- quarter(BVFMSASSFOFTVDVSSDDE_Location_table_age$FOstartDate)
	CA$Month <- month(BVFMSASSFOFTVDVSSDDE_Location_table_age$FOstartDate)
	CA$Species <- BVFMSASSFOFTVDVSSDDE_Location_table_age$SAspeciesCode
	if (!is.null(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured_sex)) 
		CA$Sex <- ifelse(is.na(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured_sex), 
		BVFMSASSFOFTVDVSSDDE_Location_table_age$SAsex, 
		BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured_sex) else CA$Sex <- BVFMSASSFOFTVDVSSDDE_Location_table_age$SAsex
	CA$Catch_category <- BVFMSASSFOFTVDVSSDDE_Location_table_age$SAcatchCategory
	CA$Landing_category <- BVFMSASSFOFTVDVSSDDE_Location_table_age$SAlandingCategory
	CA$Comm_size_cat_scale <- BVFMSASSFOFTVDVSSDDE_Location_table_age$SAcommSizeCatScale
	CA$Comm_size_cat <- BVFMSASSFOFTVDVSSDDE_Location_table_age$SAcommSizeCat
	CA$Stock <- ''												## paste0(species_name, area)?
	CA$Area <- BVFMSASSFOFTVDVSSDDE_Location_table_age$FOarea
	CA$Statistical_rectangle <- BVFMSASSFOFTVDVSSDDE_Location_table_age$FOrectangle
	CA$Sub_polygon <- BVFMSASSFOFTVDVSSDDE_Location_table_age$FOjurisdictionArea
	CA$Length_class <- BVFMSASSFOFTVDVSSDDE_Location_table_age$FMclassMeasured
	CA$Age <- BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured
	CA$Single_fish_number <- BVFMSASSFOFTVDVSSDDE_Location_table_age$BVunitName
	CA$Length_code <- BVFMSASSFOFTVDVSSDDE_Location_table_age$FMaccuracy
	CA$Aging_method <- BVFMSASSFOFTVDVSSDDE_Location_table_age$BVmethod
	CA$Age_plus_group <- ''												## ??
	CA$Otolith_weight <- ''												## ??
	CA$Otolith_side <- ''												## ??
	if (!is.null(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured_weight))   					## BVvalueMeasured
	CA$Weight <- ifelse(is.na(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured_weight), 
		NA, BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured_weight) else CA$Weight <- NA

	if (!is.null(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVmethod_maturity))  						## BVmethod
		CA$Maturity_staging_method <- ifelse(is.na(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVmethod_maturity), 
		NA, BVFMSASSFOFTVDVSSDDE_Location_table_age$BVmethod_maturity) else CA$Maturity_staging_method <- NA

	if (!is.null(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueUnitOrScale_maturity)) 					## BVvalueUnitOrScale
		CA$Maturity_scale <- ifelse(is.na(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueUnitOrScale_maturity), 
		NA, BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueUnitOrScale_maturity) else CA$Maturity_scale <- NA

	if (!is.null(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured_maturity)) 					## BVvalueMeasured
		CA$Maturity_stage <- ifelse(is.na(BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured_maturity), 
		NA, BVFMSASSFOFTVDVSSDDE_Location_table_age$BVvalueMeasured_maturity) else CA$Maturity_stage <- NA
	
#CA$Hierarchy <- BVFMSASSFOFTVDVSSDDE_Location_table_age$DEhierarchy 				## Optional

	CA[is.na(CA)] <- ''

	} else

{
HL <- NULL
HH <- NULL
SL <- NULL
CA <- NULL
}

return(list(HL = HL, HH = HH, SL = SL, CA = CA))
}


