
#' dataRDBEStoRDB
#'
#' @param x zip file containi RDBES input files
#' @param selected_upperHierarchy upper hierarchy
#'
#' @return list(HL = HL, HH = HH, SL = SL, CA = CA) RDB tables
#' @export
#'
#' @examples
#'dataRDBEStoRDB(list.rdbes, 1)
#'
dataRDBEStoRDB <- function(x, selected_upperHierarchy=1)

{

  require(dplyr)
  require(lubridate)


  #not_all_ids_na <- function(x) any(!is.na(x))

  SAtable_for_length <- subset(x$Sample, SAlowerHierarchy %in% c("A","B"))
  SAtable_for_age <- subset(x$Sample, SAlowerHierarchy %in% c("A","B","C"))

  DEtable <- subset(x$Design, DEhierarchy==selected_upperHierarchy)

  if (selected_upperHierarchy %in% c(1,2))
  {
    FTtable <- x$FishingTrip[,!(names(x$FishingTrip) %in% c("LEid","TEid", "FOid", "OSid", "SDid"))]
    FOtable <- x$FishingOperation[,!(names(x$FishingOperation) %in% c("SDid"))]
    SStable <- x$SpeciesSelection[,!(names(x$SpeciesSelection) %in% c("LEid","TEid", "FTid", "OSid"))]
    FMtable <- x$FrequencyMeasure
    VDtable <- x$VesselDetails
    VStable <- x$VesselSelection
    BVtable <- x$BiologicalVariable[,!(names(x$BiologicalVariable) %in% c("SAid"))]
    SDtable <- x$SamplingDetails
    DEtable <- x$Design


    #### for length

    FMSAtable <- merge(FMtable, SAtable_for_length, by=c("SAid"),all.x=TRUE)
    FMSASStable <- merge(FMSAtable, SStable, by=c("SSid"),all.x=TRUE)
    FMSASSFOtable <- merge(FMSASStable, FOtable, by=c("FOid"),all.x=TRUE)
    FMSASSFOFTtable <- merge(FMSASSFOtable, FTtable, by=c("FTid"),all.x=TRUE)
    FMSASSFOFTVDtable <- merge(FMSASSFOFTtable, VDtable, by=c("VDid"),all.x=TRUE)
    FMSASSFOFTVDVStable <- merge(FMSASSFOFTVDtable, VStable, by=c("VSid","VDid"),all.x=TRUE)
    FMSASSFOFTVDVSSDtable <- merge(FMSASSFOFTVDVStable, SDtable, by=c("SDid"),all.x=TRUE)
    FMSASSFOFTVDVSSDDEtable <- merge(FMSASSFOFTVDVSSDtable, DEtable, by=c("DEid"),all.x=TRUE)


    #### for CA

    BVFMSASSFOFTVDVSSDDEtable <- merge(BVtable, FMSASSFOFTVDVSSDDEtable, by=c("FMid"),all.x=TRUE)

    BVFMSASSFOFTVDVSSDDEtable_age <- subset(BVFMSASSFOFTVDVSSDDEtable, substr(BVtypeMeasured,1,3)=="Age")

    BVFMSASSFOFTVDVSSDDEtable_weight <- subset(BVFMSASSFOFTVDVSSDDEtable, substr(BVtypeMeasured,1,3)=="Wei")[,c("SAid","FMid",
                                                                                                                "BVnationalUniqueFishId","BVunitName","FMclassMeasured","BVtypeMeasured", "BVvalueMeasured", "BVvalueUnitOrScale", "BVaccuracy","BVmethod")]

    names(BVFMSASSFOFTVDVSSDDEtable_weight)[(ncol(BVFMSASSFOFTVDVSSDDEtable_weight)-4):ncol(BVFMSASSFOFTVDVSSDDEtable_weight)] <-
      paste0(names(BVFMSASSFOFTVDVSSDDEtable_weight)[(ncol(BVFMSASSFOFTVDVSSDDEtable_weight)-4):ncol(BVFMSASSFOFTVDVSSDDEtable_weight)],
             "_weight")

    BVFMSASSFOFTVDVSSDDEtable_maturity <- subset(BVFMSASSFOFTVDVSSDDEtable, substr(BVtypeMeasured,1,3)=="Mat")[,c("SAid","FMid",
                                                                                                                  "BVnationalUniqueFishId","BVunitName","FMclassMeasured","BVtypeMeasured", "BVvalueMeasured", "BVvalueUnitOrScale", "BVaccuracy","BVmethod")]

    names(BVFMSASSFOFTVDVSSDDEtable_maturity)[(ncol(BVFMSASSFOFTVDVSSDDEtable_maturity)-4):ncol(BVFMSASSFOFTVDVSSDDEtable_maturity)] <-
      paste0(names(BVFMSASSFOFTVDVSSDDEtable_maturity)[(ncol(BVFMSASSFOFTVDVSSDDEtable_maturity)-4):ncol(BVFMSASSFOFTVDVSSDDEtable_maturity)],
             "_maturity")

    BVFMSASSFOFTVDVSSDDEtable_sex <- subset(BVFMSASSFOFTVDVSSDDEtable, substr(BVtypeMeasured,1,3)=="Sex")[,c("SAid","FMid",
                                                                                                             "BVnationalUniqueFishId","BVunitName","FMclassMeasured","BVtypeMeasured", "BVvalueMeasured", "BVvalueUnitOrScale", "BVaccuracy","BVmethod")]

    names(BVFMSASSFOFTVDVSSDDEtable_sex)[(ncol(BVFMSASSFOFTVDVSSDDEtable_sex)-4):ncol(BVFMSASSFOFTVDVSSDDEtable_sex)] <-
      paste0(names(BVFMSASSFOFTVDVSSDDEtable_sex)[(ncol(BVFMSASSFOFTVDVSSDDEtable_sex)-4):ncol(BVFMSASSFOFTVDVSSDDEtable_sex)],
             "_sex")

    if (nrow(BVFMSASSFOFTVDVSSDDEtable_weight)>0) BVFMSASSFOFTVDVSSDDEtable_age <- merge(BVFMSASSFOFTVDVSSDDEtable_age,
                                                                                         BVFMSASSFOFTVDVSSDDEtable_weight, by=c("SAid","FMid",
                                                                                                                                "BVnationalUniqueFishId","BVunitName","FMclassMeasured"), all.x=TRUE)

    if (nrow(BVFMSASSFOFTVDVSSDDEtable_maturity)>0) BVFMSASSFOFTVDVSSDDEtable_age <- merge(BVFMSASSFOFTVDVSSDDEtable_age,
                                                                                           BVFMSASSFOFTVDVSSDDEtable_maturity, by=c("SAid","FMid",
                                                                                                                                    "BVnationalUniqueFishId","BVunitName","FMclassMeasured"), all.x=TRUE)

    if (nrow(BVFMSASSFOFTVDVSSDDEtable_sex)>0) BVFMSASSFOFTVDVSSDDEtable_age <- merge(BVFMSASSFOFTVDVSSDDEtable_age,
                                                                                      BVFMSASSFOFTVDVSSDDEtable_sex, by=c("SAid","FMid",
                                                                                                                          "BVnationalUniqueFishId","BVunitName","FMclassMeasured"), all.x=TRUE)

    ########################## HL ###################

    HL <- data.frame(matrix(nrow = nrow(FMSASSFOFTVDVSSDDEtable),
                            ncol=length(HL_rdb.names)), stringsAsFactors = FALSE)
    colnames(HL) <- HL_rdb.names

    HL$Record_type <- "HL"
    HL$Sampling_type <- FMSASSFOFTVDVSSDDEtable$FTsampType
    HL$Landing_country <- substr(FMSASSFOFTVDVSSDDEtable$FTarrivalLocation,1,2)  ## Is it a correct way to extract a Landing Country?
    HL$Vessel_flag_country <- FMSASSFOFTVDVSSDDEtable$SDctry
    HL$Year <- FMSASSFOFTVDVSSDDEtable$DEyear
    HL$Project <- FMSASSFOFTVDVSSDDEtable$DEsampScheme
    HL$Trip_number <- FMSASSFOFTVDVSSDDEtable$FTunitName
    HL$Station_number <- FMSASSFOFTVDVSSDDEtable$FOunitName
    HL$Species <- FMSASSFOFTVDVSSDDEtable$SAspeciesCode
    HL$Catch_category <- FMSASSFOFTVDVSSDDEtable$SAcatchCategory
    HL$Landing_category <- FMSASSFOFTVDVSSDDEtable$SAlandingCategory
    HL$Comm_size_cat_scale <- FMSASSFOFTVDVSSDDEtable$SAcommSizeCatScale
    HL$Comm_size_cat <- FMSASSFOFTVDVSSDDEtable$SAcommSizeCat
    HL$Subsampling_category <- ''								## I didn't find this field in the RDBES
    HL$Sex <- FMSASSFOFTVDVSSDDEtable$SAsex
    HL$Individual_sex <- ''
    HL$Length_class <- FMSASSFOFTVDVSSDDEtable$FMclassMeasured
    HL$Number_at_length <- FMSASSFOFTVDVSSDDEtable$FMnumberAtUnit
    #HL$Hierarchy <- FMSASSFOFTVDVSSDDEtable$DEhierarchy				## Optional

    HL[is.na(HL)] <- ''
    HL <- aggregate(Number_at_length ~ ., data=HL, function(x) sum(x,na.rm=TRUE))


    ########################## HH ###################

    HH <- data.frame(matrix(nrow = nrow(FMSASSFOFTVDVSSDDEtable), ncol=length(HH_rdb.names)), stringsAsFactors = FALSE)
    colnames(HH) <- HH_rdb.names

    HH$Record_type <- "HH"
    HH$Sampling_type <- FMSASSFOFTVDVSSDDEtable$FTsamplingType
    HH$Landing_country <- substr(FMSASSFOFTVDVSSDDEtable$FTarrivalLocation,1,2)  ## Is it a correct way to extract a Landing Country?
    HH$Vessel_flag_country <- FMSASSFOFTVDVSSDDEtable$VDflagCountry
    HH$Year <- FMSASSFOFTVDVSSDDEtable$DEyear
    HH$Project <- FMSASSFOFTVDVSSDDEtable$DEsamplingScheme
    HH$Trip_number <- FMSASSFOFTVDVSSDDEtable$FTunitName
    HH$Station_number <- FMSASSFOFTVDVSSDDEtable$FOunitName
    HH$Fishing_validity <- FMSASSFOFTVDVSSDDEtable$FOvalidity
    HH$Aggregation_level <- FMSASSFOFTVDVSSDDEtable$FOaggregationLevel
    HH$Catch_registration <- FMSASSFOFTVDVSSDDEtable$FOcatchReg
    HH$Species_registration <- ''										## Didn't find what is that in RDBES
    HH$Date <- FMSASSFOFTVDVSSDDEtable$FOstartDate
    HH$Time <- FMSASSFOFTVDVSSDDEtable$FOstartTime
    HH$Fishing_duration <- FMSASSFOFTVDVSSDDEtable$FOduration
    HH$Pos_Start_Lat_dec <- FMSASSFOFTVDVSSDDEtable$FOstartLat
    HH$Pos_Start_Lon_dec <- FMSASSFOFTVDVSSDDEtable$FOstartLon
    HH$Pos_Stop_Lat_dec <- FMSASSFOFTVDVSSDDEtable$FOstopLat
    HH$Pos_Stop_Lon_dec <- FMSASSFOFTVDVSSDDEtable$FOstopLon
    HH$Area <- FMSASSFOFTVDVSSDDEtable$FOarea
    HH$Statistical_rectangle <- FMSASSFOFTVDVSSDDEtable$FOrectangle
    HH$Sub_polygon <- FMSASSFOFTVDVSSDDEtable$FOjurisdictionArea
    HH$Main_fishing_depth <- FMSASSFOFTVDVSSDDEtable$FOfishingDepth
    HH$Main_water_depth <- FMSASSFOFTVDVSSDDEtable$FOwaterDepth
    HH$FAC_National <- FMSASSFOFTVDVSSDDEtable$FOnationalFishingActivity
    HH$FAC_EC_lvl5 <- FMSASSFOFTVDVSSDDEtable$FOmetier5
    HH$FAC_EC_lvl6 <- FMSASSFOFTVDVSSDDEtable$FOmetier6
    HH$Gear_type <- FMSASSFOFTVDVSSDDEtable$FOgear
    HH$Mesh_size <- FMSASSFOFTVDVSSDDEtable$FOmeshSize
    HH$Selection_device <- FMSASSFOFTVDVSSDDEtable$FOselectionDevice
    HH$Mesh_size_selection_device <- FMSASSFOFTVDVSSDDEtable$FOselectionDeviceMeshSize
    #HH$Hierarchy <- FMSASSFOFTVDVSSDDEtable$DEhierarchy 				## Optional

    HH[is.na(HH)] <- ''
    HH <- distinct(HH)

    ########################## SL ###################

    SL <- data.frame(matrix(nrow = nrow(FMSASSFOFTVDVSSDDEtable), ncol=length(SL_rdb.names)), stringsAsFactors = FALSE)
    colnames(SL) <- SL_rdb.names

    SL$Record_type <- "SL"
    SL$Sampling_type <- FMSASSFOFTVDVSSDDEtable$FTsamplingType
    SL$Landing_country <- substr(FMSASSFOFTVDVSSDDEtable$FTarrivalLocation,1,2)  ## Is it a correct way to extract a Landing Country?
    SL$Vessel_flag_country <- FMSASSFOFTVDVSSDDEtable$VDflagCountry
    SL$Year <- FMSASSFOFTVDVSSDDEtable$DEyear
    SL$Project <- FMSASSFOFTVDVSSDDEtable$DEsamplingScheme
    SL$Trip_number <- FMSASSFOFTVDVSSDDEtable$FTunitName
    SL$Station_number <- FMSASSFOFTVDVSSDDEtable$FOunitName
    SL$Species <- FMSASSFOFTVDVSSDDEtable$SAspeciesCode
    SL$Catch_category <- FMSASSFOFTVDVSSDDEtable$SAcatchCategory
    SL$Landing_category <- FMSASSFOFTVDVSSDDEtable$SAlandingCategory
    SL$Comm_size_cat_scale <- FMSASSFOFTVDVSSDDEtable$SAcommSizeCatScale
    SL$Comm_size_cat <- FMSASSFOFTVDVSSDDEtable$SAcommSizeCat
    SL$Subsampling_category <- ''										## I didn't find this field in the RDBES
    SL$Sex <- FMSASSFOFTVDVSSDDEtable$SAsex
    SL$Weight <- FMSASSFOFTVDVSSDDEtable$SAtotalWeightLive
    SL$Subsample_weight <- FMSASSFOFTVDVSSDDEtable$SAsampleWeightLive
    SL$Length_code <- FMSASSFOFTVDVSSDDEtable$FMaccuracy
    #SL$Hierarchy <- FMSASSFOFTVDVSSDDEtable$DEhierarchy 				## Optional

    SL[is.na(SL)] <- ''
    SL <- distinct(SL)


    ########################## CA ###################

    CA <- data.frame(matrix(nrow = nrow(BVFMSASSFOFTVDVSSDDEtable_age), ncol=length(CA_rdb.names)), stringsAsFactors = FALSE)
    colnames(CA) <- CA_rdb.names

    CA$Record_type <- "CA"
    CA$Sampling_type <- BVFMSASSFOFTVDVSSDDEtable_age$FTsamplingType
    CA$Landing_country <- substr(BVFMSASSFOFTVDVSSDDEtable_age$FTarrivalLocation,1,2)  ## Is it a correct way to extract a Landing Country?
    CA$Vessel_flag_country <- BVFMSASSFOFTVDVSSDDEtable_age$VDflagCountry
    CA$Year <- BVFMSASSFOFTVDVSSDDEtable_age$DEyear
    CA$Project <- BVFMSASSFOFTVDVSSDDEtable_age$DEsamplingScheme
    CA$Trip_number <- BVFMSASSFOFTVDVSSDDEtable_age$FTunitName
    CA$Station_number <- BVFMSASSFOFTVDVSSDDEtable_age$FOunitName
    CA$Quarter <- ifelse(as.numeric(substr(BVFMSASSFOFTVDVSSDDEtable_age$FOstartDate, 6, 7)) %in% 1:3, 1,
                         ifelse(as.numeric(substr(BVFMSASSFOFTVDVSSDDEtable_age$FOstartDate, 6, 7)) %in% 4:6, 2,
                                ifelse(as.numeric(substr(BVFMSASSFOFTVDVSSDDEtable_age$FOstartDate, 6, 7)) %in% 7:9, 3, 4)))
    CA$Month <- as.numeric(substr(BVFMSASSFOFTVDVSSDDEtable_age$FOstartDate, 6, 7))
    CA$Species <- BVFMSASSFOFTVDVSSDDEtable_age$SAspeciesCode
    if (!is.null(BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured_sex))
      CA$Sex <- ifelse(is.na(BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured_sex),
                       BVFMSASSFOFTVDVSSDDEtable_age$SAsex,
                       BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured_sex) else CA$Sex <- BVFMSASSFOFTVDVSSDDEtable_age$SAsex
    CA$Catch_category <- BVFMSASSFOFTVDVSSDDEtable_age$SAcatchCategory
    CA$Landing_category <- BVFMSASSFOFTVDVSSDDEtable_age$SAlandingCategory
    CA$Comm_size_cat_scale <- BVFMSASSFOFTVDVSSDDEtable_age$SAcommSizeCatScale
    CA$Comm_size_cat <- BVFMSASSFOFTVDVSSDDEtable_age$SAcommSizeCat
    CA$Stock <- ''												## paste0(species_name, area)?
    CA$Area <- BVFMSASSFOFTVDVSSDDEtable_age$FOarea
    CA$Statistical_rectangle <- BVFMSASSFOFTVDVSSDDEtable_age$FOrectangle
    CA$Sub_polygon <- BVFMSASSFOFTVDVSSDDEtable_age$FOjurisdictionArea
    CA$Length_class <- BVFMSASSFOFTVDVSSDDEtable_age$FMclassMeasured
    CA$Age <- BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured
    CA$Single_fish_number <- BVFMSASSFOFTVDVSSDDEtable_age$BVunitName
    CA$Length_code <- BVFMSASSFOFTVDVSSDDEtable_age$FMaccuracy
    CA$Aging_method <- BVFMSASSFOFTVDVSSDDEtable_age$BVmethod
    CA$Age_plus_group <- ''												## ??
    CA$Otolith_weight <- ''												## ??
    CA$Otolith_side <- ''												## ??
    if (!is.null(BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured_weight))   					## BVvalueMeasured
      CA$Weight <- ifelse(is.na(BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured_weight),
                          NA, BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured_weight) else CA$Weight <- NA

    if (!is.null(BVFMSASSFOFTVDVSSDDEtable_age$BVmethod_maturity))  						## BVmethod
      CA$Maturity_staging_method <- ifelse(is.na(BVFMSASSFOFTVDVSSDDEtable_age$BVmethod_maturity),
                                           NA, BVFMSASSFOFTVDVSSDDEtable_age$BVmethod_maturity) else CA$Maturity_staging_method <- NA

    if (!is.null(BVFMSASSFOFTVDVSSDDEtable_age$BVvalueUnitOrScale_maturity)) 					## BVvalueUnitOrScale
      CA$Maturity_scale <- ifelse(is.na(BVFMSASSFOFTVDVSSDDEtable_age$BVvalueUnitOrScale_maturity),
                                  NA, BVFMSASSFOFTVDVSSDDEtable_age$BVvalueUnitOrScale_maturity) else CA$Maturity_scale <- NA

    if (!is.null(BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured_maturity)) 					## BVvalueMeasured
      CA$Maturity_stage <- ifelse(is.na(BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured_maturity),
                                  NA, BVFMSASSFOFTVDVSSDDEtable_age$BVvalueMeasured_maturity) else CA$Maturity_stage <- NA

    #CA$Hierarchy <- BVFMSASSFOFTVDVSSDDEtable_age$DEhierarchy 				## Optional

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

