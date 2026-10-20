masterTable <- function(RDBESObject)
{
  
  hierarchy <- unique(RDBESObject$DE$DEhierarchy)
  
  if (any(hierarchy %in% c(1, 5)))
    
  {
    
    RDBESObject1_5 <- filterRDBESDataObject(RDBESObject,
                                            fieldsToFilter = c("DEhierarchy"),
                                            valuesToFilter = c("1", "5"),
                                            killOrphans = TRUE)
    
    RDBESObject1_5$BV <- RDBESObject1_5$BV%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), "")))))) 
    RDBESObject1_5$FM <- RDBESObject1_5$FM%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))  
    RDBESObject1_5$SA <- RDBESObject1_5$SA%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))
    RDBESObject1_5$SS <- RDBESObject1_5$SS%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))
    RDBESObject1_5$SL <- RDBESObject1_5$SL%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))
    RDBESObject1_5$FO <- RDBESObject1_5$FO%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))
    RDBESObject1_5$FT <- RDBESObject1_5$FT%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))
    RDBESObject1_5$VD <- RDBESObject1_5$VD%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))
    RDBESObject1_5$VS <- RDBESObject1_5$VS%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))
    RDBESObject1_5$SD <- RDBESObject1_5$SD%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))
    RDBESObject1_5$DE <- RDBESObject1_5$DE%>%dplyr::select(-(ends_with("id") & where(~ all(is.na(na_if(as.character(.x), ""))))))
    
    sample_master_table <- 
      RDBESObject1_5$SA%>%
      merge(RDBESObject1_5$SS, by=c("SSid"),all.x=TRUE)%>%
      merge(RDBESObject1_5$SL, by=c("SLid"),all.x=TRUE)%>%
      merge(RDBESObject1_5$FO, by=c("FOid"),all.x=TRUE)%>%
      merge(RDBESObject1_5$FT, by=c("FTid"),all.x=TRUE)%>%
      merge(RDBESObject1_5$VD, by=c("VDid"),all.x=TRUE)%>%
      merge(RDBESObject1_5$VS, by=c("VSid","VDid"),all.x=TRUE)%>%
      merge(RDBESObject1_5$SD, by=c("SDid"),all.x=TRUE)%>%
      merge(RDBESObject1_5$DE, by=c("DEid"),all.x=TRUE)
    
    length_master_table <- RDBESObject1_5$FM%>%
      merge(sample_master_table, by=c("SAid"), all.x = TRUE)
    
    bio_master_table <- RDBESObject1_5$BV%>%
      merge(length_master_table, by=c("FMid"), all.x=TRUE)%>%
      mutate(n = 1)
    
    
    length_master_table <- 
      uncount(length_master_table, length_master_table$FMnumAtUnit)%>%
      mutate(n = 1)
    
  }
  
  return(list(lengthData = length_master_table, 
              biologicalData = bio_master_table, 
              sampleData = sample_master_table))
  
}
