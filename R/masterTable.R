#### input: RDBES object
#### output: list of the following master tables: sampleData (all tables in RDBES object hierarchy,
#### from SA to DE, joined together), lengthData (all tables in RDBES object hierarchy,
#### from FM to DE, joined together), biologicalData (all tables in RDBES object hierarchy,
#### from BV to DE, joined together).
#### Note: the returned lengthData is an expanded table of FM according to the column FMnumAtUnit,
#### so that 1 line = 1 length measured individual

masterTable <- function(RDBESObject)
{

  require(dplyr)
  require(purrr)
  require(RDBEScore)

  
  hierarchy <- unique(RDBESObject$DE$DEhierarchy)
  
  tablesList <- RDBESObject[getTablesInRDBESHierarchy(hierarchy)]
  
  
  ### delete all "id"-columns with NA
  tablesListCleaned <- lapply(tablesList, function(tbl) {
    tbl %>% 
      dplyr::select(-(
        ends_with("id") & 
          where(function(x) all(is.na(na_if(as.character(x), ""))))
      ))
  })
  
  
  ### number of SA table
  SAindex <- which(names(tablesListCleaned) == "SA")
  
  ### sample master table
  sample_master_table <- suppressMessages(tablesListCleaned[1:SAindex]%>%
    rev()%>%
    reduce(left_join))
  
  ### length_master_table
  length_master_table <- tablesListCleaned$FM%>%
    merge(sample_master_table, by=c("SAid"), all.x = TRUE)
  
  ### bio_master_table
  bio_master_table <- tablesListCleaned$BV%>%
    merge(length_master_table, by=c("FMid"), all.x=TRUE)%>%
    mutate(n = 1)
  
  ### expand table: 1 line = 1 individual
  length_master_table <- 
    uncount(length_master_table, length_master_table$FMnumAtUnit)%>%
    mutate(n = 1)
  
  
  return(list(lengthData = length_master_table, 
              biologicalData = bio_master_table, 
              sampleData = sample_master_table))
  
}
