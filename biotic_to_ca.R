biotic_to_ca <- function(CA, HL, by=list(), lengthColumn = "Length", bioLength = lengthColumn,
                         catchLength = lengthColumn, catchNumber="NumberAtLength") {
  origCA <- CA;
  colnames(origCA)[which(names(origCA) == bioLength)] <- "Length"
  #Column handling
  by <- unlist(by);
  #Count CA records
  CAdist <- as.data.frame(table(CA[c(by,bioLength)]),stringsAsFactors = FALSE);
  names(CAdist) <- c(by, "Length", "CA_N");
  CAdist$Length <- as.numeric(CAdist$Length);
  CAdist$CA_N <- as.numeric(CAdist$CA_N);
  #Select HL records
  HLdist <- HL[c(by,catchLength, catchNumber)];
  names(HLdist) <- c(by, "Length", "HL_N");
  HLdist$Length <- as.numeric(HLdist$Length);
  HLdist$HL_N <- as.numeric(HLdist$HL_N);
  #Combine records
  dist <- merge(CAdist, HLdist);
  #NA cleanup
  #Calculate differences
  dist$diff <- dist$HL_N - dist$CA_N;
  #Generate CA records
  resultKey <- c(by, "Length");
  for(row in 1:nrow(dist)) {
    if(dist[row,"diff"] == 0) {
      next;
    }
    diff <- dist[row,"diff"];
    for(r in 1:diff) {
      origCA[nrow(origCA)+1,] <- NA;
      origCA[nrow(origCA), resultKey] <- dist[row,resultKey];
    }
  }
  return(origCA);
};
