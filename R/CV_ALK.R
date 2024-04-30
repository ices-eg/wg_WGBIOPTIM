#' Function to calculate CV on ALK
#'
#' @param DF dataframe with lengths in the firts column and ages on the other columns
#'
#' @return CV by age class and total CV
#' @export
#'
#' @examples CV_ALK(example_ageLength)

CV_ALK <- function (DF) {


  DF <- DF[, colnames(DF) != "total_ages"]

  DFrame=data.frame(ages=colnames(DF)[-c(1,ncol(DF))],CV=rep(999,(ncol(DF)-2)), Nipi=rep(999,(ncol(DF)-2)), VarNipi=rep(999,(ncol(DF)-2)))
  ages=colnames(DF)[-c(1,ncol(DF))]
  for (ag in 1:length(ages)) {
    age_temp=DF[,c(1,(1+ag))]
    age_temp[which(is.na(age_temp[,2])),2] <- 0
    age_temp$Ni=DF$total_lengths
    if (ncol(DF)>3){
      age_temp$ni=rowSums(DF[,-c(1,ncol(DF))],na.rm=T)
    } else {
      age_temp$ni= DF[,2]
    }
    age_temp$nipi=age_temp[,2]
    age_temp$pi=age_temp$nipi/age_temp$ni
    age_temp$Nipi=age_temp$Ni*age_temp$pi
    age_temp$Varpi = (age_temp$pi*(1-age_temp$pi)) / age_temp$ni
    age_temp$VarNipi=(age_temp$Ni^2) * age_temp$Varpi
    DFrame$Nipi[ag]=sum(age_temp$Nipi,na.rm=T)
    DFrame$VarNipi[ag]=sum(age_temp$VarNipi,na.rm=T)
    DFrame$CV[ag]=sum(age_temp$VarNipi,na.rm=T)^0.5/sum(age_temp$Nipi,na.rm=T) *100

  }

  DFrame$total_CV = ifelse(sum( DFrame$Nipi) == 0, 0, (sum( DFrame$VarNipi)^0.5)/sum( DFrame$Nipi) ) * 100

  return(DFrame[,c(1,2,5)])
}
