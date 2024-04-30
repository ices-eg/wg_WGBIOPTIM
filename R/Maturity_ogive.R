

#' Maturity ogive parameters estimation (L25, L50, L75)
#'
#' @param data  Dataframe with information from the simulations, containning thefollowing variables: length, maturity (0- immature; 1  - mature); ID_sim (identification of the number of the simulation run); type (number of individuals selected in the current simulation).
#'
#' @return L25 (length at which 25% of the individuals are mature), L50 (length at which 50% of the individuals are mature), L75 (length at which 75% of the individuals are mature)
#' @export
#'
#' @examples
#' Maturity_ogive(simulation_biodata)
#'
#' @importFrom stats aggregate binomial glm median sd
#'
Maturity_ogive <- function (data){

  sim<-unique(data$ID_sim)  ## identification of the number of the simulation run
  results<-matrix(nrow=length(sim),ncol=6)
  for(nb in 1: length(sim))
  {
    glm1 <- glm(factor(maturity)~Lt,family=binomial,data=data[data$ID_sim==nb,])
    Lmat <- signif(dose.p (glm1, p = c(0.25, 0.50, 0.75)), digits = 3)
    results[nb,1]<-unique(data$year[data$ID_sim==nb]) #year
    results[nb,2]<-as.numeric(Lmat[[1]]) #L25 - length at which 25% of the individuals are mature
    results[nb,3]<-as.numeric(Lmat[[2]]) #L50 - length at which 50% of the individuals are mature
    results[nb,4]<-as.numeric(Lmat[[3]]) #L75 - length at which 75% of the individuals are mature
    results[nb,5]<-nb ##variable corresponding to the number of the simulation run (ID_sim)
    results[nb,6]<-unique(data$type) ##numOtolitsPerClass
  }
  colnames(results)<-c("year","L25","L50","L75","ID_sim","type") ## type is the simulation stratification option, e.g. number of individuals selected by                                                                      ## length class (numOtolitsPerClass)
  return(results)
}
