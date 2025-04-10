#' Parameters of the von Bertalanffy growth model
#' 
#' This function fits the von Bertalanffy growth model to a given dataset using nonlinear least squares optimization.
#'
#' @param dataV A data frame containing the data. It must include columns for the response variable (lengthClass) and the predictor variable (age). Ensure there are no NA values in these columns.
#'
#' @param lengthClass The name of the column in dataV representing the observed length at age.
#' @param age The name of the column in dataV representing the age of the organism.
#' 
#' @param Linf Initial value for the asymptotic length parameter need to be provided.
#' @param K Initial value for the growth coefficient parameter need to be provided.
#' @param t0 Initial value for the theoretical age at length zero need to be provided.
#' @param maxiter Maximum number of iterations for the nonlinear least squares fitting process. Default is 10,000.
#'
#' @return An object of class nls representing the fitted model. This object contains the optimized values for the parameters Linf, K, and t0, along with other details about the model fit.



fitVonBertalanffy <- function(dataV, lengthClass, age, Linf, K, t0, maxiter = 10000) {
  
  # von Bertalanffy growth model formula
  vbTypical <- as.formula(paste(lengthClass, "~ Linf * (1 - exp(-K * (", age, " - t0)))"))
  
  
  # initial parameter values
  svTypical <- list(Linf = Linf, K = K, t0 = t0)
  
  # control parameters for the nls function
  control <- nls.control(maxiter = maxiter)
  
  # nonlinear least squares fitting
  fitTypical <- nls(vbTypical, data = dataV, start = svTypical, control = control)
  
  return(fitTypical)
}
