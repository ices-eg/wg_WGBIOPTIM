

#' Assign the sampled individuals as immature or mature according with the maturity stage
#'
#' @param dataBio Dataframe with the varable maturity stage included, this should a numeric variable
#' @param param.stage_mature This is the first maturity stage of the mature fish
#'
#' @return A dataframe with the variable "mat_stg" included, with only 0s and 1s. The 0 if the individual is immature and as 1 if the individual is mature.
#' @export
#'
#' @examples
#' mature(dataBio)
mature <- function(dataBio, param.stage_mature) {
  ifelse(dataBio$mat_stg<param.stage_mature, 0, ifelse(dataBio$mat_stg>=param.stage_mature,1,NA))
  return(dataBio)
}


