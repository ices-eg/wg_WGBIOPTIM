#' MWCV (mean weigthed CV) calculation
#' @param df1 data frame of sampled data by length class in CA format (RDB) (individual measurements)
#' @param variable "lenCls"
#' @return MWCV
#' @export
#' @importFrom stats median sd
#' @examples
#'MWCV(example_samples,"lenCls")
MWCV<-function(df1, variable){


  x <- df1[,variable]

  if(is.vector(x) | is.factor(x[variable])) {x<-data.frame(dummy=1,x); colnames(x)[2]<-variable}

  # mean weighed CV
  sigma_i<-sqrt(nrow(x)*as.matrix(prop.table(table(x[,variable]))*(1-prop.table(table(x[,variable])))))
  cv_i <- sigma_i / (nrow(x)*as.matrix(prop.table(table(x[,variable]))))
  MWCV<-round(sum(sigma_i)/nrow(x)*100,1)

  # output

  return(MWCV)
}
