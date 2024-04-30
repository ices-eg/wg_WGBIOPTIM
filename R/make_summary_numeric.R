#' Summary statistics calculation: mean length, se, median, min, max,  n classes sampled
#' @param df1 data frame of sampling data
#' @param variable "lenCls"
#' @param a coefficient of length-weight relationship
#' @param b coefficient of length-weight relationship
#' @return table reporting the different estimates
#' @export
#' @examples make_summary_numeric(example_samples,"lenCls",a=0.0006,b=3)
#'
make_summary_numeric<-function(df1, variable, a, b){
  # from the script of Nuno Prista 2017

  x <- df1[,variable]

   if(is.vector(x) | is.factor(x[variable])) {x<-data.frame(dummy=1,x); colnames(x)[2]<-variable}
   # summary
  n_indiv = nrow(x)
  #NAs_x<-sum(is.na(x[,variable]))
  mean_x<-mean(x[,variable], na.rm=TRUE)
  stand_err_mean_x<-sd(x[,variable], na.rm=TRUE)/sqrt(length(x[,variable]))
  median_x<-median(x[,variable], na.rm=TRUE)
  min_x<- min(x[,variable], na.rm=TRUE)
  max_x<-max(x[,variable], na.rm=TRUE)
  n_class_sampled_x<-length(unique(x[!is.na(x[,variable]),variable]))
  #ATTT!!!!
  if(variable == "lenCls") estim_weight_sample <- sum(exp(a)*x$lenCls^b) else {estim_weight_sample<-NA}
  # mean weighed CV

  # output
  DF<-data.frame(var=variable, n = n_indiv, estim_weight = estim_weight_sample,
             mean=mean_x, se = stand_err_mean_x, cv=round(stand_err_mean_x/mean_x*100,1),
             min = min_x, median = median_x, max = max_x, n_class_sampled=n_class_sampled_x)
return(DF)
  }
