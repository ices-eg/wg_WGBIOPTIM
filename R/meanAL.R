#' Mean age-at-length
#'
#' @param df A data frame containing a column named `Length` and age columns starting with "Age".
#'
#' @return A data frame containing `Length` and `mean_age` columns.
#'
#' @export
#'
#' @examples
#' mean_age_at_length(example_ageLength)


mean_age_at_length <- function(df) {
  mean_ages <- numeric(nrow(df))                                 
  
  # Determine the number of age columns
  age_columns <- colnames(df)[grepl("^Age", colnames(df))]       
  num_ages <- length(age_columns)
  
  # calculate mean age at length for each row (length class)
  for (i in 1:nrow(df)) {                                       
    # Extract ages and their counts
    ages <- as.numeric(df[i, age_columns])
    counts <- ages * (0:(num_ages - 1))
    
    # Calculate mean age
    if (sum(ages) > 0) {                                        
      mean_ages[i] <- sum(counts) / sum(ages)
    } else {
      mean_ages[i] <- NA
    }
  }
  
  # Create a new data frame with Length and mean_age
  result_df <- data.frame(Length = df$Length, mean_age = mean_ages)
  
  return(result_df)
}