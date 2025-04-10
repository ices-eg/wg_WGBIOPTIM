#' Calculate CV per length class and total CV
#'
#' @param df A data frame containing a column named `Length` and variable columns starting with the specified prefix.
#' @param variable A string specifying the prefix of the variable columns to calculate CV for (e.g., "Age", "Weight").
#'
#' @return A data frame containing `Length`, CV per length class, and total CV.
#'
#' @export
#'
#' @examples
#' cv_length(example_ageLength, "Age")


cv_length <- function(df, variable) {
  # vectors to store CV per length class
  cv_per_length <- numeric(nrow(df))
  
  # Determine the number of variable columns
  variable_columns <- colnames(df)[grepl(paste0("^", variable), colnames(df))]
  
  # Initialize a data frame to store intermediate calculations
  DFrame <- data.frame(Length = df$Length, CV = rep(999, nrow(df)), Nipi = rep(999, nrow(df)), VarNipi = rep(999, nrow(df)))
  
  # Loop through each row to calculate CV per length class
  for (i in 1:nrow(df)) {
    # Extract variable values
    values <- as.numeric(df[i, variable_columns])
    values[is.na(values)] <- 0
    
    # Calculate ni and nipi
    if (length(variable_columns) > 1) {
      ni <- rowSums(df[i, variable_columns], na.rm = TRUE)
    } else {
      ni <- df[i, variable_columns]
    }
    nipi <- values
    
    # Calculate proportions and weighted values
    pi <- nipi / ni
    Ni <- df$total_lengths[i]
    Nipi <- Ni * pi
    Varpi <- (pi * (1 - pi)) / ni
    VarNipi <- (Ni^2) * Varpi
    
    # Update DFrame
    DFrame$Nipi[i] <- sum(Nipi, na.rm = TRUE)
    DFrame$VarNipi[i] <- sum(VarNipi, na.rm = TRUE)
    DFrame$CV[i] <- sqrt(sum(VarNipi, na.rm = TRUE)) / sum(Nipi, na.rm = TRUE) * 100
  }
  
  # Calculate total CV
  DFrame$total_CV <- ifelse(sum(DFrame$Nipi) == 0, 0, sqrt(sum(DFrame$VarNipi)) / sum(DFrame$Nipi) * 100)
  
  # Create a new data frame with Length, CV per length, and total CV
  result_df <- data.frame(Length = df$Length, CV_per_length = DFrame$CV)
  result_df$total_CV <- DFrame$total_CV[1]  # Total CV is the same for all rows
  
  return(result_df)
}