validateSampleStructure <- function(originalModesAntimodes, 
                                    newModesAntimodes, 
                                    epsilon, 
                                    theta_min, 
                                    theta_prop, 
                                    x_col, 
                                    original_data = NULL, 
                                    test_data = NULL, 
                                    bw = NULL,
                                    count_col = NULL) 
{

  if (is.null(newModesAntimodes) || is.null(originalModesAntimodes) || 
      nrow(newModesAntimodes) == 0 || nrow(originalModesAntimodes) == 0 ||
      !("type" %in% names(newModesAntimodes)) || !("type" %in% names(originalModesAntimodes))) {
    return(FALSE) 
  }
  

  orig_modes <- originalModesAntimodes[originalModesAntimodes$type == "mode", , drop = FALSE]
  new_modes  <- newModesAntimodes[newModesAntimodes$type == "mode", , drop = FALSE]
  
  orig_antimodes <- originalModesAntimodes[originalModesAntimodes$type == "antimode", , drop = FALSE]
  new_antimodes  <- newModesAntimodes[newModesAntimodes$type == "antimode", , drop = FALSE]
  
  ### check for number modes/antimodes
  if (nrow(orig_modes) != nrow(new_modes) || nrow(orig_antimodes) != nrow(new_antimodes)) {
    return(FALSE)
  }
  
  ### check for epsilon
  
  if (nrow(orig_modes) > 0) {
    orig_x <- as.numeric(orig_modes[[x_col]])
    new_x  <- as.numeric(new_modes[[x_col]])
    if (max(abs(orig_x - new_x), na.rm = TRUE) > epsilon*bw) 
      return(FALSE)
  }
  if (nrow(orig_antimodes) > 0) {
    orig_x_anti <- as.numeric(orig_antimodes[[x_col]])
    new_x_anti  <- as.numeric(new_antimodes[[x_col]])
    if (max(abs(orig_x_anti - new_x_anti), na.rm = TRUE) > epsilon*bw) 
      return(FALSE)
  }
  
  ### check for thetas
  if (nrow(orig_modes) > 0 && nrow(orig_antimodes) > 0) {

    orig_df <- data.frame(x = c(as.numeric(orig_modes[[x_col]]), 
                                as.numeric(orig_antimodes[[x_col]])), 
                          y = c(orig_modes$value, orig_antimodes$value))
    orig_df <- orig_df[order(orig_df$x), ]
    
    new_df  <- data.frame(x = c(as.numeric(new_modes[[x_col]]), 
                                as.numeric(new_antimodes[[x_col]])), 
                          y = c(new_modes$value, new_antimodes$value))
    new_df  <- new_df[order(new_df$x), ]
    
    orig_amps <- abs(diff(orig_df$y))
    new_amps  <- abs(diff(new_df$y))
    
    amp_ratios <- new_amps / orig_amps
    amp_ratios_clean <- amp_ratios[is.finite(amp_ratios)]
    if (length(amp_ratios_clean) == 0) 
      return(FALSE)
    

    if (min(amp_ratios_clean, na.rm = TRUE) < theta_min) 
      return(FALSE)
    if ((max(amp_ratios_clean, na.rm = TRUE) - min(amp_ratios_clean, na.rm = TRUE)) > theta_prop) 
      return(FALSE)
    
  } else if (nrow(orig_modes) > 0 && nrow(orig_antimodes) == 0) {
    
    peak_ratios <- new_modes$value / orig_modes$value
    if (min(peak_ratios, na.rm = TRUE) < theta_min) return(FALSE)
    

    if (length(peak_ratios) > 1) {
      if ((max(peak_ratios, na.rm = TRUE) - min(peak_ratios, na.rm = TRUE)) > theta_prop) 
        return(FALSE)
    }
    

    if (!is.null(original_data) && !is.null(test_data) && !is.null(count_col)) {
      orig_total <- sum(original_data[[count_col]], na.rm = TRUE)
      new_total  <- sum(test_data[[count_col]], na.rm = TRUE)
      
      orig_shares <- orig_modes$value / orig_total
      new_shares  <- new_modes$value / new_total
      
      if (max(abs(orig_shares - new_shares), na.rm = TRUE) > theta_prop) 
        return(FALSE)
      

      orig_nonzero <- sum(original_data[[count_col]] > 0, na.rm = TRUE)
      new_nonzero  <- sum(test_data[[count_col]] > 0, na.rm = TRUE)
      if (new_nonzero / orig_nonzero < theta_min) return(FALSE)
    }
  } else {
    
    return(FALSE)
  }
  
  return(TRUE)
}
