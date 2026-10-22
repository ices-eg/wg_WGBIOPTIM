#### determining stable modes/antimodes
#### dataTable: data, already aggregated, 
#### e.g. by year, or year/quarter, or yer/quarter/metier, etc.
#### dataTable have to include: 
####  - column with pattern "length" or "age"
####  - column with pattern "n_" ("n_individuals" or something like this)


findLocalModesAnimodes <- function(dataTable, allowed_shift = NULL)
{
  
  if (is.null(dataTable$bw)) 
    binwidth <- min(diff(sort(unique(dataTable[[var_name]]))), na.rm = TRUE) else 
      binwidth <- unique(dataTable$bw)
  
  var_name <- dataTable%>%select(matches("^(age|length)")) %>% names()
  var_n    <- dataTable%>%select(matches("^(n_)")) %>% names()
  count <- var_n
  
  group_vars <- setdiff(names(dataTable), c(var_name, count, "bw"))
  
  dataTable <- dataTable[order(dataTable[[var_name]]), ]
  
  vector_orig <- as.numeric(dataTable[[count]])
  max_orig <- max(vector_orig, na.rm = TRUE)
  
## Allowed shift between modes determined from smoothing algorithm and real modes
  if (is.null(allowed_shift)) {
    class_step <- min(diff(sort(unique(dataTable[[var_name]]))), na.rm = TRUE)
    allowed_shift <- class_step * 2
  }
  
  # Filter Goley-Savitzky (p=3 - order of polynoms, n=5 - "moving window")
  if (nrow(dataTable) >= 5) {
    vector_smoothed <- as.numeric(signal::sgolayfilt(vector_orig, p = 3, n = 5))
    vector_smoothed <- pmax(vector_smoothed, 0) 
  } else {
    vector_smoothed <- vector_orig
  }
  max_smoothed <- max(vector_smoothed, na.rm = TRUE)
  
  ### Mode is considered as a mode if its hight > 5% from largest peak - should it be as parameter?
  extract_peaks <- function(v, is_orig, max_n, filter_neg) {
    p_modes <- findpeaks(v, nups = 1, ndowns = 1)
    if (filter_neg && !is.null(p_modes)) {
      p_modes <- p_modes[p_modes[, 1] > 0.05 * max_n, , drop = FALSE]
    }
    
    p_antimodes <- findpeaks(-v, nups = 1, ndowns = 1, 
                             zero = "+", minpeakheight = -Inf)
    
    res_list <- list()
    
    if (!is.null(p_modes) && nrow(p_modes) > 0) {
      valid_idx <- p_modes[, 2] <= nrow(dataTable) & p_modes[, 2] > 0
      if (any(valid_idx)) {
        sub_dataTable <- dataTable[p_modes[valid_idx, 2], , drop = FALSE]
        sub_dataTable$val_x <- sub_dataTable[[var_name]]
        sub_dataTable$val_y <- p_modes[valid_idx, 1]
        sub_dataTable$type  <- "mode"
        res_list[[length(res_list) + 1]] <- sub_dataTable
      }
    }
    
    if (!is.null(p_antimodes) && nrow(p_antimodes) > 0) {
      valid_idx <- p_antimodes[, 2] <= nrow(dataTable) & p_antimodes[, 2] > 0
      if (any(valid_idx)) {
        sub_dataTable <- dataTable[p_antimodes[valid_idx, 2], , drop = FALSE]
        sub_dataTable$val_x <- sub_dataTable[[var_name]]
        sub_dataTable$val_y <- -p_antimodes[valid_idx, 1]
        sub_dataTable$type  <- "antimode"
        res_list[[length(res_list) + 1]] <- sub_dataTable
      }
    }
    
    if (length(res_list) == 0) return(data.frame())
    res_block <- do.call(rbind, res_list)
    
    res_block <- res_block[!is.na(res_block$type) & !is.na(res_block$val_x), ]
    if (nrow(res_block) == 0) return(data.frame())
    
    res_block <- res_block[order(res_block$val_x), ]
    
    rle_info <- rle(res_block$type)
    res_block$type_run <- rep(seq_along(rle_info$lengths), rle_info$lengths)
    
    unique_runs <- unique(res_block$type_run)
    final_rows <- list()
    
    for (run in unique_runs) {
      run_data <- res_block[res_block$type_run == run, , drop = FALSE]
      if (run_data$type[1] == "antimode") {
        min_val <- min(run_data$val_y, na.rm = TRUE)
        min_data <- run_data[run_data$val_y == min_val, , drop = FALSE]
        center_idx <- ceiling(nrow(min_data) / 2)
        final_rows[[length(final_rows) + 1]] <- min_data[center_idx, ]
      } else {
        final_rows[[length(final_rows) + 1]] <- run_data[1, ]
      }
    }
    
    res_block <- do.call(rbind, final_rows)
    
    if (nrow(res_block) > 0) {
      while(nrow(res_block) > 0 && res_block$type[1] == "antimode") { 
        res_block <- res_block[-1, , drop = FALSE] 
      }
      while(nrow(res_block) > 0 && res_block$type[nrow(res_block)] == "antimode") { 
        res_block <- res_block[-nrow(res_block), , drop = FALSE] 
      }
    }
    return(res_block)
  }
  
  ### all original and smoothed modes together
  extrema_orig     <- extract_peaks(vector_orig,  
                                    is_orig = TRUE,  max_orig,
                                    filter_neg = TRUE)
  extrema_smoothed <- extract_peaks(vector_smoothed, 
                                    is_orig = FALSE, 
                                    max_smoothed, 
                                    filter_neg = FALSE)
  
  if (nrow(extrema_orig) == 0 || nrow(extrema_smoothed) == 0 || 
      is.null(extrema_orig$type) || is.null(extrema_smoothed$type)) { 
    return(data.frame()) 
  }
  
  # searching stable (robust) modes/antimodes
  stable_list <- list()
  for (j in 1:nrow(extrema_orig)) {
    orig_row <- extrema_orig[j, ]
    
    match_indices <- which(
      extrema_smoothed$type == orig_row$type &
        extrema_smoothed$val_x >= (orig_row$val_x - allowed_shift) &
        extrema_smoothed$val_x <= (orig_row$val_x + allowed_shift)
    )
    
    if (length(match_indices) > 0) {
      matches <- extrema_smoothed[match_indices, , drop = FALSE]
      matches$shift <- abs(matches$val_x - orig_row$val_x)
      best_match <- matches[which.min(matches$shift), ]
      
      out_row <- orig_row[, group_vars, drop = FALSE]
      out_row[[var_name]] <- orig_row$val_x
      out_row$value      <- orig_row$val_y 
      out_row$type       <- orig_row$type
      
      stable_list[[length(stable_list) + 1]] <- out_row
    }
  }
  
  if (length(stable_list) == 0) return(data.frame())
  
  result <- do.call(rbind, stable_list)
  result <- result[order(result[[var_name]]), ]
  
  # --- 
  while(nrow(result) > 0 && result$type[1] == "antimode") { 
    result <- result[-1, , drop = FALSE] 
  }
  while(nrow(result) > 0 && result$type[nrow(result)] == "antimode") { 
    result <- result[-nrow(result), , drop = FALSE] 
  }
  
  if (nrow(result) == 0) return(data.frame())
  
  clean_list <- list()
  clean_list[[1]] <- result[1, , drop = FALSE]
  
  if (nrow(result) > 1) {
    for (k in 2:nrow(result)) {
      last_added <- clean_list[[length(clean_list)]]
      current_row <- result[k, , drop = FALSE]
      
      if (current_row$type == last_added$type) {
        if (current_row$type == "mode") {
          if (current_row$value > last_added$value) {
            clean_list[[length(clean_list)]] <- current_row
          }
        } else {
          if (current_row$value < last_added$value) {
            clean_list[[length(clean_list)]] <- current_row
          }
        }
      } else {
        clean_list[[length(clean_list) + 1]] <- current_row
      }
    }
  }
  
  final_result <- do.call(rbind, clean_list)
  
  while(nrow(final_result) > 0 && final_result$type[1] == "antimode") { 
    final_result <- final_result[-1, , drop = FALSE] 
  }
  while(nrow(final_result) > 0 && final_result$type[nrow(final_result)] == "antimode") { 
    final_result <- final_result[-nrow(final_result), , drop = FALSE] 
  }
  
  return(final_result)
}
