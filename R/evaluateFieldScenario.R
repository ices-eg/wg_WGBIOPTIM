
evaluateScenario <- function(dataTableList, lambda = 0.25) {
  
  dataTableOriginal <- dataTableList$original
  dataTableMinimal  <- dataTableList$minimal
  dataTableTested   <- dataTableList$tested
  
  var_name <- dataTableOriginal %>% select(matches("^(age|length)")) %>% names()
  var_n    <- dataTableOriginal %>% select(matches("^(n_)")) %>% names()
  
  if (grepl("age", var_name)) variable = "age" else variable = "length"
  
  if (var_n != "n_individuals") {
    dataTableOriginal$n_individuals <- dataTableOriginal[[var_n]]
    dataTableMinimal$n_individuals  <- dataTableMinimal[[var_n]]
    dataTableTested$n_individuals   <- dataTableTested[[var_n]]
  }
  
  epsilon     <- unique(dataTableMinimal$epsilon)
  theta_min   <- unique(dataTableMinimal$theta_min)
  theta_prop <- unique(dataTableMinimal$theta_prop)
  
  technical_cols <- c("species", var_name, "n_individuals", "data", "bw")
  facet_vars <- setdiff(names(dataTableOriginal), technical_cols)
  if (length(facet_vars) == 0) facet_vars <- NULL
  
  current_bw      <- as.numeric(unique(dataTableOriginal$bw))
  current_species <- unique(dataTableOriginal$species)
  

  current_facets <- list()
  if (!is.null(facet_vars)) {
    current_facets <- lapply(facet_vars, function(v) unique(dataTableOriginal[[v]]))
    names(current_facets) <- facet_vars
  }
  
  # 1. Создаем мастер-ось X из оригинала
  classes <- dataTableOriginal %>% 
    select(all_of(c("species", var_name, setdiff(facet_vars, "species")))) %>% 
    distinct()
  

  align_and_get_counts <- function(target_df, classes, var_name, var_n) {
    aligned <- classes %>% left_join(target_df, by = intersect(names(classes), names(target_df)))
    aligned[[var_n]][is.na(aligned[[var_n]])] <- 0
    aligned <- aligned[order(aligned[[var_name]]), ]
    return(aligned)
  }
  
  df_orig_raw <- align_and_get_counts(dataTableOriginal, classes, var_name, "n_individuals")
  df_min_raw  <- align_and_get_counts(dataTableMinimal,  classes, var_name, "n_individuals")
  df_test_raw <- align_and_get_counts(dataTableTested,   classes, var_name, "n_individuals")
  
  x_coords <- as.numeric(df_orig_raw[[var_name]])
  if (length(x_coords) == 0 || any(is.na(x_coords))) {

    x_coords <- as.numeric(df_orig_raw[, 1]) 
  }
  

  sum_n_orig <- sum(df_orig_raw$n_individuals, na.rm = TRUE)
  sum_n_min  <- sum(df_min_raw$n_individuals, na.rm = TRUE)
  sum_n_test <- sum(df_test_raw$n_individuals, na.rm = TRUE)
  

  p_orig <- if (sum_n_orig > 0) df_orig_raw$n_individuals / sum_n_orig else 0
  p_min  <- if (sum_n_min > 0)  df_min_raw$n_individuals / sum_n_min else 0
  p_test <- if (sum_n_test > 0) df_test_raw$n_individuals / sum_n_test else 0
  

  cdf_orig <- cumsum(p_orig)
  cdf_min  <- cumsum(p_min)
  cdf_test <- cumsum(p_test)
  
  k_classes <- length(cdf_orig)
  dist_orig_to_min  <- sum(abs(cdf_orig - cdf_min)) / k_classes   
  dist_orig_to_test <- sum(abs(cdf_orig - cdf_test)) / k_classes  
  

  max_allowed_corridor <- dist_orig_to_min + lambda * sqrt(dist_orig_to_min)
  

  tolerance_floor      <- 0.05 
  #max_allowed_corridor <- max(max_allowed_corridor, tolerance_floor)
  
  is_inside_l1_corridor <- (dist_orig_to_test <= max_allowed_corridor)
  
  # --- 4. БИОЛОГИЧЕСКИЕ ШТРАФЫ (PENALTIES) ---
  df_orig_raw$p_pdf <- p_orig
  df_test_raw$p_pdf <- p_test
  orig_peaks <- findLocalModesAnimodes(df_orig_raw)
  test_peaks <- findLocalModesAnimodes(df_test_raw)
  
  penalty_topo  <- 0.0
  penalty_eps   <- 0.0
  penalty_theta <- 0.0
  
  if (nrow(orig_peaks) > 0 && nrow(test_peaks) > 0) {
    orig_modes <- orig_peaks[orig_peaks$type == "mode", , drop = FALSE]
    test_modes <- test_peaks[test_peaks$type == "mode", , drop = FALSE]
    orig_antis <- orig_peaks[orig_peaks$type == "antimode", , drop = FALSE]
    test_antis <- test_peaks[test_peaks$type == "antimode", , drop = FALSE]
    
### penalty 1
    delta_m <- abs(nrow(orig_modes) - nrow(test_modes))
    delta_a <- abs(nrow(orig_antis) - nrow(test_antis))
    if (delta_m > 0 || delta_a > 0) {
      penalty_topo <- (delta_m)^2 + (delta_a)^2
    }
    
### penalty 2
    max_allowed_shift <- epsilon * current_bw
    if (nrow(orig_modes) == nrow(test_modes) && nrow(orig_modes) > 0) {
      drift <- max(abs(as.numeric(orig_modes[[var_name]]) - as.numeric(test_modes[[var_name]])), na.rm = TRUE)
      if (drift > max_allowed_shift) penalty_eps <- drift - max_allowed_shift
    }
    
### penalty 3
    if (nrow(orig_modes) == nrow(test_modes) && nrow(orig_antis) == nrow(test_antis) && nrow(orig_antis) > 0) {
      orig_df <- data.frame(x = c(as.numeric(orig_modes[[var_name]]), as.numeric(orig_antis[[var_name]])), y = c(orig_modes$value, orig_antis$value)) %>% arrange(x)
      test_df <- data.frame(x = c(as.numeric(test_modes[[var_name]]), as.numeric(test_antis[[var_name]])), y = c(test_modes$value, test_antis$value)) %>% arrange(x)
      
      amp_ratios <- abs(diff(test_df$y)) / abs(diff(orig_df$y))
      amp_ratios_clean <- amp_ratios[is.finite(amp_ratios)]
      
      if (length(amp_ratios_clean) > 0) {
        if (min(amp_ratios_clean) < theta_min) penalty_theta <- penalty_theta + (theta_min - min(amp_ratios_clean))
        if ((max(amp_ratios_clean) - min(amp_ratios_clean)) > theta_prop) penalty_theta <- penalty_theta + ((max(amp_ratios_clean) - min(amp_ratios_clean)) - theta_prop)
      }
    }
  }
  
  total_penalties   <- penalty_topo + penalty_eps + penalty_theta
  is_scenario_valid <- is_inside_l1_corridor && (total_penalties == 0)
  
  return(list(
    variable = variable,
    species = current_species,
    facet_vars = current_facets,
    n_orig = sum_n_orig,                
    n_min  = sum_n_min,                           
    n_test = sum_n_test, 
    p_orig = p_orig,  
    p_min  = p_min,
    p_test = p_test,
    x_values = x_coords,
    is_valid = is_scenario_valid,
    inside_l1_corridor = is_inside_l1_corridor,
    d_min = dist_orig_to_min,                      
    corridor_length = max_allowed_corridor,        
    distance_to_tested = dist_orig_to_test,        
    total_penalties_score = total_penalties,
    penalty_topo = penalty_topo,
    penalty_eps = penalty_eps,
    penalty_theta = penalty_theta
  ))
}
