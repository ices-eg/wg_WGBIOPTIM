
minimalReferenceSample <- function(data,
                                epsilon = 0, 
                                theta_min = 0.1, 
                                theta_prop = 0.5, 
                                max_steps = 20)
{
  var_name <- data%>%select(matches("^(age|length)"))%>%names()
  var_n <- data %>% select(matches("^(n_)")) %>% names()
  
  binwidth <- unique(data$bw)
  

  originalModesAntimodes <- findLocalModesAnimodes(data)
  

  if (nrow(originalModesAntimodes) == 0) {
    
    data$epsilon <- epsilon
    data$theta_min <- theta_min
    data$theta_prop <- theta_prop
    
    return(data)
  }
  
  original_counts <- data[[var_n]]
  best_counts <- original_counts
  
  # --- ЭТАП 1: БИНАРНЫЙ ПОИСК ---
  low <- 0.0
  high <- 1.0
  
  ### binary search algorithm
  
  for (step in 1:max_steps) 
    {
    
    mid <- (low + high) / 2
    

    target_N <- round(sum(original_counts, na.rm = TRUE) * mid)
    
    if (target_N < 5) {
      low <- mid
      
      data$epsilon <- epsilon
      data$theta_min <- theta_min
      data$theta_prop <- theta_prop
      
      next
    }
    
    p_orig <- original_counts / sum(original_counts, na.rm = TRUE)
    p_orig[is.na(p_orig)] <- 0
    
    test_data <- data
    
    B <- 1000
    
    set.seed(42) 
    
    ### monte carlo simulations
    monte_carlo_matrix <- rmultinom(n = B, size = target_N, prob = p_orig)
    
    ### average of monte carlo simulations
    mean_sampled_vector <- rowMeans(monte_carlo_matrix)
    
    final_realistic_vector <- round(mean_sampled_vector)
    
    
    
    diff_N <- target_N - sum(final_realistic_vector)
    if (diff_N != 0 && length(final_realistic_vector) > 0) {

      max_idx <- which.max(final_realistic_vector)
      final_realistic_vector[max_idx] <- final_realistic_vector[max_idx] + diff_N
    }
    
    test_data[[var_n]] <- final_realistic_vector
    
    newModesAntimodes <- findLocalModesAnimodes(test_data)
    
    is_valid <- validateSampleStructure(originalModesAntimodes = originalModesAntimodes, 
                                        newModesAntimodes      = newModesAntimodes, 
                                        epsilon                = epsilon, 
                                        theta_min              = theta_min, 
                                        theta_prop            = theta_prop, 
                                        x_col                  = var_name,
                                        original_data          = data,
                                        test_data              = test_data,
                                        bw                     = binwidth,
                                        count_col              = var_n)
    if (is_valid) {
      best_counts <- test_data[[var_n]]
      high <- mid  
    } else {
      low <- mid   
    }
  }
  
  
  data_changed <- data
  data_changed[[var_n]] <- best_counts
  
  data_changed$epsilon <- epsilon
  data_changed$theta_min <- theta_min
  data_changed$theta_prop <- theta_prop

  
  return(data_changed)
}
