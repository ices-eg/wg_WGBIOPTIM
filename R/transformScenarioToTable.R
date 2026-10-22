transformScenarioToTable <- function(scenarioList) 
  {
  if (!is.null(scenarioList$variable)) scenarioList <- list(scenarioList)
  
  lapply(scenarioList, function(res) {
    base_row <- data.frame(
      variable              = res$variable,
      species               = res$species,
      N_original            = res$n_orig, 
      N_minimal             = res$n_min,
      N_tested              = res$n_test, 
      is_valid              = res$is_valid,
      inside_l1_corridor    = res$inside_l1_corridor,
      d_min_area            = res$d_min,
      corridor_max_area     = res$corridor_length,
      distance_to_tested    = res$distance_to_tested,
      total_penalties_score = res$total_penalties_score,
      penalty_topo          = res$penalty_topo,
      penalty_eps           = res$penalty_eps,
      penalty_theta         = res$penalty_theta,
      stringsAsFactors      = FALSE
    )
    if (length(res$facet_vars) > 0) {
      facet_vars_df <- as.data.frame(res$facet_vars, stringsAsFactors = FALSE)
      base_row <- bind_cols(facet_vars_df, base_row)
    }
    
    return(base_row)
  }) %>% bind_rows()
}
