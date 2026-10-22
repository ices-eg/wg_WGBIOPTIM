dataAggregation <- function(dataList, 
                            variable = "length", 
                            aggregating_vars = NULL, 
                            species_vec = NULL, 
                            binwidth_vec = NULL)
  
{

  if (variable == "age")
  {
  
  var_name <- "age"
  variable_table <- subset(dataList$biologicalData, BVtypeMeas=="Age")%>%
    mutate(variable_defined = BVvalueMeas)
  
  } else
  {
  var_name <- "length_class_cm"
  variable_table <- dataList$lengthData%>%
    mutate(variable_defined = floor(FMclassMeas/10))
  }


variable_table <- variable_table%>%
  mutate(
    year = DEyear,
    month = as.numeric(substr(FOendDate, 6, 7)),
    quarter = case_when(is.na(month) ~ NA,
                        month %in% 1:3 ~ 1,
                        month %in% 4:6 ~ 2,
                        month %in% 7:9 ~ 3,
                        TRUE ~ 4),
    metier6 = FOmetier6,
    metier5 = FOmetier5,
    species = SAspeCodeFAO)

if (is.null(species_vec)) species_vec <- sort(unique(variable_table$species))


variable_table <- variable_table%>%
  rename(!!sym(var_name) := variable_defined)


if (is.null(aggregating_vars)) group_vars <- c("species", var_name) else 
  group_vars <- c("species", aggregating_vars, var_name)



########### BIN

if (is.null(binwidth_vec)) binwidth_vec <- c(rep(1, length(species_vec)))

if (length(species_vec) != length(binwidth_vec)) {
  if (length(binwidth_vec)==1)  binwidth_vec <- c(rep(binwidth_vec, length(species_vec))) else
    binwidth_vec <- c(rep(binwidth_vec[1], length(species_vec)))
}

bin_lookup <- setNames(binwidth_vec, species_vec)


#################### Aggregation #################

data_agg_list <- list()


for (i in 1:length(species_vec)) 
  {
  
  current_species <- species_vec[i]
  current_bw <- as.numeric(bin_lookup[current_species])
  
  variable_table_sub <- variable_table%>%dplyr::filter(species == current_species)
  
  if (nrow(variable_table_sub) == 0) next 
  
  if (current_bw > 1) 
    {
    new_col_name <- paste0(var_name, "_bin")
    group_vars_bin <- c(setdiff(group_vars, var_name), new_col_name)
    
    apply_binning <- function(data) 
      {
      if (is.null(data) || nrow(data) == 0) return(NULL)
      if (grepl("length", var_name)) {
        data <- data%>%mutate(!!new_col_name := floor(!!sym(var_name) / current_bw) * current_bw)
      } else {
        data <- data%>%mutate(!!new_col_name := !!sym(var_name))
      }
      data %>%
        summarise(n_individuals = sum(n, na.rm = TRUE), .by = any_of(group_vars_bin)) %>%
        arrange(across(all_of(group_vars_bin))) %>%
        rename(!!var_name := !!sym(new_col_name))
    }
    
    variable_table_aggregated <- apply_binning(variable_table_sub)%>%
      mutate(bw = current_bw)
  
    } else
    {
    variable_table_aggregated <- variable_table_sub%>%
    summarise(n_individuals = sum(n, na.rm=TRUE), .by = any_of(group_vars))%>%
    arrange(across(all_of(group_vars)))%>%
    mutate(bw = current_bw)
    }
      
  data_agg_list[[i]] <- variable_table_aggregated
  
  } 

return(dataAggregated = data_agg_list)

}
