dataDistributionVisualization <- function(dataTableOriginal, 
                                          dataTableMinimal = NULL, 
                                          dataTableTested = NULL,
                                          color = NULL)
{
  
  ############### Labeller ###################
  
  smart_labeller <- function(labels) {
    res <- label_both(labels)
    lapply(res, function(x) stringr::str_wrap(x, width = 12))
  }
  
  
  
  
  plot_folder <- paste0(getwd(), "/figures")
  if (!dir.exists(plot_folder)) { dir.create(plot_folder) }
  
  var_name <- dataTableOriginal%>%select(matches("^(age|length)")) %>% names()
  var_n    <- dataTableOriginal%>%select(matches("^(n_)")) %>% names()
  
  if (grepl("age", var_name)) variable = "age" else variable = "length"
  
  if (var_n != "n_individuals") {
    dataTableOriginal$n_individuals <- dataTableOriginal[[var_n]]

    if (!is.null(dataTableMinimal)) dataTableMinimal$n_individuals <- dataTableMinimal[[var_n]]
    
    if (!is.null(dataTableTested)) dataTableTested$n_individuals  <- dataTableTested[[var_n]]

  } 

  technical_cols <- c("species", var_name, "n_individuals", "data", "bw")
  facet_vars <- setdiff(names(dataTableOriginal), technical_cols)
  if (length(facet_vars) == 0) facet_vars <- NULL
  
    current_bw <- as.numeric(unique(dataTableOriginal$bw))
    current_species <- unique(dataTableOriginal$species)
    if (is.null(color)) current_color <- "grey" else current_color <- color
    current_border_color <- current_color
    
    border_colors <- c("dark_frame" = "black")
    border_colors[current_species] <- current_color 
    
    dataTableOriginal$data <- "original"
    if (!is.null(dataTableMinimal)) dataTableMinimal$data  <- "reference"
    if (!is.null(dataTableTested)) dataTableTested$data  <- "tested"
    
    combined_data <- bind_rows(dataTableOriginal, dataTableMinimal, dataTableTested)
    combined_data$data <- factor(combined_data$data)
    combined_data$border_group <- ifelse(combined_data$data == "tested", "dark_frame", current_species)
    
    combined_data <- droplevels(combined_data)
    
    #if (grepl("age", var_name)) variable = "age" else variable = "length"
    
    x_label <- if(grepl("age", var_name)) {
      if(current_bw > 1) paste0("Age, year, binwidth ", current_bw) else "Age, year"
    } else {
      if(current_bw > 1) paste0("Length class, binwidth ", current_bw, " cm") else "Length class, cm"
    }
    
    title_suffix <- if(current_bw > 1) paste0(" (binwidth = ", current_bw, ")") else ""
    
    valid_x_values <- as.numeric(combined_data[[var_name]])
    if (length(valid_x_values) > 0 && any(!is.na(valid_x_values))) {
      max_var <- max(valid_x_values, na.rm = TRUE)
    } else {
      max_var <- 100
    }
    
    
    ################## function ##########
    
    group_vars <- c("data", facet_vars)
    
    label_data <- combined_data %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(total_n = sum(n_individuals, na.rm = TRUE), .groups = "drop") %>%
      mutate(data = case_when(data=="original" ~ "orig.",
                              data=="reference" ~ "ref.",
                              TRUE ~ "test."))%>%
      mutate(label_text = paste0(data, ": n = ", total_n)) %>%
      group_by(across(all_of(facet_vars))) %>% 
      summarise(
        label_block = paste(label_text, collapse = "\n"), 
        .groups = "drop"
      )
    
    
    ############ PLOT #############
    
    
    p <- ggplot(combined_data, aes(x = factor(.data[[var_name]]), y = n_individuals, 
                                           alpha = data)) +
      geom_bar(
              aes(color = ifelse(data == "tested", "dark_frame", species)),
               fill = current_color,
               stat = "identity", position = "identity", linewidth = 0.5) +
      scale_color_manual(values = border_colors, guide = "none") +
      #scale_color_manual(values = border_colors, guide = "none")  +
      scale_alpha_manual(values = c("original" = 0.15, "reference" = 0.90, "tested" = 0.40),
                         name = "data"
                            ) +
      labs(title = ifelse(is.null(dataTableMinimal) && is.null(dataTableTested),
                          paste0("Original data, ", variable, title_suffix),
                          ifelse(!is.null(dataTableMinimal) && is.null(dataTableTested),
                                 paste0("Original data vs reference distribution-preserving data, ", variable, title_suffix),
                                 paste0("Original data vs reference distribution-preserving data vs tested data, ", variable, title_suffix))),
           subtitle = current_species) +
      xlab(x_label) +
      ylab("Number sampled individuals") +

      theme(
        plot.title = element_text(face = "bold", size = 17),
        plot.subtitle = element_text(hjust = 0.5, size = 16, face = "bold", color = current_color),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(size = 13, face = "bold"),
        axis.title.y = element_text(size = 13, face = "bold"),
        legend.title = element_text(face="bold", size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.4, "cm"),
        strip.text = element_text(face="bold", size = 14),
        strip.clip = "off",
        panel.grid.minor = element_blank()
      ) +
      scale_x_discrete(breaks = as.character(seq(0, max_var, by = 5)), drop = FALSE) + 
      
      geom_text(
        data = label_data,
        aes(x = Inf, y = Inf, label = label_block), # Both x and y set to Inf for top right
        inherit.aes = FALSE,
        hjust = 1.1,                                # Pulls the text slightly inside from the right border
        vjust = 1.5,                                # Pulls the text slightly inside from the top border
        size = 4, 
        fontface = "bold.italic", 
        color = "black"
      )
    
    
    
    if (!is.null(facet_vars) && length(facet_vars)==1)
    {
      p <- p + facet_grid(rows = as.formula(paste(facet_vars, "~ .")), 
                          scales = "free_x", 
                          labeller = smart_labeller, drop = TRUE)
      
    }  else if (!is.null(facet_vars) && length(facet_vars) > 1)
    {
      if ("year" %in% facet_vars) {
        rows.plot <- "year"
        cols <- setdiff(facet_vars, "year")
      } else {
        rows.plot <- facet_vars[1]
        cols <- facet_vars[-1]
      }
      cols.plot <- if (length(cols) == 1) cols else paste(cols, collapse = " + ")
      p <- p + facet_grid(as.formula(paste(rows.plot, "~", cols.plot)), 
                          labeller = smart_labeller,
                          scales = "free_x", 
                          drop = TRUE)
      
      
    }
    
    
  filename <- ifelse(is.null(dataTableMinimal) && is.null(dataTableTested), 
               paste0("hist_plot_", variable, "_", current_species, "_original_bw_", current_bw, ".png"),
                ifelse(!is.null(dataTableMinimal) && is.null(dataTableTested),
                  paste0("hist_plot_", variable, "_", current_species, "_original_vs_reference_bw_", current_bw, ".png"),
                  paste0("hist_plot_", variable, "_", current_species, "_original_vs_reference_vs_tested_bw_", current_bw, ".png")))
       
        suppressWarnings({
           ggplot2::ggsave(file.path(plot_folder, filename), 
                         plot = p, width = 15, height = 10, dpi = 300)
       print(p)
       })

}



