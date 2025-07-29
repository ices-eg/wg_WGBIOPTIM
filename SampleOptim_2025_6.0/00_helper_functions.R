library(RColorBrewer)
library(colorblindcheck)
library(ggplot2)


colors <- brewer.pal(7, "RdYlBu") #color blind friendly palette

# Define your theme
custom_theme <- theme_classic() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = 10))


#'----
#'fread
#' Get the separator used in a file
#'
#' This helper function reads the first line of a specified file and determines
#' whether the separator is a comma or a semicolon based on the number of
#' elements split by each character.
#'
#' @param fname A character string representing the file name or path.
#' @return A character string indicating the separator used in the file.
getSep <- function(fname) {
  # Read the first line of the file specified by 'fname'
  linha <- readLines(con = fname, n = 1L, ok = TRUE, warn = TRUE, encoding = "unknown",
    skipNul = TRUE
  )

  # Split the line by comma and store the result in 'a'
  a <- unlist(strsplit(linha, ","))

  # Split the line by semicolon and store the result in 'b'
  b <- unlist(strsplit(linha, ";"))

  # Compare the lengths of 'b' and 'a' to determine the separator
  if (length(b) > length(a)) {
    return(";")  # Return semicolon if it has more elements
  } else {
    return(",")  # Return comma otherwise
  }
}

#' https://www.rdocumentation.org/packages/data.table/versions/1.17.0/topics/fread
#' The separator between columns. Defaults to the character in the set [,\t |;:] that separates the sample of rows into the most number of lines with the same number of fields. Use NULL or "" to specify no separator; i.e. each line a single character column like base::readLines does.



#======================================================

# https://analysisfunction.civilservice.gov.uk/policy-store/data-visualisation-colours-in-charts/#section-9
discrete_palette <- c("#12436D", "#28A197", "#801650", "#F46A25", "#3D3D3D",
                      "#A285D1")

blue_theme <- function() {
  theme(
    # add border 1)
    panel.border = element_rect(colour = "blue", fill = NA, linetype = 2),
    # color background 2)
    panel.background = element_rect(fill = "aliceblue"),
    # modify grid 3)
    panel.grid.major.x = element_line(colour = "steelblue", linetype = 3, size = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y =  element_line(colour = "steelblue", linetype = 3, size = 0.5),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "steelblue", face = "italic", family = "Times New Roman"),
    axis.title = element_text(colour = "steelblue", family = "Times New Roman"),
    axis.ticks = element_line(colour = "steelblue"),
    # legend at the bottom 6)
    legend.position = "bottom"
  )
}

#https://best-practice-and-impact.github.io/afcharts/articles/cookbook.html#grouped-bar-chart
# Define a custom theme function
sample_theme <- function() {
  ggplot2::theme_classic() +
    theme(
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.line = element_line(size = 0.5),
      axis.text = element_text(size = 10)
    )
}


#=================================================

plot_box <- function(data = data_samplebio, x_col, y_col = "length_class") {
  if (!(x_col %in% names(data))) {
    stop(paste("Error: Column", x_col, "does not exist in the dataset."))
  }
  if (!(y_col %in% names(data))) {
    stop(paste("Error: Column", y_col, "does not exist in the dataset."))
  }

  if (!is.numeric(data[[y_col]])) {
    stop(paste("Error: Column", y_col, "is not numeric."))
  }

  boxplot(data[[y_col]] ~ data[[x_col]],
          xlab = paste(x_col),
          ylab = paste(y_col),
          main = paste0(y_col, " by ", x_col),
          col = "lightblue")  # Rotate x-axis labels)
}

#================================================




#------------------------------------------------
#set overall theme
theme_Publication <- function(base_size=10, base_family="Arial") {

  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.text = element_text(),
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.spacing  = unit(0, "cm"),
           legend.title = element_text(face="italic"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0")
   ))

}








standard.heatmap = function(mat, numbers = TRUE, digits = 2, cex = 0.8, 
                            abs = FALSE, main = "") {
  
  col.l = colorRampPalette(c('ivory', 'tomato'))(30)
  
  my.padding = list(
    layout.heights = list(
      top.padding = 0,
      main.key.padding = 0,
      key.axis.padding = 0,
      axis.xlab.padding = 0,
      xlab.key.padding = 0,
      key.sub.padding = 0),
    layout.widths = list(
      left.padding = 0,
      key.ylab.padding = 0,
      ylab.axis.padding = 0,
      axis.key.padding = 0,
      right.padding = 0),
    panel.background = list(
      col = "grey90"
    )
  )
  
  if (abs)
    transformed = abs(mat)
  else
    transformed = mat
  
  levelplot(transformed, col.regions = col.l, colorkey = FALSE,
            scales = list(x = list(rot = 90)), xlab = "", ylab = "", main = main,
            panel = function(y, x, z, ...) {
              
              panel.levelplot(y = y, x = x, z = z, ...)
              
              if (numbers) {
                
                if (abs)
                  rounded = as.character(round(mat[cbind(x, y)], digits))
                else
                  rounded = as.character(round(z, digits))
                rounded[is.na(rounded)] = ""
                ltext(x = x, y = y, labels = rounded, cex = cex)
                
              }#THEN
              
            }, par.settings = my.padding
  )
  
}#STANDARD.HEATMAP

standard.histogram = function(var, legend) {
  
  histogram(var,
            scales = list(tck = c(1, 0)),
            xlab = legend, ylab = "percent of total", col = "skyblue",
            panel = function(...) {
              
              panel.grid(h = -1, v = 0)
              panel.histogram(...)
              
            })
  
}#STANDARD.HISTOGRAM

standard.xyplot = function(formula, xlab, ylab, regression = FALSE) {
  
  xyplot(formula, data = data, pch = 19, 
         scales = list(tck = c(1, 0)),
         xlab = xlab, ylab = ylab,
         panel = function(...) {
           
           panel.grid(h = -1, v = -1)
           panel.xyplot(..., col = "skyblue")
           if (regression)
             panel.smoother(..., col = "skyblue3", lwd = 3, span = 1/2, level = 0.95)
           
         })
  
}#STANDARD.XYPLOT

prepare.map = function(data, fun) {
  
  latitudes = sort(unique(data[, "lat"]))
  longitudes = sort(unique(data[, "lon"]), decreasing = FALSE)
  map = matrix(0, nrow = length(longitudes), ncol = length(latitudes), 
               dimnames = list(longitude = longitudes, latitude = latitudes))
  
  for (i in seq_along(longitudes))
    for (j in seq_along(latitudes)) {
      
      location = subset(data, (lon == longitudes[i]) & (lat == latitudes[j]))
      
      if (nrow(location) == 0)
        map[i, j] = NA
      else
        map[i, j] = fun(location)
      
    }#FOR
  
  map[is.nan(map)] = NA
  
  return(map)
  
}#PREPARE.MAP




#'---------------------------------
#'
#'library(ggplot2)

create_custom_plot <- function(data, x_col, y_col, fill_col, plot_theme = theme_minimal()) {
  ggplot(data, aes_string(x = x_col, y = y_col, fill = fill_col)) +
    geom_bar(stat = "identity", position = position_dodge2(width = 1, preserve = "single")) +
    labs(title = "Monthly Sample Counts by Port",
         x = "Month",
         y = "Number of Samples") +
    #scale_x_continuous(breaks = 1:12) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +  # Ensure x-axis has consistent breaks and labels
    #facet_grid(year ~ .) +
    facet_wrap(~ year, ncol = 1, labeller = label_both) +  # Arranges years vertically with labels
    plot_theme +
    theme(legend.position = "bottom") +
    scale_fill_viridis_d()
}

#'-------------------------------


#----------------------------------

#' semester
#' Determine the semester based on the month.
#' This function takes a numeric month (1-12) and returns the corresponding semester.
#'
#' @param month A numeric value representing the month (1-12).
#' @return An integer indicating the semester (1 or 2).
#' @examples
#' semester(1)  # Returns 1
#' semester(7)  # Returns 2
semester <- function(month) {
  
  # error handling: check if 'month' is not numeric or is outside the range of 1 to 12
  if (!is.numeric(month) || any(month < 1)==TRUE || any(month > 12)==TRUE) {
    stop("Input must be numeric values between 1 and 12.")
  }
  return(ceiling(month / 6))
}

#--------------------------------------

#' quarter
#' Determine the quarter based on the month.
#' This function takes a numeric month (1-12) and returns the corresponding quarter.
#'
#' @param month A numeric value representing the month (1-12).
#' @return An integer indicating the semester (1, 2, 3 or 4).
#' @examples
#' quarter(1)  # Returns 1
#' quarter(5)  # Returns 2

quarter <- function(month){
  
  # error handling: check if 'month' is not numeric or is outside the range of 1 to 12
  if (!is.numeric(month) || any(month < 1)==TRUE || any(month > 12)==TRUE) {
    stop("Input must be numeric values between 1 and 12.")
  }
  
  return(ceiling(month/3))
}
#=======================================
