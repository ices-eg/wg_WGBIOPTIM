#######################################################################################################
#######################################################################################################
#######################################################################################################
##
##
##   Biological sampling optimization (Script "SampleOptim")
##   Developed by: Patricia Goncalves (patricia@ipma.pt)
##   Last version development period: June 2021
##   Version: v3.1
##
##   Reference:
##   Gonçalves, Patrícia 2019. "SampleOptim" a data analysis R-tool to optimize fish sampling for
##   biological parameters as input on fish stock assessment.
##
##
##
#######################################################################################################
#######################################################################################################
#######################################################################################################

#--------------------------------------------------
#' In this script:
#'
#' load data
#' rename variables
#' missing values
#' basic descriptive statistics
#' plot: outliers, correlations

#' changes:
#' - working directory to location of this file
#' - rename variables: english + consistency with pdf document
#' - minor change in path definition to make it more universal
#'
#' - pattern of missing variables
#' - other plots
#' - would be good: targeted pre-diagnostics to see if later assumptions hold
#'
#' QUESTION: if this is the first script, should it act as 'main' with
#'           selection of all the options upfront (replacement or not, etc)

#--------------------------------------------------


# Set working directory to source file location
library("rstudioapi")
setwd(dirname(getActiveDocumentContext()$path))

library(janitor) # data cleaning and standardisation of names
library(colorblindcheck) # accessibility eg palette_check(rainbow_pal, plot = TRUE)
library(visdat) # to visualise patterns of missing data
library(naniar) # patterns of missing data
library(RColorBrewer)
library(here)

##Packages:
library(FSA)
library(FSAdata)
library(nlstools)
library(reshape)
library(ggplot2)
library(ggthemes)
library(cvTools)
library(dplyr)
library(robustbase)
library(MASS)
library(psyphy)
library(boot)
library(RCurl)



########################################################################################################
#### Files path and function source:
output_dir <- "output" # remove '/' and later use file.path (to work in various operating systems)
#create directory if it doesn't exist
if (!dir.exists(output_dir)) {dir.create(output_dir)}

# helper functions, eg themes for plots
source('00_helper_functions.R')

#########################################################################################################

#### 1. Data Preliminary analysis: (Exploratory analysis)

#### missing data: patterns and percentage


vis_miss(data_samplebio)

#count total missing values in each column
sapply(data_samplebio, function(x) sum(is.na(x)))

#number of missings in each column, broken down by a categorical variable
gg_miss_fct(x = data_samplebio, fct = port_name)
gg_miss_fct(x = data_samplebio, fct = length_class)

data_samplebio %>%
  group_by(port_name) %>%
  miss_var_summary()

gg_miss_var(data_samplebio, show_pct = TRUE)

miss_var_summary(data_samplebio)





#### Summary

summary(data_samplebio)
str(data_samplebio)

table(data_samplebio$length_class, data_samplebio$month) ##summary of the number of individuals by length class and month

plot(data_samplebio$length_class, data_samplebio$month,
     xlab = "Length class", ylab="Month of sampling")

# alternative to check for outliers


# use plot_box function in 0_helper_functions.R
plot_box(data_samplebio, 'month', "length_class")

plot_box(data_samplebio, 'port_name', "individual_wg")

plot_box(data_samplebio, 'month', "age")

#----------------


### Number of samples by Port, year and month
nsamples_year_month <- data_samplebio |>
  group_by(port_name, month, year) |>
  summarise(nb_samples = n()) |>
  ungroup()
#write.table(nsamples_year_month, "numbersamples_summary_WHB.csv",sep=SEP)

n2 <- data_samplebio |>
  group_by(port_name, month, year) |>
  count()



# bar plot, include stack option
create_custom_plot(nsamples_year_month, "month", "nb_samples", "port_name", theme_Publication())



#### Length classes of the samples by Port, year and month
lengthclass_samples_year_mes1 <- data_samplebio %>%
  group_by(port_name, month, year) %>%
  count(length_class)
#write.table(lengthclass_samples_year_mes, "numberlengthclasses_samples_summary_WHB.csv",sep=SEP)

lengthclass_samples_year_mes <- data_samplebio |>
  group_by(port_name, year, month, length_class) |>
  summarise(nb_length_class = n()) |>
  ungroup()



########Figure a - length distribution samples by Port by year and month
Port <- unique(data_samplebio$port_name) ## list of Ports names
year <- sort(unique(lengthclass_samples_year_mes$year)) ##list of years on the samples data

for(bb in 1:length(year)){

  # convert the month column from number to abbreviation (1 for Jan, 2 for Feb, etc.)
  lengthclass_samples_year_mes$month <- factor(lengthclass_samples_year_mes$month,
                                               levels = 1:12,
                                               labels = month.abb)


  plota <- ggplot(lengthclass_samples_year_mes[lengthclass_samples_year_mes$year==year[bb],],
                  aes(x=length_class, y=nb_length_class, colour = port_name)) + 
    xlab("length") + 
    ylab("number of individuals")+
    geom_line()  +
    facet_wrap(~factor(month)) +
    sample_theme()

  dev.copy(png, paste0(output_dir, year[bb],"_length_distribution_samples_Port_year",".png"))
  print(plota)
  dev.off()
}





#### Age of the samples by Port, year and month

age_samples_year_mes <- data_samplebio |>
  group_by(port_name, year, month, age) |>
  summarise(nb_age = n()) |>
  ungroup() |>
  mutate(month = factor(month, levels = 1:12, labels = month.abb))


########Figure b - age distribution samples by Port by year and month
#Port<-unique(data_samplebio$port_name) ## list of Ports names
year <- sort(unique(age_samples_year_mes$year)) ##list of years on the samples data

for(bb in 1:length(year)){
  
  age_data <- age_samples_year_mes[age_samples_year_mes$year == year[bb],] |>
    mutate(age = as.numeric(as.character(age)))
  
  plotb <- age_data |>
    ggplot(aes(x = age, y = nb_age, colour = port_name)) +
    geom_line() +
    xlab("age") +
    ylab("number of individuals") +
    facet_wrap(~factor(month)) +
    custom_theme
  
  # Specify width and height for the PNG
  dev.copy(png, filename = paste0(output_dir, year[bb], "_age_distribution_samples_Port_year", ".png"), width = 800, height = 600)
  

  # dev.copy(png, paste(output_dir, year[bb], "_age_distribution_samples_Port_year", ".png", sep = ""))
   print(plotb)
   dev.off()
}

# Convert age to a factor
age_samples_year_mes$age <- as.factor(age_samples_year_mes$age)

m <- age_samples_year_mes$age%>%unique()%>%length()
custom_colors <- brewer.pal(n = m, name = "Paired")

# Create a bar plot
ggplot(age_samples_year_mes, aes(x = month, y = nb_age, fill = age)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Age Distribution by Month", x = "Month", y = "Count of Age Samples") +
  theme_Publication() +
  scale_fill_manual(values = custom_colors) +  # Apply the custom color palette
  facet_wrap(~ port_name)


#  alternative line plot
age_samples_year_mes |> 
  mutate(age = as.numeric(as.character(age))) |>
  ggplot(aes(x = month, y = nb_age, color = as.factor(age), group = age)) +
  geom_line(size = 1) +  # Adjust line thickness as needed
  labs(title = "Age Distribution by Month", x = "Month", y = "Count of Age Samples") +
  theme_Publication() +
  scale_color_manual(values = custom_colors) +  # Apply the custom color palette
  facet_wrap(~ port_name)

#create_custom_plot(age_samples_year_mes, "month", "n", "age", theme_classic())


## Figure 1 - Length distribution by year
years <- year[year %in% unique(data_samplebio$year)]
#lines_plot<-round(length(years)/2)
#par(mfrow=c(1,1))
for(nb in 1: length(years)){
  fig1 <- hist(data_samplebio$length_class[data_samplebio$year == years[nb]],
               xlab = "length",
               ylab = "number of individuals",
               main = years[nb])
  dev.copy(png, paste(output_dir, years[nb],
                      "_length_distribution",
                      ".png", 
                      sep = ""))
  dev.off()
}


###Figure 2 - Age distribution by year

for(nb in 1: length(years)){
  fig2 <- hist(data_samplebio$age[data_samplebio$year == years[nb]],
              xlab = "age",
              ylab = "number of individuals",
              main = years[nb])
  dev.copy(png, paste(output_dir, years[nb], "_age_distribution", ".png", sep=""))
  dev.off()
}

ggplot(data_samplebio[data_samplebio$year == years[nb], ], aes(x = factor(year), y = age)) +
  geom_violin(fill = "lightblue") +
  labs(x = "Year", y = "Age", title = paste("Violin Plot for Year", years[nb])) +
  theme_minimal()



#########################################################################################################
#########################################################################################################
#########################################################################################################


# Sample data creation
set.seed(123)
years <- rep(2000:2005, each = 100)
values <- rnorm(600, mean = rep(1:6, each = 100), sd = 0.5)
data <- data.frame(year = years, value = values)

# Density plot with facets
ggplot(data, aes(x = value, fill = as.factor(year))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ year) +
  labs(title = "Density Plots of Values Over Years", x = "Value", y = "Density") +
  theme_minimal()


# Combined density plot
ggplot(data, aes(x = value, color = as.factor(year), fill = as.factor(year))) +
  geom_density(alpha = 0.3) +
  labs(title = "Combined Density Plot of Values Over Years", x = "Value", y = "Density") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()



#------------------------------------
# categorical variables. For instance, is sex undetermined because fish is not mature?
ggplot(data = data_samplebio) +
  geom_count(mapping = aes(x = maturity_stage, y = sex))

ggplot(data = data_samplebio) +
  geom_count(mapping = aes(x = maturity_stage, y = port_name))


# continuous variables
ggplot(data = data_samplebio) +
  geom_point(mapping = aes(x = age, y = length_class), alpha = 1 / 100)


# bin age to acts as a categorical variable.
ggplot(data = data_samplebio, mapping = aes(x = age, y = length_class)) +
  geom_boxplot(mapping = aes(group = cut_width(age, 1)))



#https://rdrr.io/rforge/geo/man/ices.html
library(corrplot)

# Remove rows with NA values and filter out constant columns
data_cleaned <- data_samplebio %>%
  dplyr::select(where(is.numeric)) %>%
  na.omit() %>%
  select_if(~ sd(.) != 0)  # Remove constant columns

# Calculate the correlation matrix
corr <- cor(data_cleaned)


corrplot(cor(data_cleaned), method = "ellipse",
         title = "method = 'ellipse'")

