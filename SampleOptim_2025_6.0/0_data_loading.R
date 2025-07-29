#############################################################
## This is the data loading script that is part of:
## Biological sampling optimization (R tool "SampleOptim")
##
## Tasks:
## - loads data
## - checks that the data contains the right information/columns
## - renames columns for consistency with subsequent scripts
## - recodes missing variables (check format?)
#############################################################

library(dplyr)
#--------------------

###Biological sample data (Applied to a period of years)
data <- read.table("SampleOptim_input_fileformat_example.csv", sep = ";",
                   header = TRUE)


# create copy to rename/transform
data_samplebio <- data 

# rename columns  ----------------------------------------

# check that we have the necessary columns
glimpse(data_samplebio) 


#------------------
# fish identifier: if the relevant column names is "ID_BIO_FISH"
# check that there is an unique and different ID_BIO_FISH for each row
check_id <- length(unique(data_samplebio$ID_BIO_FISH)) == nrow(data_samplebio)


# if there is no previous fish identifier, create one assuming that each row...
if(!check_id){
data_samplebio <- data_samplebio |>
  mutate(ID_BIO_FISH = row_number())
}#
#-------------------



# Define the name vector for renaming
# name_vector <- c(old_name1 = "new_name1", ...)
name_vector <- c(ID_BIO_FISH = "id_bio_fish",
                 date = "date", 
                 Month = "month", 
                 Year = "year",
                 Species = "fao_code",
                 Port = "port_name", 
                 LengthClass = "length_class",
                 Weight = "individual_wg",
                 Sex = "sex",
                 MaturityStage = "maturity_stage", 
                 Age = "age")


# Identify columns to rename
to_rename <- names(data_samplebio) %in% names(name_vector)

# Renaming columns ----
names(data_samplebio)[to_rename] <- name_vector[names(data_samplebio)[to_rename]]


# format columns --------------------

data_samplebio <- data_samplebio |>
  janitor::clean_names() |> # standardise names (already in renaming)
  dplyr::mutate(
    across(c(id_bio_fish, year, length_class, age), as.integer),
    month = ifelse(month %in% 1:12, as.integer(month), NA),
    individual_wg = as.numeric(individual_wg),
    date = as.character(date), # could use custom function or "lubridate" to guess date format, but no need
    sex = factor(sex, levels = c("F", "M", "I")),
    port_name = factor(port_name), # could use external file with valid ports
    maturity_stage = as.character(maturity_stage) # !!!
  )


#---------

data_samplebio <- data_samplebio |>
  dplyr::mutate(across(where(is.character), ~ na_if(.x, ""))) # replace empty values by NA?
