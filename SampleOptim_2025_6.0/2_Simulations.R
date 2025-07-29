#######################################################################################################
#######################################################################################################
#######################################################################################################
##
##   Biological sampling optimization (Script "SampleOptim")
##   Developed by: Patricia Goncalves (patricia@ipma.pt)
##   Last version development period: may 2023
##   Version: v4.1
##   Updated: 18 april 2024
##
##   Reference:
##   Gonçalves, Patrícia 2021. "SampleOptim" a data analysis R-tool to optimize fish sampling for
##   biological parameters as input on fish stock assessment.
##
##   version:
#######################################################################################################
#######################################################################################################
#######################################################################################################

#'  changes:
#'  - paths (relative, more general)
#'  - remove creation of data_samplebio (previously done, we can load)
#'    - general question: how to organise scripts? Main?
#'  - params$a instead of separate params.a (keep it as option)
#'  - not done: should include quick range validation?


library(stringr)
library(tidyr) # for instance, to separate extra_otolith column

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
library(here)


# source scripts and define variables -----------------------
source("00_helper_functions.R") 
# Function for randomly select the samples by length class
source("sample_selection_function.R") 

set.seed(2019) #any value
SEP <- ","
VERBOSE <- FALSE
path1 <- paste0(output_dir)  # The common output directory

#############################################
############   READ INPUT SETTINGS FROM FILE ###################################

fname <- "input_params.csv" # replace T by Q to denote Quarter
sep_is <- getSep(fname) # or use data.table::fread


###### Input settings file
stt <- read.table(fname, sep = sep_is, header = FALSE, stringsAsFactors = FALSE) %>%
  t() %>%
  as.data.frame(stringsAsFactors = FALSE, row.names = FALSE) %>%
  setNames(.[1, ]) %>%
  slice(-1) # how generic this should be?


#======
#File Header; how general is this?
#Variable name;Mandatory;Variable.value;Definition
# in pdf documentation: Names.of.variables Mandatory Variable.options Default Definition
NAME = 0
MANDATORY = grep('mandatory', stt$`Variable name`, ignore.case = TRUE)
VALUE = grep('value', stt$`Variable name`, ignore.case = TRUE)
DEFINITION = grep('definition', stt$`Variable name`, ignore.case = TRUE)
#======

# Define logical TRUE/FALSE columns
col_logic <- c("age_only", "port", "distuniporto")

# Define numeric columns
col_num <- c("sex_ratio", "min_lc", "max_lc", "interval_lc", "min_age", "max_age",
             "min_otol", "max_otol", "interval_otol", "linf", "k", "t0",
             "year_start", "year_end", "stage_mature", "n")

# Extract parameters
param <- stt %>%
  slice(VALUE) %>% # filter by the variable value (or choose variable name = value)
  dplyr::rename_with( ~ str_remove(tolower(.), "^num_|^log_")) %>%
  dplyr::mutate(extra_otol = str_split(extra_otol, " ")) %>%  # Split into a list of numbers
  tidyr::unnest_wider(extra_otol, names_sep = "_") %>% # unnest into separate columns
  dplyr::mutate(across(c(all_of(col_num), matches("extra_otol")), ~ as.numeric(.)), across(all_of(col_logic), ~ as.logical(.))) %>%
  dplyr::rename(
    timeStrata = time_strata,
    # rename for consistency with previous code
    Linf = linf,
    K = k,
    numSubSamples = n
  )
#-------------------------


# previous param.species, now param$species
# if I wanted the same format I could do
#for (col in names(param)) {
# assign(paste0("param.", col), param[[col]], envir = .GlobalEnv)
#}


## Check the input names, types and values for the simulations
str(param)

# validate input ------------------

if ((param$min_lc > param$max_lc) |
    (param$min_age > param$max_age) |
    (param$min_otol > param$max_otol) |
    (param$year_start > param$year_end) |
    (param$numSubSamples <= 0) |
    (!str_detect(param$species, "^[A-Z]{3}$")) |# fao_code with 3 letters?
    (!param$timeStrata %in% c('Y', 'S', 'Q'))
    # logical codes are TRUE or FALSE
    # acceptable ranges
    # (param$min_lc < 0 | param$max_lc > 100) |  # Example range for min_lc and max_lc
    #  (param$min_age < 0 | param$max_age > 120) |  # Example range for min_age and max_age
    # (param$min_otol < 0 | param$max_otol > 50) # Example range for min_otol and max_otol
) { print("ERROR in input parameter table")}


# -------
# parameters set

# extra set of otoliths
extra <- param %>%
  dplyr::select(matches('extra_otol'))

otolitSet <- as.numeric(c(
  seq(param$min_otol, param$max_otol, by = param$interval_otol),
  unlist(extra)
))

#set length class set
class_length_set <- seq(param$min_lc, param$max_lc, param$interval_lc)

# von Bertalanffy growth model parameters used as starting values
Linf <- param$Linf # asymptotic size
K <- param$K # growth coefficient
t0 <- param$t0 # theoretical age when size is zero


timeInterval <- param$timeStrata #Time interval (Q="quarter"), for the otoliths selection by length class
sexRatio <- param$sex_ratio #0.5
numSubSamples <- param$numSubSamples # number of simulations (bootstrap runs)

#newage, define a step to create set of intermediate ages
newage <- seq(param$min_age, param$max_age, 0.1) ###set age distribution vector for predictions
newage <- data.frame(newage)
colnames(newage) <- "age"


year_init <- param$year_start
year_last <- param$year_end
yearSet <- c(year_init:year_last)



#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

#### 2.1 Data preparation to run the simulations:
###REMOVE NAs in age data (Note: I only use the individuals where age has been attributed)

data_samplebio_copy <- data_samplebio # create a copy in case I want to revert to
### Age classes (just removing date, why?)
data_samplebio <- data_samplebio |>
  dplyr::select(-date)


#-----------------================================
# this should be in the first or second script
# QUESTION: why this chunk of code?
# here we use data that was processed in 1. If it is to be standalone we may output it
# 1- If we remove NA in age and month, we can only get error
#    messages in either id_bio_fish (check in previous scripts) or length_class.
# 2- and if we do get error messages, what do we do about it?

# 3- is it just to check that they were removed? do it in script 1?

data_samplebio <- data_samplebio[!is.na(data_samplebio$age), ]
data_samplebio <- data_samplebio[!is.na(data_samplebio$month), ]


na_message_col <- c('age', 'month', 'length_class')

na_messages <- data_samplebio %>%
  summarise(across(all_of(na_message_col), ~ ifelse(
    any(is.na(.)), paste("ERROR->", cur_column(), "HAS NA\n"), ""
  ))) %>%
  unlist() %>%
  na.omit()

if (length(na_messages) > 0) {
  message(na_messages)
}


#--------------------
summary(data_samplebio)

# use only years that are both in data_samplebio and input_parameters
valid_years <- unique(data_samplebio$year)
anos <- intersect(yearSet, valid_years)

########################################################################################################
########################################################################################################
########################################################################################################
####2.2  Set function for von Bertalanffy growth model parameters (Linf, K and t0)
########################################################################################################

# Filter data for a specific year.
# Sample data based on specified parameters.
# Fits a non-linear model to the sampled data.
# Handles errors and warnings during the model fitting process.
# Returns a list of results for each iteration.

#' Process Fish Data for a Given Year
#'
#' This function processes fish data for a specified year, performing sampling and fitting a non-linear model.
#' It depends on the input parameters (defined above in list "param")
#' @param this_year An integer representing the year for which data is processed.
#' @param numOtolitsPerClass
#' @param porto
#' @param distUniPorto
#' @return A list containing the results of the sampling and model fitting, including coefficients and predicted values.

#' 
vonberByYear <- function(this_year, numOtolitsPerClass, porto, distUniPorto) {
  print(paste0("Vonber: ", this_year))
  
  # apply to each subSample and measure time to complete (param numSubSamples)
  system.time(vvv <- lapply(1:numSubSamples, function(i, this_year) {
    
    print(paste("V:", this_year, ":", numOtolitsPerClass, ":", i, sep = ""))
    
    # randomly selected "id_bio_fish" in each length class
    dataV <- temporalSample(
      data_samplebio[data_samplebio$year == this_year, ],
      numOtolitsPerClass, # number to sample
      # from param
      class_length_set, # classes for sampling, from param
      sexRatio, # desired ratio of females to males, in param
      timeInterval, #in param, time period for sampling (eg: "1S", "2S", "1T", "2T", "3T", "4T", "A", "S", "T")
      # choices
      porto, # sample by ports --->> decision to make up front!
      distUniPorto, # use uniform distribution for ports (default is TRUE)!!!! CHANGE
      verbose = VERBOSE # print information
    ) # dataV
    
    #--------------------------
    # Initialize variables for fitting the model
    fitTypical <- NULL # Placeholder for the fitted model
    skip <- FALSE  # only one <,  Flag to indicate if the fitting should be skipped
    
    #------
    # Attempt to fit a non-linear model
    tryCatch({
      fitTypical <- nls(vbTypical, data = dataV, start = svTypical, control)
    },
    error = function(e) { print(e); skip <<- TRUE}, ## in case it didn't converge, <<
    warning = function(w) { print(w); skip <<- TRUE},
    message = function(m) { print(m); skip <<- TRUE},
    finally = function(m) { print("will skip")}
    ) # tryCatch
    
    #-------
    # Check if fitting was successful or if it should be skipped
    if (skip  | is.null(fitTypical)) {
      return(list( year = this_year, data = dataV, coef = NA, Linf = NA, K = NA, t0 = NA, predict_values = NA, n = i, j = numOtolitsPerClass ) # kist
      ) # return
    } # end of if
    
    #-------------------------
    # Extract coefficients from the fitted model
    coef <- summary(fitTypical)$coefficients # Get coefficients
    Linf <- summary(fitTypical)$coefficients[[1]]  # Asymptotic length
    K <- summary(fitTypical)$coefficients[[2]]     # Growth coefficient
    t0 <- summary(fitTypical)$coefficients[[3]]    # Theoretical age at length zero
    
    # Predict values based on the fitted model
    predict_values <- predict(fitTypical, newage = newage) # generate predictions
    
    # return list
    return(list(year = this_year, data = dataV, coef = coef, Linf = Linf, K = K, t0 = t0,
                predict_values = predict_values, n = i, j = numOtolitsPerClass) )
    
  }, this_year))[[3]] #  Extract the third element from the list returned by lapply
  
  return(vvv)
}#



#########################################################################################################
#########################################################################################################
#########################################################################################################
####            **  START SIMULATIONS  **    ###########################################################
#########################################################################################################
### Note: repeat the code between lines 100 and 388 for each "numOtolitsPerClass" (number of otoliths by length class)
#########################################################################################################
#########################################################################################################
### Quarter/Sample
## SR=1:1
# numOtolitsPerClass=1  (Conditions:porto=FALSE, distUniPorto=FALSE)
# numOtolitsPerClass - is the number of otoliths selected by length class
#########################################################################################################
#########################################################################################################


for (numOtolitsPerClass in otolitSet) {
  
  #-------
  # Define common suffix for several plots and tables
  path2 <- paste0(numOtolitsPerClass, "_", timeInterval, ".png")  #
  # fig3_length_path <- file.path(path1, paste0("Fig3_length_", path2))
  path3 <- paste0(numOtolitsPerClass, ".png")
  
  #--------
  ##### Initial values set and fixed by species
  svTypical <- list(Linf = Linf, K = K, t0 = t0) ##Initial parameters values for the growth curve
  vbTypical <- length_class ~ Linf * (1 - exp(-K * (age - t0))) ##von Bertallanfy growth model
  
  control <- nls.control(maxiter = 10000)
  
  # Apply vonberByYear function to relevant years for given numOtolitsPerClass
  vonber <- sapply(anos, vonberByYear, 
                   numOtolitsPerClass = numOtolitsPerClass, 
                   porto = FALSE, distUniPorto = FALSE)
  
  ######################################
  ######################################
  #extracting the variables from each of the 100 samples (subsamples)
  #==================================
  #' Extract Variables from Dataset
  #'
  #' This function extracts specific variables from a given dataset and returns them as a data frame.
  #'
  #' @param x A list containing various attributes including year, data, Linf, K, t0, n, etc.
  #' @return A data frame containing the extracted variables: year, age, length_class, Linf, K, t0, n, sex, maturity_stage, individual_wg, and month.
  #' @examples
  #' # Example usage:
  #' data <- list(year = 2021, data = list(age = c(1, 2), length_class = c(10, 20), sex = c("M", "F"), maturity_stage = c(1, 2), individual_wg = c(5, 10), month = c(1, 2)), Linf = 50, K = 0.1, t0 = 0, n = 100)
  #' result <- extract_variables(data)
  #' print(result)
  #' 
  extract_variables <- function(x, id) {
    # Create a data frame to store the extracted variables
    data.frame(
      ID_sim = id,  # Include the index of the element
      year = x$year,  # Extract the 'year' variable from the input
      Linf = x$Linf,
      K = x$K,
      t0 = x$t0,
      n = x$n,
      age = x$data$age,  # Extract the 'age' variable from the nested 'data' list
      Lt = x$data$length_class,
      sex = x$data$sex,
      mat_stg = x$data$maturity_stage,
      wt = x$data$individual_wg,
      month = x$data$month
    )
  }#
  
  extracted_data <- do.call(rbind, lapply(seq_along(vonber), function(i) extract_variables(vonber[[i]], i)))

 
  #a<-unlist(vb_k)
  #b<-unique(extracted_data$K)
  #identical(a,b)
  #--------------------------------
  
  ###create matrix with data of the variables of each one of the 100 samples (sub-samples)
  
  dados_lt_age <- extracted_data %>%
    dplyr::select(Lt, age, year, month, ID_sim) %>%
    mutate(quarter = as.factor(quarter(month)), 
           type = numOtolitsPerClass)
  
  
  write.table(dados_lt_age,
              paste0(output_dir,"dados_lt_age_", numOtolitsPerClass,".csv"),
              sep = SEP, row.names = FALSE)
  
  #'---------------------------------------------
  ### Predict
  
  # Extract and bind data from the vonber list
  total <- do.call(rbind, lapply(vonber, function(x) {
    cbind(seq_along(x$predict_values), x$n, x$predict_values, x$j, x$year)
  }))
  
  colnames(total) <- c("ID_ind", "ID_sim", "pred_lt", "type", "year")
  
  vb_predict_melt <- total
  
  write.table(vb_predict_melt,
              paste0(output_dir, "vb_predict_melt_", numOtolitsPerClass, ".csv"),
              sep = SEP, row.names = FALSE)
  
  #-------------------------------------
  ## Bio data
  
  dados_bio <- extracted_data %>%
    dplyr::select(Lt, age, year, wt, month, ID_sim) %>%
    mutate(quarter = as.factor(quarter(month)),
           type = numOtolitsPerClass)
  
  write.table(dados_bio,
              paste0(output_dir, "dados_bio_", numOtolitsPerClass, ".csv"),
              sep = SEP, row.names = FALSE)
  
  
  #####
  if (param$age_only == FALSE) {
  
      ##Data of length, age, sex, maturity stage, weight, month and year
    dados_bio <- extracted_data %>%
      dplyr::select(Lt, age, year, sex, mat_stg, wt, month, ID_sim) %>%
      mutate(mat_stg = as.numeric(mat_stg),
             maturity = case_when(
               mat_stg < param$stage_mature | is.na(mat_stg) ~ 0,
               mat_stg >= param$stage_mature ~ 1,
               TRUE ~ NA_real_), 
             quarter = as.factor(quarter(month)), 
             type = numOtolitsPerClass)
    
    
    write.table(dados_bio,
                paste0(output_dir, "dados_bio_", numOtolitsPerClass, ".csv"), sep = SEP, row.names = FALSE)
    
    
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ###############################################################################
    ### MATURITY OGIVE
    ### DETERMINE: L25, L50, L75
    ##  Confidence intervals by year
    ##
    ## NOTE: SUBSET 1? QUARTER (SPAWNING SEASON)
    ##############################################################################
    ##############################################################################
    
    years <- unique(dados_bio$year)
    
    table_mature <- function(data = dados_bio) {
      results <- data.frame(year = numeric(), L25 = numeric(), L50 = numeric(), L75 = numeric(), ID_sim = integer(), type = character())
      
      for (nb in unique(data$ID_sim)) {
        # Fit a generalized linear model for maturity based on length 'Lt' for the current simulation ID
        glm1 <- glm(factor(maturity) ~ Lt, family = binomial, data = subset(data, ID_sim == nb))
        Lmat <- signif(dose.p(glm1, p = c(0.25, 0.50, 0.75)), digits = 3)
        
        results <- rbind(results, data.frame(
          year = unique(data$year[data$ID_sim == nb]),
          L25 = as.numeric(Lmat[1]),
          L50 = as.numeric(Lmat[2]),
          L75 = as.numeric(Lmat[3]),
          ID_sim = nb,
          type = unique(data$type)
        ))
      } # end loop over ID_sim
      
      return(results)
    } # table_mature
    
    #
    table_mo <- table_mature(data = dados_bio) ##Data from the whole year
    write.table(table_mo,
                paste0(output_dir, "table_res_mo_", numOtolitsPerClass, ".csv"), sep = SEP,row.names = FALSE)
    
  } #if params$age_only FALSE
  
  ####Figure 3 - compare length and age distributions by year (by simulations)
  ##Note: the length distribution did not change by simulations
  ##(because the selection is based in the number of otoliths by length class)
  ################ Data from simulations - figures (length and age distribution)
  
  Fig3_length <-
    file.path(paste0( output_dir, "Fig3_length_", numOtolitsPerClass,"_",timeInterval,".png"))
  
  png(file = Fig3_length)
  Fig3_length <- ggplot(dados_lt_age, aes(x = Lt, colour = factor(ID_sim))) +
    geom_density(show.legend = FALSE) + 
    facet_wrap( ~ year, ncol = 2) + 
    theme_classic()
  print(Fig3_length)
  dev.off()
  
  
  Fig3_age <-
    file.path(paste0(output_dir, "Fig3_age_", numOtolitsPerClass, "_", timeInterval, ".png"))
  
  png(file = Fig3_age)
  Fig3_age <- ggplot(dados_lt_age, aes(x = age, colour = factor(ID_sim))) +
    geom_density(show.legend = FALSE) + 
    facet_wrap( ~ year, ncol = 2) + 
    theme_classic()
  print(Fig3_age)
  dev.off()
  
  
  #==============================
  ####Determine mean length at age - original data and by simulation (for each year)
  
  # if I want the SD too
  
  # Creating a combined summary table
  table_original_combined <- data_samplebio %>%
    group_by(age, year) %>%
    summarize(
      m_lt = mean(length_class),
      sd_lt = sd(length_class),
      .groups = 'drop'
    ) %>%
    mutate(data = "original", 
           ID_sim = 0,
           type = numOtolitsPerClass)
  
  
  table_original <- table_original_combined %>% 
    dplyr::select(age, year, m_lt, data, ID_sim, type)
  
  

  #==============================
  table_simul_combined <- dados_lt_age %>%
    group_by(age, year, ID_sim) %>%
    summarize(
      m_lt = mean(Lt),
      sd_lt = sd(Lt),
      data = "simulations",
      type = numOtolitsPerClass) 
  
  
  table_simul <- table_simul_combined |>
    dplyr::select(age, year, m_lt, data, ID_sim, type) ###organize columns
  
  ######
  
  
  ## Combine original data and simulations in one table (data frame)
  table_original_simul <- merge(table_original, table_simul, all = TRUE)
  
  year_simul <- c(year_init:year_last)
  
  table_original_simulsub <- table_original_simul[table_original_simul$year %in% year_simul, ]
  
  write.table(table_original_simulsub,
              paste0(output_dir, "table_original_simulsub_mla_", numOtolitsPerClass,".csv"),
              sep = SEP, row.names = FALSE
  )
  
  #############
  ###Compare mean length at age from distributions of original data with the data from simulations
  Fig4_length <-
    file.path(paste0(output_dir, "Fig4_length_", numOtolitsPerClass, "_", timeInterval,".png"))
  png(file = Fig4_length)
  Fig4_length <- ggplot(table_original_simulsub,
                        aes(x = factor(age), # 
                            y = m_lt,
                            fill = factor(type))) +
    geom_bar(stat = "identity", position = position_dodge()) + 
    facet_wrap( ~ year, ncol = 2) + 
    xlab("Age") +
    ylab("Mean length (cm)") +
    theme_classic()
  print(Fig4_length)
  dev.off()
  
  
  Fig4_age <-
    file.path(paste0(output_dir, "Fig4_age_", numOtolitsPerClass,"_", timeInterval,".png"))
  png(file = Fig4_age)
  Fig4_age <- ggplot(table_original_simulsub,
                     aes(
                       x = factor(age),
                       y = m_lt,
                       colour = factor(type)
                     )) +
    geom_boxplot() +
    facet_wrap( ~ year, ncol = 2) + 
    xlab("Age") +
    ylab("Mean length (cm)") +
    theme_classic()
  print(Fig4_age)
  dev.off()
  
  
  #=======================
  # sd should be done at the same time as the rest
  #=======================
  ####Determine sd (length) at age - original data and by simulation (for each year)
  table_originalsd <- table_original_combined %>%
    dplyr::select(age, year, sd_lt, data, ID_sim, type)
  
  table_simulsd <- table_simul_combined %>%
    dplyr::select(age, year, sd_lt, data, ID_sim, type)  # Organize columns
  
  ######
  ## Combine data original and simulations in one table (data frame)
  table_original_simulsd <- merge(table_originalsd, table_simulsd, all =
                                    TRUE)
  
  table_original_simulsubsd <- table_original_simulsd[table_original_simulsd$year %in% year_simul, ]
  write.table(table_original_simulsubsd,
              paste0(output_dir,"table_original_simulsub_sd_", numOtolitsPerClass,".csv"),
              sep = SEP, row.names = FALSE)
  
  #############
  ###Compare standard deviation of length at age from distributions of original data with the data from simulations
  Fig5_length <-
    file.path(paste0(output_dir,"Fig5_length_",numOtolitsPerClass,"_",timeInterval,".png"))
  png(file = Fig5_length)
  Fig5_length <- ggplot(table_original_simulsubsd,
                        aes(x = factor(age), # length?
                            y = sd_lt,
                            fill = factor(type))) +
    geom_bar(stat = "identity", position = position_dodge()) + 
    facet_wrap( ~ year, ncol = 2) + 
    xlab("Age") +
    ylab("Standard deviation of length (cm)") +
    theme_classic()
  print(Fig5_length)
  dev.off()
  
  
  Fig5_age <-
    file.path(paste0(output_dir,"Fig5_age_",numOtolitsPerClass,"_", timeInterval,".png"))
  png(file = Fig5_age)
  Fig5_age <- ggplot(table_original_simulsubsd,
                     aes(x = factor(age),
                         y = sd_lt,
                         colour = factor(type))) +
    geom_boxplot() +
    facet_wrap( ~ year, ncol = 2) + 
    xlab("Age") +
    ylab("Standard deviation of length (cm)") +
    theme_classic()
  print(Fig5_age)
  dev.off()
  
  
  ################################################################
  #################################################################
  ### Growth parameters from the von Bertallanfy model by year
  
  trimestre_simulvb <- extracted_data %>% 
    dplyr::distinct(Linf, ID_sim, K, ID_sim, t0, ID_sim, year) %>%
    mutate(type = numOtolitsPerClass)
  
  
  write.table(
    trimestre_simulvb,
    paste0(output_dir,"results_simulvbgm_",numOtolitsPerClass,".csv"),
    sep = SEP,row.names = FALSE)
  
  
  ####For all the years, to compare the VGBGM parameters between years
  ####Figure 6 - Summary of parameters by year for the full set of simulations (n=100), by numOtolitsPerClass (number of selected otoliths)
  Fig6_K_VBGM <-
    file.path(paste0(output_dir,"Fig6_K_VBGM_",numOtolitsPerClass,"_", timeInterval,".png"))
  png(file = Fig6_K_VBGM)
  fig6_K <- ggplot(trimestre_simulvb, aes(x = factor(year), y = K)) +
    geom_boxplot() + 
    xlab("year") + 
    sample_theme() 
  print(fig6_K)
  dev.off()
  
  
  #'----
  Fig6_t0_VBGM <- file.path(
    paste0(output_dir, "Fig6_t0_VBGM_", numOtolitsPerClass, "_", timeInterval, ".png"))
  
  png(file = Fig6_t0_VBGM)
  fig6_t0 <- ggplot(trimestre_simulvb, aes(x = factor(year), y = t0)) +
    geom_boxplot() + 
    xlab("year") + 
    sample_theme() # in 00_helper_functions.R
  print(fig6_t0)
  dev.off()
  
  
  Fig6_Linf_VBGM <-
    file.path(paste0(output_dir,"Fig6_Linf_VBGM_",numOtolitsPerClass,"_",timeInterval,".png"))
  
  png(file = Fig6_Linf_VBGM)
  fig6_Linf <- ggplot(trimestre_simulvb, aes(x = factor(year), y = Linf)) +
    geom_boxplot() + 
    xlab("year") + 
    sample_theme() 
  print(fig6_Linf)
  dev.off()
  
  #### Von bertallanfy growth model
  svTypical <- list(Linf = Linf, K = K, t0 = t0) ##Initial growth parameters
  vbTypical <- Lt ~ Linf * (1 - exp(-K * (age - t0))) ##von Bertallanfy growth model
  control <- nls.control(maxiter = 10000)
  
  ###############################################################################
  ###############################################################################
  ### Determining mean square error between sim and obs (cost=mspe, mape, rtmspe)
  ### mspe(sim, obs, na.rm=TRUE)
  ##############################################################################
  ##############################################################################
  
  #=======================
  table_results <- function(data = dados_lt_age) {
    years <- unique(data$year)
    results <- matrix(NA, nrow = length(years), ncol = 8)
    colnames(results) <- c("year", "Linf", "k", "t0", "mspe", "mape", "rtmspe", "type")
    
    for (nb in seq_along(years)) {
      year_data <- data[data$year == years[nb], ]
      
      fitTypical <- tryCatch({
        nls(vbTypical, data = year_data, start = svTypical, control)
      }, error = function(e) {
        message("Error in fitting model for year ", years[nb], ": ", e$message)
        return(NULL)
      }, warning = function(w) {
        message("Warning in fitting model for year ", years[nb], ": ", w$message)
        return(NULL)
      })
      
      if (is.null(fitTypical)) {
        results[nb, ] <- c(as.numeric(years[nb]), rep(NA, 6), unique(data$type))
        next
      }
      
      # Fit with different costs: mspe, map, rtmspe
      metrics <- sapply(c(mspe, mape, rtmspe), function(cost) {
        cvFit(fitTypical, Lt ~ Linf * (1 - exp(-K * (age - t0))), 
              data = year_data, y = year_data$age, cost = cost, k = 10)$cv
      })
      
      results[nb, ] <- c(as.numeric(years[nb]), 
                         as.numeric(coef(fitTypical)), 
                         metrics, 
                         unique(data$type))
    }
    return(results)
  }
  
  
  table_res <- table_results(data = dados_lt_age)
  
  
  write.table( table_res,
               paste0(output_dir,"table_res_stat_",numOtolitsPerClass,".csv"),
               sep = SEP,row.names = FALSE)
  
  
} # for numOtolitsPerClass in otolitSet

######################################## END Simulations ##################################################
###########################################################################################################

