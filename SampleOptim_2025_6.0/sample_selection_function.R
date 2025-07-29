#' This script contains the simulation data operating model
#'It defines the following functions
#' - bySex: Selects a sample of fish based on sex ratio and total numbers specified.
#' - byClassbySex: select individuals by class length and sex (if sexR = TRUE)
#' - byPortbyClassbySex: select individuals by port, length class and sex (sexR=T) 
#' - temporalSample
#'-----------------------------------------------------

## bySex

#' Selects a sample of fish based on sex ratio and total numbers specified.
#'
#' @param parte A data frame containing fish data with a column for sex.
#' @param sexR A boolean indicating whether to choose based on sex ratio.
#' @param numTotF Total number of females to select.
#' @param numTotM Total number of males to select.
#' @return A vector of selected fish IDs.

bySex <- function(parte, sexR, numTotF, numTotM){

  # Check if the input data frame has the required columns
  required_columns <- c("sex", "id_bio_fish")
  if (!all(required_columns %in% colnames(parte))) {
    stop("The input data frame 'parte' must contain the columns: sex, id_bio_fish.")
  }

  #If empty data frame, output a message and no further will execute
  if (nrow(parte) == 0) {  # Check if the input data frame has any rows.
    message("The input data frame is empty. Returning NULL.")
    return(NULL)
  }

  #------------------------------------------------------
  # the following will execute if the data frame is not empty
  q_s1 <- NULL;  # Initialize the output variable to store selected fish IDs.
  if(sexR){ # if select by sexRatio

    #------
    # choose females
    numF <- numTotF; # Set the number of females to select.
    parteF <- parte[parte$sex == "F", ]; # Filter the data frame for females.
    if (dim(parteF)[1] > numF) {  # Check if there are enough females.
      q_s1_F <- sample(parteF$id_bio_fish, numF, replace = FALSE)  # Sample females.
    } else {
      q_s1_F <- parteF$id_bio_fish;  # Select all females if not enough.
    }
    numF <- length(q_s1_F);  # Update the count of selected females.

    # -------
    # Choose males.
    numM <- numTotM;  # Set the number of males to select.
    parteM <- parte[parte$sex == "M", ];  # Filter the data frame for males.
    if (dim(parteM)[1] > numM) {  # Check if there are enough males.
      q_s1_M <- sample(parteM$id_bio_fish, numM, replace = FALSE)  # Sample males.
    } else {
      q_s1_M <- parteM$id_bio_fish;  # Select all males if not enough.
    }
    numM <- length(q_s1_M);  # Update the count of selected males.

    q_s1 <- c(q_s1_F, q_s1_M);  # Combine selected females and males.

    # -------
    # Choose indeterminate to complete the specified number of quantity.
    numI <- numTotF + numTotM - numF - numM;  # Calculate the number of indeterminate fish needed.
    if (numI > 0) {  # Check if indeterminate fish are needed.
      parteI <- parte[parte$sex == "I", ];  # Filter the data frame for indeterminate fish.
      if (dim(parteI)[1] > numI) {  # Check if there are enough indeterminate fish.
        q_s1_I <- sample(parteI$id_bio_fish, numI, replace = FALSE)  # Sample indeterminate fish.
      } else {
        q_s1_I <- parteI$id_bio_fish;  # Select all indeterminate fish if not enough.
      }
      q_s1 <- c(q_s1, q_s1_I);  # Combine selected indeterminate fish.
    }
  #----------------------------------------------------
  } else {  # Choose without specifying sex ratio.
    quantidade <- numTotF + numTotM;  # Calculate total quantity needed.
    if (nrow(parte) > quantidade) {  # Check if there are enough fish.
      q_s1 <- sample(parte$id_bio_fish, quantidade, replace = FALSE);  # Sample fish without sex consideration.
      if (quantidade != length(q_s1)) {  # Check for sampling errors.
        print(paste("ERROR in quantity:", quantidade, "different from:", length(q_s1)));  # Print error message.
      }
    } else {
      q_s1 <- parte$id_bio_fish;  # Select all fish if not enough.
    }
  } # choose without specifying sex ratio.


  return(q_s1); # Return the vector of selected fish IDs.
} # bySex


#----------------------------------------------
#' Sample fish based on sex ratio or not 
#'
#' @param parte A data frame containing fish data with columns 'sex' and 'id_bio_fish'.
#' @param sexR A boolean indicating whether to sample based on sex ratio, in temporalSample.
#' @param numTotF The total number of female fish to sample.
#' @param numTotM The total number of male fish to sample.
#' @param replace_choice A logical value indicating whether to sample with replacement. Default is FALSE.
#'
#' @return A vector of sampled fish identifiers, id_bio_fish.
#'
bySex_T <- function(parte, sexR, numTotF, numTotM, replace_choice = FALSE) {

  #---------------------
  # input validation:
  # If empty data frame, output a message and return NULL
  if (nrow(parte) == 0) {
    message("The input data frame is empty. Returning NULL.")
    return(NULL)
  }
  
  # Check if the input data frame has the required columns
  required_columns <- c("sex", "id_bio_fish")
  if (!all(required_columns %in% colnames(parte))) {
    stop("The input data frame 'parte' must contain the columns: sex, id_bio_fish.")
  }

  #--------------------
  # Key helper function to sample fish based on required count

  sampleFish <- function(data, count, replace_choice = replace_choice) {
    if (nrow(data) <= count) {
      return(data$id_bio_fish)
    }
    return(sample(data$id_bio_fish, count, replace = replace_choice))
  }
  
  #--------------------

  if (sexR) { # if sampling should consider sex ratio
    sampled_F <- sampleFish(parte[parte$sex == "F", ], numTotF) # sample female
    sampled_M <- sampleFish(parte[parte$sex == "M", ], numTotM) # sample male
    sampled_all <- c(sampled_F, sampled_M)  # Combine sampled female and male fish IDs
    # if after sampling we still didn't reach the desired number
    numI <- numTotF + numTotM - length(sampled_F) - length(sampled_M)
    if (numI > 0) { # if there are fish with indeterminate sex, sample them
      sampled_all <- c(sampled_all, sampleFish(parte[parte$sex == "I", ], numI))
    }
  } else { # Choose without specifying sex ratio
    total_count <- numTotF + numTotM
    sampled_all <- sampleFish(parte, total_count)
  }

  #-------------
  # Return the vector of selected fish IDs.
  return(sampled_all)
}


#-----------------------------------------------

##
# byClassbySex
#
# select individuals by class length and sex (if sexR = TRUE)
#
byClassbySex <- function(everybody, classes, sexR, numTotF, numTotM, verbose = FALSE){

  idRes <- NULL;
  for(cl in classes){   # for each length class
    parte <- everybody[everybody$length_class == cl, ];
    q_s1 <- bySex(parte, sexR, numTotF, numTotM) # sample by sex within the class
    idRes <- c(idRes, q_s1)
  }

  if(verbose) print(paste("Number of individuals: ", length(idRes)))

  return(idRes);
}


byClassbySex_T <- function(everybody, classes, sexR, numTotF, numTotM, verbose = FALSE) {
  
  idRes <- unlist(lapply(classes, function(cl) {
    parte <- everybody[everybody$length_class == cl, ]
    bySex(parte, sexR, numTotF, numTotM)  # sample by sex within the class
  }))
  
  return(idRes)
}

#----------------------------------------------

##
# byPortbyClassbySex
#
# select individuals by port, length class and sex (sexR=T) (numTotF - total number of females, numTotM - total number of males)
# A given (numTotF+numTotM) number of individuals by port of each length class (distUniPorto=FALSE) (maximum value = numeroPortos*numeroClasses*quantidade)

# distUniPorto = TRUE:  individuals for each length are selected  (maximum value = numeroClasses*quantidade)
#  with an uniform distribution between ports
#
# todos - includes individuals selected according to the selected temporal pattern
# distUniPorto - determines how "porto" affects the selected samples
# conjClasses - set of classes
#

#' Distributes quantities across ports based on specified classes and sex ratio.
#'
#' @param todos A data frame containing port data.
#' @param classes A vector of classes to categorize the distribution.
#' @param sexR Logical, defined in "amostraTemporal", indicating if to sample by sex
#' @param sexRatio A numeric value representing the ratio of females to males.
#' @param quantidade Total quantity to be distributed across ports.
#' @param distUniPorto Logical indicating if the distribution should be uniform.
#' @param verbose Logical indicating if detailed output should be printed.
#'
#' @return A vector of distributed quantities for each port.
#'
byPortbyClassbySex <- function(todos, classes, sexR, sexRatio, quantidade, distUniPorto, verbose = FALSE){

  idRes <- NULL; # Initialize an empty vector to store results
  conjPortos <- unique(todos$port_name); # unique port names
  numPortos <- length(conjPortos) # number of unique ports

  if(distUniPorto){ # if uniform distribution across ports is required

    # quantity to be assigned to each port based on uniform distribution
    quantidadePorPorto <- trunc(dunif(seq(1, 1, length = numPortos)) / numPortos * quantidade)

    # remaining quantity to be distributed
    falta <- quantidade - sum(quantidadePorPorto)
    # vector to randomly assign the remaining quantity to ports
    resto <- c(rep(1, falta), rep(0, (numPortos - falta)))

    # Update the quantity assigned to each port with the random distribution
    quantidadePorPorto <- quantidadePorPorto + sample(resto, numPortos)
  } else {
    # If not uniformly distributing, assign the total quantity to each port
    quantidadePorPorto <- array(quantidade, numPortos)
  }


  # Number of females and males based on the sex ratio
  numF <- ceiling(quantidadePorPorto*sexRatio)
  numM <- quantidadePorPorto - numF

  nporto <- 1 # Initialize counter for the port index
  for(p in conjPortos){ # Loop through each port to process the distribution
    if(verbose) print(paste("Port: ", p))
    pPorto <- todos[todos$port_name == p, ] # filter data by current port
    # Call the function to distribute quantities by classes and sex
    q_s1 <- byClassbySex(pPorto, classes, sexR,numF[nporto], 
                              numM[nporto], verbose = verbose)
    idRes <- c(idRes, q_s1) # Append the results to the idRes vector
    nporto <- nporto + 1
  }

  return(idRes)
}



############################ LAST VERSION ###############################

#' Sample Temporal Data
#'
#' This function samples temporal data based on specified parameters such as quantity, sex ratio, and time period.
#'
#' @param tab A data frame containing the data to be sampled.
#' @param quantidade The total number of individuals to select by set of lenght classes
#' @param conjClasses A vector of length classes to be considered for sampling.
#' @param sexRatio numeric value indicating desired proportion of females to males (ex: < 0: selection doesn't depend on sex, 1 - 100% females, 0 - 100% males, 0.2 - 20% females 80% males)
#' @param tm A character string indicating the time period for sampling (options: "1S", "2S", "1T", "2T", "3T", "4T", "A", "S", "T").
#' @param porto A logical value indicating whether to sample by ports (default is TRUE).
#' @param distUniPorto A logical value indicating whether to use distance from ports (default is TRUE).
#' @param verbose A logical value indicating whether to print additional information during execution (default is FALSE).
#'
#' @return # array with randomly selected "id_bio_fish" for the specified parameters in each length class


temporalSample <- function(tab, quantidade, conjClasses, sexRatio = NaN, tm = "Q", porto=TRUE, distUniPorto=TRUE, verbose=F){

  # Add semester and quarter columns to the input data frame
  tab$semester <- semester(tab$month);
  tab$quarter <- quarter(tab$month);
  

  original <- tab;
  idRes <- NULL

  # Filter the data based on the specified time period (tm)
  tab <-switch(tm,
              "1S" = tab[tab$semester==1,],   # first semester
              "2S" = tab[tab$semester==2,],   # second semester
              "1Q" = tab[tab$quarter==1,],  # first quarter
              "2Q" = tab[tab$quarter==2,],  # second quarter
              "3Q" = tab[tab$quarter==3,],  # third quarter
              "4Q" = tab[tab$quarter==4,],  # fourth quarter
              "Y"  = tab,                     # year
              "S"  = tab,                     # semester
              "Q"  = tab                      # quarter
  );

  
  # Check if the filtered data is empty and print an error message if so
  if(is.null(tab)){
    print(paste("ERROR: Time period was not correctly entered: <", tm, ">"));
    return (NULL);
  }


  if(verbose) print(paste("Time period: ", tm));
  if(verbose) print(paste("Number of individuals: ",dim(tab)[1]));


  sexR <- TRUE;  # Initialize sex ratio flag
  numTotF <- quantidade;  # Total number of females to sample
  numTotM <- 0;  # Total number of males to sample

  if(!is.nan(sexRatio) & (sexRatio >= 0) & (sexRatio <= 1)){
    # Calculate the number of females and males to sample based on sex ratio
    numTotF<-ceiling(quantidade * sexRatio);
    numTotM<-quantidade -numTotF;
    if(verbose) print(paste("SexRatio: ", sexRatio, "Females:", numTotF, "Males:", numTotM));
  } else {
    sexR <- FALSE;  # Set sex ratio flag to FALSE if not applicable
    if(verbose) print("SexRatio: Ignore");
  }

  if(verbose) print(paste("Number of classes: ", length(conjClasses)));

  # Sampling logic for semester
  if(tm == "S"){

    ano <- unique(tab$year);
    for(a in ano){  # Loop through each year
      for(s in 1:2){  # Loop through each semester
        if(verbose) print(paste("year:", a, " ", s, "S", sep=''))
        parte<-tab[tab$year == a & tab$semester == s, ]
        if(porto)
          q_s1 <- byPortbyClassbySex(parte, conjClasses, sexR, sexRatio, quantidade, distUniPorto, verbose=verbose)
        else
          q_s1 <- byClassbySex(parte, conjClasses, sexR, numTotF, numTotM, verbose = verbose)
        idRes <- c(idRes, q_s1);  # Append results to the ID list
      }
    }# loop over year
  }else if(tm == "Q"){
  # Sampling logic for quarter
    ano <- unique(tab$year);
    for(a in ano){
      for(t in 1:4){
        if(verbose) print(paste("year:", a, " ", t, "Q", sep=''))
        parte<-tab[tab$year == a & tab$quarter == t,]
        if(porto){
          q_s1 <- byPortbyClassbySex(parte,conjClasses,sexR,sexRatio,quantidade,distUniPorto, verbose=verbose)
        }else{
          q_s1 <- byClassbySex(parte, conjClasses,sexR,numTotF,numTotM, verbose=verbose)
          }
        idRes <- c(idRes, q_s1)
      }
    }#loop over year
  }else if(tm == "Y"){
  # Sampling logic for year
    ano <- unique(tab$year);
    for(a in ano){
      if(verbose) print(paste("year:",a,sep=''))
      parte <- tab[tab$year == a,]
      if(verbose) print(dim(parte))
      if(porto){
        q_s1 <- byPortbyClassbySex(parte, conjClasses, sexR, sexRatio, quantidade, distUniPorto, verbose=verbose)
      }else{
        q_s1 <- byClassbySex(parte, conjClasses, sexR, numTotF, numTotM, verbose=verbose)
      }  
      if(verbose) print(length(q_s1))
      idRes <- c(idRes, q_s1);
    }
  }else{ #tab includes choice among tm= "1S","2S","1T","2T","3T" or "4T"
    ano<-unique(tab$year);
    for(a in ano){
      if(verbose) print(paste("ano:",a," ",tm,sep=''))
      parte<-tab[tab$year==a,]
      if(verbose) print(dim(parte))
      if(porto){
        q_s1 <- byPortbyClassbySex(parte,conjClasses,sexR,sexRatio,quantidade,distUniPorto, verbose=verbose)
      }else{
        q_s1 <- byClassbySex(parte, conjClasses,sexR,numTotF,numTotM, verbose=verbose)
      }  
      idRes<-c(idRes, q_s1);
    }
  }
  
  return(original[original$id_bio_fish %in% idRes,])
}

#------------------------------------
#' Sample Temporal Data
#'
#' Samples temporal data based on specified parameters including time period, sex ratio, and other options
#' calls custom functions "porPortosPorClassesPorSexo" or "porClassesPorSexo"
#'
#' @param tab A data frame containing the data to be sampled.
#' @param quantidade The total number of samples to be taken.
#' @param conjClasses A set of classes to be used for processing.
#' @param sexRatio A numeric value representing the ratio of females to males.
#' @param tm A character string indicating the time period (e.g., "Q" for quarters, 1S for first semester).
#' @param porto A logical value indicating whether to include port data.
#' @param distUniPorto A logical value indicating whether to use uniform distribution for ports.
#' @param verbose A logical value indicating whether to print additional information.
#'
#' @return A data frame containing the sampled data.
#'
temporalSample_T <- function(tab, quantidade, conjClasses, sexRatio = NaN, tm = "Q", porto = TRUE, distUniPorto = TRUE, verbose = FALSE) {

  # Add semester and quarter columns to the data frame (just call them Q and S)
  tab$semester <- semester(tab$month)
  tab$quarter <- quarter(tab$month)

  original <- tab  # Store the original data frame for final output

  # Define the time period to loop over (semester, quarter, or year)
  loop_column <- if (grepl("S", tm)) "semester" else if (grepl("Q", tm)) "quarter" else NULL
  
  #=====================================================
  #' here I would only filter tab if there are digits, otherwise chose all
  # How many time periods to loop over 
  period_count <- ifelse(grepl("\\d", tm), 
                         # If digits are found, extract and convert them to numeric
                         as.numeric(gsub("[^0-9]", "", tm)),
                         # if no digits found, then 2 if semesters, 4 if quarter, 1 if year
                          switch(loop_column, semester = 2, quarter = 4, 1))

  #---------------------
  # Filter the data based on the specified time period (tm)
  if (!is.null(loop_column)) {tab <- tab[tab[[loop_column]] %in% (1:period_count), ]}
  
  #=====================================================
  #---------------------
  #' why this?
  if (nrow(tab) == 0) {
    print(paste("ERROR: Time period was not correctly entered: <", tm, ">"))
    return(NULL)
  }

  if (verbose) print(paste("Time period: ", tm, " | Number of individuals: ", nrow(tab)))
  #---------------------
  # Calculate the number of females and males based on sex ratio
  if(!is.nan(sexRatio) & (sexRatio >= 0) & (sexRatio <= 1)){
    sexR <- TRUE
    # Calculate the number of females and males to sample based on sex ratio
    numTotF <- ceiling(quantidade * sexRatio)
    numTotM <- quantidade - numTotF
    if(verbose) print(paste("SexRatio: ", sexRatio, "Females:", numTotF, "Males:", numTotM));
  } else {
    sexR <- FALSE
    numTotF <- quantidade  # Total number of females to sample
    numTotM <- 0  # Total number of males to sample
    if(verbose) print("SexRatio: Ignore")
  }

  #-----------------------
  # helper function to process data based on year, period and porto -----------
  process_data <- function(year, period) {
    parte <- tab[tab$year == year & period, ]
    if (porto) {
      return(byPortbyClassbySex(parte, conjClasses, sexR, sexRatio, quantidade, distUniPorto, verbose = verbose))
    } else {
      return(byClassbySex(parte, conjClasses, sexR, numTotF, numTotM, verbose = verbose))
    }
  }

  # loop over time period ----------------------------
  all_years <- unique(tab$year)
  idRes <- NULL
  for (a in all_years) {
    if (verbose) print(paste("year:", a))
    if (!is.null(loop_column)) {
      for (p in 1:period_count) {
        if (verbose) print(paste(" ", p, ifelse(loop_column == "semester", "S", "Q"), sep = ''))
        idRes <- c(idRes, process_data(a, tab[[loop_column]] == p))
      }
    } else {
      idRes <- c(idRes, process_data(a, TRUE))
    }
  }#end loop over time period

  #------------------
  # return
  return(original[original$id_bio_fish %in% idRes, ])
}
#--------------------------
