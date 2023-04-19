#=================
# Demo_create_CA_from_RDBES
# Nuno Prista, SLU Aqua, Sweden, @WGBIOPTIM (from developments at WGRDBES-RAISE&TAF)
#=================

# ============================================
# Development notes

# fun_create_XX_from_RDBES need to be developed to other hierarchies
# all require pre-loading of create_RDBES_Id_table
# some require pre-loading of RDBES_codeLists
# ATT some specifications are SLU specific (may need generalization, change to RDB colnames, etc.)

# ============================================


rm(list=ls())

# -------------------------
# load packages
# -------------------------

library(remotes)
install_github("ices-tools-dev/RDBEScore")
library(RDBEScore)
library(data.table)

# -------------------------
# load functions
# -------------------------

source("000_Funs/fun_create_TR_from_RDBES.r")
source("000_Funs/fun_create_HH_from_RDBES.r")
source("000_Funs/fun_create_SL_from_RDBES.r")
source("000_Funs/fun_create_HL_from_RDBES.r")
source("000_Funs/fun_create_CA_from_RDBES.r")
source("000_Funs/createTableOfRDBESIds.r") # used in RDBEScore but not yet exported


# -------------------------
# Import and handle RDBES objects
# -------------------------


# Data import and loading [Inputs]
list.files("001_Inputs")
zipFiles <- c("001_Inputs/HCS_H5_2022_09_30_052520.zip")
RDBESdataObj <- importRDBESDownloadData(zipFiles)
validateRDBESDataObject(RDBESdataObj, verbose = FALSE)

myFields <- c("DEstratumName")
myValues <- c("CombGears_20_Cod_Q3")

# myFields <- c('FOarea')
# myValues <- c('27.3.a.21')

RDBESdataObj <- filterRDBESDataObject(RDBESdataObj, 
						fieldsToFilter = myFields,
                        valuesToFilter = myValues )

RDBESdataObj <- findAndKillOrphans(objectToCheck = RDBESdataObj, verbose = FALSE)

# -------------------------
# Generates RDB objects
# -------------------------

newTR<-create_TR_from_RDBES(x=RDBESdataObj, hierarchy = 5)
newHH<-create_HH_from_RDBES(x=RDBESdataObj, hierarchy = 5)
newSL<-create_SL_from_RDBES(x=RDBESdataObj, hierarchy = 5, dir_codelists = "001_Inputs/")
newHL<-create_HL_from_RDBES(x=RDBESdataObj, hierarchy = 5, dir_codelists = "001_Inputs/")
newCA<-create_CA_from_RDBES(x=RDBESdataObj, hierarchy = 5)
