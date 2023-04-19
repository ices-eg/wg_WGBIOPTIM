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

# Import and tweak dummy RDBES data
	# note: at present the conversion tool assumes lower hierarchy C (all individuals biologically sampled)
list.files("001_Inputs")
RDBESdataObj <- createRDBESDataObject(rdbesExtractPath = "001_Inputs/h1_v_1_19_13") # h1 dummy data 
RDBESdataObj <- createRDBESDataObject(rdbesExtractPath = "001_Inputs/h5_v_1_19_13") # h5 dummy data 
if(unique(RDBESdataObj$BV$BVtypeMeas)=="Age") RDBESdataObj$BV$BVtypeMeas<-"LengthTotal"

validateRDBESDataObject(RDBESdataObj, verbose = FALSE)

RDBESdataObj <- findAndKillOrphans(objectToCheck = RDBESdataObj, verbose = FALSE)

# -------------------------
# Generates RDB objects
# -------------------------

# run for hierarchy 1 example
newTR<-create_TR_from_RDBES(x=RDBESdataObj, hierarchy = 1)
newHH<-create_HH_from_RDBES(x=RDBESdataObj, hierarchy = 1)
newSL<-create_SL_from_RDBES(x=RDBESdataObj, hierarchy = 1, dir_codelists = "001_Inputs/")
newHL<-create_HL_from_RDBES(x=RDBESdataObj, hierarchy = 1, dir_codelists = "001_Inputs/")
newCA<-create_CA_from_RDBES(x=RDBESdataObj, hierarchy = 1)

# run for hierarchy 5 example
newTR<-create_TR_from_RDBES(x=RDBESdataObj, hierarchy = 5)
newHH<-create_HH_from_RDBES(x=RDBESdataObj, hierarchy = 5)
newSL<-create_SL_from_RDBES(x=RDBESdataObj, hierarchy = 5, dir_codelists = "001_Inputs/")
newHL<-create_HL_from_RDBES(x=RDBESdataObj, hierarchy = 5, dir_codelists = "001_Inputs/")
newCA<-create_CA_from_RDBES(x=RDBESdataObj, hierarchy = 5)
