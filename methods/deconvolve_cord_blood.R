#!/usr/bin/env R

# Author: Sean Maden
#
# Get cord blood cell type estimates using FlowSorted.Blood.EPIC::estimateCellCounts2()
#
#

# dependencies
libv <- c("FlowSorted.Blood.EPIC", "FlowSorted.CordBlood.450k", "HDF5Array")
sapply(libv, library, character.only = TRUE)


#----------
# load data
#----------
# get h5se rg
fname <- "remethdb_epic-hm850k_h5se_rg_1669220613_0-0-3"
fpath <- file.path("recount-methylation-files", "compilations", fname)
rg <- loadHDF5SummarizedExperiment(fpath)

#-----------------------------------------
# deconvolution, with cord blood reference
#-----------------------------------------
ref <- "CordBlood"
est <- estimateCellCounts2(rg[,c(1:10)], compositeCellType = ref)



