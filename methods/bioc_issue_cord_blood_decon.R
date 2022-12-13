#----------
# load data
#----------
# dependencies
libv <- c("minfi", "FlowSorted.Blood.EPIC", "FlowSorted.CordBlood.450k", 
          "FlowSorted.CordBloodCombined.450k", "minfiData")
sapply(libv, library, character.only = T)
# load example data
rg <- get(data(RGsetEx))

#---------------------------------------
# 1. Use compositeCellType = "CordBlood"
#---------------------------------------
# attempt 1A: passing CordBlood to compositeCellType argument
compct <- "CordBlood"
estimateCellCounts2(rg, compositeCellType = compct)
# returns:
# > Error in estimateCellCounts2(rg, compositeCellType = compct) :
# > all elements of argument 'cellTypes' needs to be part of the reference 
# > phenoData columns 'CellType' (containg the following elements: '')

# attempt 1B: passing cellTypes from reference data
# get cell types from reference dataset
cb <- get(data(FlowSorted.CordBlood.450k))
ctv <- unique(cb$CellType)
estimateCellCounts2(rg, compositeCellType = compct, cellTypes = ctv)
# returns:
# > Error in p[trainingProbes, ] : subscript out of bounds

#-----------------------------------------------
# 2. Use compositeCellType = "CordBloodCombined"
#-----------------------------------------------
BiocManager::install("CordBloodCombined")

# attempt 2A: passing CordBlood to compositeCellType argument
compct <- "CordBloodCombined"
estimateCellCounts2(rg, compositeCellType = compct)

#-----------------------------------------------
# 3. Use compositeCellType = "CordBloodNorway"
#-----------------------------------------------
# attempt 3A: passing CordBlood to compositeCellType argument
compct <- "CordBloodNorway"
estimateCellCounts2(rg, compositeCellType = compct)
# Error in estimateCellCounts2(rg, compositeCellType = compct) :
#  all elements of argument 'cellTypes' needs to be part of the reference 
# phenoData columns 'CellType' (containg the following elements: '')

# attempte 3B: passing cellTypes from reference data
cb <- get(data(FlowSorted.CordBloodNorway.450k))
ctv <- unique(cb$CellType)
estimateCellCounts2(rg, compositeCellType = compct, cellTypes = ctv)

#------------------------------
# 4. Use CordTissueAndBlood
#------------------------------
BiocManager::install("FlowSorted.CordTissueAndBlood.450k")

compct <- "CordTissueAndBlood"
estimateCellCounts2(rg, compositeCellType = compct)
