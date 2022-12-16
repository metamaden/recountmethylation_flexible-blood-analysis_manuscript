libv <- c("HDF5Array", "minfi", "FlowSorted.CordBlood.450k")
sapply(libv, library, character.only = T)

rg.epic.fname <- "remethdb_h5se-rg_epic_0-0-2_1589820348"
rg.epic <- loadHDF5SummarizedExperiment(rg.epic.fname)

cb <- get(data("FlowSorted.CordBlood.450k"))
ctv <- unique(cb$CellType)


rgf <- RGChannelSet(Green = as.matrix(getGreen(rg.epic[,seq(5)])),
                    Red = as.matrix(getRed(rg.epic[,seq(5)])),
                    annotation = annotation(rg.epic))

estimateCellCounts(rgf, cellTypes = ctv, compositeCellType = "CordBlood")