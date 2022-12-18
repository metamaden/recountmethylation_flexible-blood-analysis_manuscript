#!/usr/bin/env R

# Author: Sean Maden
#
# Use an updated deconvolution function to estimate cord blood cell types, 
# including enucleated red blood cells.
#

libv <- c("minfi", "minfiData", "FlowSorted.CordBlood.450k", "HDF5Array")
sapply(libv, library, character.only = TRUE)
source("estimateCellCounts.R")

#----------
# load data
#----------
# rg sets
rg.450k.fname <- "remethdb_h5se-rg_hm450k_0-0-2_1607018051"
rg.epic.fname <- "remethdb_h5se-rg_epic_0-0-2_1589820348"
rg.450 <- loadHDF5SummarizedExperiment(rg.450k.fname)
rg.epic <- loadHDF5SummarizedExperiment(rg.epic.fname)
# metadata
md.fname <- "mdf-epic-hm450k-merge_normal-blood_0-0-3.rda"
md <- get(load(md.fname))
md.cb <- md[md$blood.subgroup=="cord_blood",]

#gr.fname <- "gr-gseadj_h5se_hm450k-epic-merge_0-0-3"
#gr <- loadHDF5SummarizedExperiment(gr.fname)
#gr.cb <- gr[,gr$blood.subgroup=="cord_blood"]
#rm(gr)

#----------------
# est blood cells
#----------------
# get params
plat <- "IlluminaHumanMethylation450k"
bc.type <- "CordBlood"
rg.cb <- get(data(FlowSorted.CordBlood.450k))
bcv <- unique(rg.cb$CellType)

# est hm450k
rgf <- rg.450[,colnames(rg.450) %in% rownames(md.cb)]
rgf <- RGChannelSet(Green = as.matrix(getGreen(rgf)), 
                    Red = as.matrix(getRed(rgf)),
                    annotation = annotation(rgf))
est1 <- minfi::estimateCellCounts(rgf, compositeCellType = "CordBlood",
                           cellTypes = bcv, referencePlatform = plat)
save(est1, file = "cordblood-est_hm450k.rda")

# est epic
rgf <- rg.epic[,colnames(rg.epic) %in% rownames(md.cb)]
rgf <- RGChannelSet(Green = as.matrix(getGreen(rgf)), 
                    Red = as.matrix(getRed(rgf)),
                    annotation = annotation(rgf))
dim(rgf)
rgf <- convertArray(rgf, outType = plat)
dim(rgf)
est2 <- minfi::estimateCellCounts(rgf, compositeCellType = "CordBlood",
                                  cellTypes = bcv, referencePlatform = plat)
save(est2, file = "cordblood-est_epic.rda")

# save full estimates
cb.est <- as.data.frame(rbind(est1, est2))
cb.est.fname <- "cordblood-celltype-est_merged-hm450k-epic.rda"
save(cb.est, file = cb.est.fname)

#-----------------------------
# analyze estimated cell types
#-----------------------------
summary(rowSums(cb.est))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9755  1.0166  1.0252  1.0305  1.0393  1.1712

#-----------------------
# analyze nRBC estimates
#-----------------------
summary(cb.est$nRBC)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00000 0.02263 0.05319 0.07493 0.61884

# correlation with other cell types
ct <- cor(cb.est, method = "spearman")
