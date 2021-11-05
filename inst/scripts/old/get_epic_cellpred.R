#!/usr/bin/env R

# Author: Sean Maden

# Get blood cell type fraction predictions for EPIC arrays.

library(minfi)
rg.fname <- "remethdb_h5se_rg_epic-hm850k_merged_1621537799-1589820348_0-0-3"
rg <- HDF5Array::loadHDF5SummarizedExperiment(rg.fname)

save.dpath <- file.path("."); save.fname <- "celltypepred_epic.rda"
save.fpath <- file.path(save.dpath, save.fname)

num.samp <- 1; seqv <- seq(1, nrow(rg), num.samp); t1 <- Sys.time()

for(si in seqv){
  message("Starting at index ", si, "...")
  start.index <- si; end.index <- ifelse(si + num.samp - 1 > nrow(rg),
                                         nrow(rg), si + num.samp -1)
  rgf <- rg[, start.index:end.index] # sex and cell preds use raw red/grn signals
  rgf.matrix <- minfi::RGChannelSet(Green = as.matrix(minfi::getGreen(rgf)),
                                    Red = as.matrix(minfi::getRed(rgf)),
                                    annotation = BiocGenerics::annotation(rgf))
  cellpredi <- try(minfi::estimateCellCounts(rgf.matrix))
  if(class(cellpredi) == "try-error"){ # make new null matrix in case of error
    cellpredi <- as.data.frame(matrix(nrow = num.samp, ncol = 6))
    rownames(cellpredi) <- colnames(rgf.matrix)
    colnames(cellpredi) <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
    cellpredi$gsm <- colnames(rgf.matrix)}
  message("Finished cell type predictions.")
  cellpredi.df <- as.data.frame(cellpredi);cellpredi.df$gsm <- rownames(cellpredi.df)
  message("Saving cell type predictions.")
  appendvar <- ifelse(start.index == 1, FALSE, TRUE)
  data.table::fwrite(cellpredi.df, file = save.fpath, append = appendvar, 
                     sep = ",", row.names = F)
  message("Finished index ", si, ", time elapsed: ", Sys.time() - t1)
}

cellpred <- data.table::fread(file = save.fpath, sep = ",", data.table = F)


# error messages
# 1 samp -- Error in matrix(0, nSubj, nCol) : non-numeric matrix extent
# 2 samp -- Error in matrix(0, nSubj, nCol) : non-numeric matrix extent

# with hm450k anno


cellpred <- minfi::estimateCellCounts(rg)

# with epic anno
minfi::estimateCellCounts(rg,
                          referencePlatform="IlluminaHumanMethylationEPIC")


