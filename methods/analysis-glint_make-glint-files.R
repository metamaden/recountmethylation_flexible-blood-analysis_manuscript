#!/usr/bin/env R

# Author: Sean Maden
#
# Prep bval and covariates matrix for GLINT, EPISTRUCTURE calculations. This 
# algorithm predicts population structure from HM450K DNAm (Rahmani et al 2017).
# This scripts writes DNAm/Beta-values for the the genotype-explanatory probes 
# provided in the manuscript supplement.
#
#

library(HDF5Array)
library(minfi)
library(data.table)
library(xlsx)

#----------
# load data
#----------
save.dpath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_final_results")
# get only the informative probes from Rahmani et al 2017 supp table
cgv.fname <- "13072_2016_108_MOESM3_ESM.xlsx"
cgv.fpath <- file.path(save.dpath, cgv.fname)
cgv <- read.xlsx(cgv.fpath, sheetIndex = 1)
cgv <- c(colnames(cgv), as.character(cgv[,1]))
# load grset
gr.fname <- "gr-gseadj_h5se_hm450k-epic-merge_0-0-3" # "gr-noob_h5se_hm450k-epic-merge_0-0-3"
gr.fpath <- file.path(save.dpath, gr.fname)
gr <- loadHDF5SummarizedExperiment(gr.fpath)
gr <- gr[rownames(gr) %in% cgv,]; dim(gr) # 4616 12242

#---------------------
# save files for glint
#---------------------
sepsym = "\t"; ext = "txt"

# write covariates matrix
mdf.fpath <- file.path(save.dpath, paste0("covariates-glint_4-blood-subgroups_2-platforms.",ext))
mdf <- as.data.frame(colData(gr))
# filter metadata
cnv <- c("platform", "predage", "predsex", "blood.subgroup", "predcell.CD8T", 
       "predcell.CD4T", "predcell.NK", "predcell.Bcell", "predcell.Mono", 
       "predcell.Gran")
mdf <- mdf[,cnv]
# format variables
mdf[,"predsex"] <- ifelse(mdf[,"predsex"] == "M", 1, 0)
mdf[,"platform"] <- ifelse(mdf[,"platform"] == "hm450k", 1, 0)
mdf[,"blood.subgroup"] <- ifelse(mdf[,"blood.subgroup"] == "PBMC", 1,
       ifelse(mdf[,"blood.subgroup"] == "whole_blood", 2,
              ifelse(mdf[,"blood.subgroup"] == "cord_blood", 3, 0)))
# write new covariates
mcnv <- matrix(nrow = 0, ncol = ncol(mdf)); colnames(mcnv) <- colnames(mdf)
fwrite(mcnv, file = mdf.fpath, sep = sepsym, row.names = F, col.names = T, append = F)
fwrite(mdf, file = mdf.fpath, sep = sepsym, row.names = T, col.names = F, append = T)

# write dnam matrix
bval.fname <- paste0("bval-gseadj-glint_4-blood-subgroups_2-platforms.", ext)
bval.fpath <- file.path(save.dpath, bval.fname)
# get filtered bval table; replace missing values with medians
mbval <- t(apply(as.matrix(getBeta(gr)),1,function(ri){
       ri[is.na(ri)] <- median(ri,na.rm=T);return(round(ri, 4))}))
rownames(mbval) <- rownames(gr)
# write table colnames, values
mcnv <- matrix(nrow = 0, ncol = ncol(gr)); colnames(mcnv) <- colnames(gr)
fwrite(mcnv, file = bval.fpath, sep = sepsym, row.names = F, col.names = T, append = F)
fwrite(mbval, file = bval.fpath, sep = sepsym, row.names = T, col.names = F, append = T)

#rl <- readLines(bval.fpath, n = 2)
#rl <- readLines(bval.fpath, n = 2)
