#!/usr/bin/env R

# Author: Sean Maden
#
# Get results for the Inoshita et al 2015 sex DMP validation in
# whole blood and PBMCs, compare, etc.

library(minfi); library(sva); library(gridExtra); 
library(HDF5Array); library(gaston); library(ggrepel)
library(methyPre); library(UpSetR); library(ggdist)
library(ggpubr)

#----------
# load data
#----------
# get all metadata for 2 platforms
md.fname <- "mdf-epic-hm450k-merge_normal-blood_0-0-3.rda"
md.fpath <- file.path(md.fname)
md <- get(load(file.path(md.fpath)))

# get study samples
study.acc <- "GSE67393"
mdf <- md[md$gse == study.acc,]

# ttest validation results
tdf.wb.fname <- "ttest-df-results-mvalfit_wholeblood-2platforms_inoshita-2015-validate.rda"
tdf.pbmc.fname <- "ttest-df-results-mvalfit_pbmc-2platforms_inoshita-2015-validate.rda"
tdf.wb <- get(load(tdf.wb.fname))
tdf.pbmc <- get(load(tdf.pbmc.fname))
# format tdf
for(c in c(2:4)){
  tdf.wb[,c] <- as.numeric(tdf.wb[,c])
  tdf.pbmc[,c] <- as.numeric(tdf.pbmc[,c])}

# load search index results
sit.fpath <- "sitest_results_blood-groups-2platforms.csv"
sit <- data.table::fread(sit.fpath, sep = ",", header = T, data.table = F)
sit <- sit[,c(2:ncol(sit))]
colnames(sit) <- c("gsm", "k=1000", "platform")

# get sample metadata
gseid <- "GSE67393"
md.filt <- md$gse == gseid | md$blood.subgroup == "whole_blood" | md$blood.subgroup == "PBMC"
mde <- md[md.filt,]
mde[mde$gse == gseid,]$blood.subgroup <- "Inoshita_et_al_2015"
lvlv <- c("Inoshita_et_al_2015", "whole_blood", "PBMC")
mde$blood.subgroup <- factor(mde$blood.subgroup, levels = lvlv)

# sample summary by sex
for(lvli in lvlv){
  mdf <- mde[mde$blood.subgroup == lvli,]
  message("group: ", lvli)
  message("females: ", length(which(mdf$predsex == "F")))
  message("males: ", length(which(mdf$predsex == "M")))
}

# group: Inoshita_et_al_2015
# females: 51
# males: 62

# group: whole_blood
# females: 3942
# males: 2924

# group: PBMC
# females: 397
# males: 230

# summarize ages
for(lvli in lvlv){
  mdf <- mde[mde$blood.subgroup == lvli,]
  message("group: ", lvli)
  message("age median: ", median(mdf$predage))
  message("age sd: ", sd(mdf$predage))
}
# group: Inoshita_et_al_2015
# age median: 47.402
# age sd: 12.2435232818108
# group: whole_blood
# age median: 35.709
# age sd: 21.2697547955506
# group: PBMC
# age median: 18.131
# age sd: 19.3089231855776