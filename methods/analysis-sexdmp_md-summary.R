#!/usr/bin/env R

# Author: Sean Maden
#
# Analyze sex DMP validation results from PBMC and whole blood. Get DMP annotation 
# data and summarize metadata terms for samples from validation of Inoshita et al 
# 2015 sex DMP study.
#
#

libv <- c("minfi", "sva", "gridExtra", "HDF5Array", "gaston", "ggrepel",
          "methyPre", "UpSetR", "ggdist", "ggpubr", "minfiData", "minfiDataEPIC")
sapply(libv, library, character.only = TRUE)

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

# load grant dmp sets
# load grant et al 2021 dmps
grant.csv.fname <- "grant-2021_wb-sex-dmp.csv"
grant.dmp <- read.csv(grant.csv.fname)
# inoshita et al 2015 dmps
dmp.study.fname <- "table1_sex-dmp_inoshita-2015.csv"
dmp.study <- read.csv(dmp.study.fname)

# get example rgsets
data("RGsetEx"); data("RGsetEPIC")

#-------------------
# supp table -- dmps
#-------------------
# get pbmc, wb dmps as top 1k most significant
tdf.pbmc <- tdf.pbmc[order(tdf.pbmc$ttest.pnom),]
tdf.wb <- tdf.wb[order(tdf.wb$ttest.pnom),]
# get cgid vectors of dmp sets
dmp.pbmc <- tdf.pbmc[seq(1000),1]
dmp.wb <- tdf.wb[seq(1000),1]
dmp.ino <- dmp.study$NAME
dmp.grant <- grant.dmp$X
# make new df
cgv <- unique(c(dmp.pbmc, dmp.wb, dmp.study$NAME))
cgdf <- data.frame(cgid = cgv, stringsAsFactors = F)
# get cpg categories
cgdf$is.whole.blood.dmp <- cgdf$is.pbmc.dmp <- cgdf$is.inoshita.2015.dmp <-
  cgdf$is.grant.2021.dmp <- FALSE
cgdf[cgdf[,1] %in% dmp.pbmc,]$is.pbmc.dmp <- TRUE
cgdf[cgdf[,1] %in% dmp.wb,]$is.whole.blood.dmp <- TRUE
cgdf[cgdf[,1] %in% dmp.ino,]$is.inoshita.2015.dmp <- TRUE
cgdf[cgdf[,1] %in% dmp.grant,]$is.grant.2021.dmp <- TRUE
# get the dnam directions (male - female), and ttest pvalues
cgdf$dnam.direction.whole.blood.dmp <- cgdf$dnam.direction.pbmc.dmp <-
  cgdf$dnam.direction.inoshita.2015.dmp <- cgdf$ttest.pvalue.whole.blood <-
  cgdf$ttest.pvalue.pbmc <- NA
for(ri in seq(nrow(cgdf))){
  cgid <- cgdf[ri,1]; message(cgid)
  if(cgid %in% dmp.wb){
    tdf.wb.filt <- tdf.wb[tdf.wb$cgid == cgid,,drop = FALSE]
    dnam.dir <- tdf.wb.filt$mean.male - tdf.wb.filt$mean.female
    cgdf$dnam.direction.whole.blood.dmp[ri] <- round(dnam.dir, 3)
    pvali <- as.numeric(format(tdf.wb.filt[2], scientific = T, digits = 3))
    cgdf$ttest.pvalue.whole.blood[ri] <- pvali}
  if(cgid %in% dmp.pbmc){
    tdf.pbmc.filt <- tdf.pbmc[tdf.pbmc$cgid == cgid,]
    dnam.dir <- tdf.pbmc.filt$mean.male - tdf.pbmc.filt$mean.female
    cgdf[ri,]$dnam.direction.pbmc.dmp <- round(dnam.dir, 3)
    pvali <- as.numeric(format(tdf.pbmc.filt[2], scientific = T, digits = 3))
    cgdf[ri,]$ttest.pvalue.pbmc <- pvali}
  if(cgid %in% dmp.ino){
    dmpi <- dmp.study[dmp.study$NAME == cgid,]
    dnam.dir <- as.numeric(dmpi[12])
    cgdf[ri,]$dnam.direction.inoshita.2015.dmp <- round(dnam.dir, 2)}
}

# append probe annotations
# bind anno for unique probes
anno1 <- getAnnotation(RGsetEx)
anno2 <- getAnnotation(RGsetEPIC)
anno1 <- anno1[,c(19, 24, 25, 26)]
anno2 <- anno2[,c(19, 22, 23, 24)]
anno2 <- anno2[!rownames(anno2) %in% rownames(anno1),]
anno <- rbind(anno1, anno2)
anno <- anno[rownames(anno) %in% cgdf$cgid,]
anno <- anno[order(match(rownames(anno), cgdf$cgid)),]
identical(rownames(anno), cgdf$cgid) # TRUE
cgdf$relation.to.island <- anno$Relation_to_Island
cgdf$ucsc.gene.names <- anno$UCSC_RefGene_Name
cgdf$ucsc.gene.groups <- anno$UCSC_RefGene_Group
cgdf$ucsc.gene.acc <- anno$UCSC_RefGene_Accession

# save
st.fname <- "st_sex-dmp_all-dmp-info"
write.csv(as.matrix(cgdf), file = paste0(st.fname, ".csv"))
write.table(as.matrix(cgdf), file = paste0(st.fname, ".tsv"))
save(cgdf, file = paste0(st.fname, ".rda"))

#-----------------------
# sample summary by sex
#-----------------------
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

#---------------
# summarize ages
#---------------
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
# age median: NA
# age sd: 21.2697547955506

# group: PBMC
# age median: 5.241
# age sd: 19.3089231855776