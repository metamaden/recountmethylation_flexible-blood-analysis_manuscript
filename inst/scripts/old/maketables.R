#!/usr/bin/env R

# Author: Sean Maden
# 
# Make tables for manuscript

#----------------
# cpg probe stats
#----------------
# load test data
library(minfiData); library(minfiDataEPIC)
data(RGsetEx); data(RGsetEPIC)

# get the probe annotations
anno.hm450k <- getAnnotation(RGsetEx); anno.epic <- getAnnotation(RGsetEPIC)
cnv <- c("chr", "pos", "strand", "Name", "AddressA", "AddressB", "Islands_Name", 
         "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")
annof.hm450k <- anno.hm450k[,cnv]; annof.epic <- anno.epic[,cnv]

# get the merged annotations
anno.out <- annof.hm450k[!annof.hm450k$Name %in% annof.epic$Name,]; dim(anno.out)
anno.merged <- rbind(annof.epic, anno.out); dim(anno.merged)

# detp summaries
dp.hm450k <- get(load("detp-sstats_hm450k-blood-groups.rda"))
dp.epic <- get(load("detp-sstats_epic-blood-groups.rda"))
colnames(dp.hm450k) <- paste0(colnames(dp.hm450k), ".hm450k")
colnames(dp.epic) <- paste0(colnames(dp.epic), ".epic")
rownames(dp.hm450k) <- dp.hm450k$cgid.all
dp.hm450k <- dp.hm450k[,!grepl("cgid", colnames(dp.hm450k))]
rownames(dp.epic) <- dp.epic$cgid.all
dp.epic <- dp.epic[,!grepl("cgid", colnames(dp.epic))]

# anovas
lff <- list.files(); lff <- lff[grepl("anova-df", lff)]
lav <- lapply(lff, function(x){as.data.frame(do.call(rbind, get(load(x))),
              stringsAsFactors = F)})
names(lav) <- paste0(gsub(".*blood-|.rda", "", lff), ".hm450k")
lav <- lapply(names(lav), function(x){avi <- lav[[x]]; 
  colnames(avi) <- c("cgid", paste0(x,".xssq"), paste0(x,".pval")); avi})
pav <- do.call(cbind, lav)
pav <- pav[!duplicated(pav[,1]),]; rownames(pav) <- pav[,1]
pav <- pav[,!grepl("cgid", colnames(pav))]

# anovas -- faster method
av <- get(load("anova-df_epic_blood-all.rda"))
vvar <- gsub(" ", "", rownames(av[[1]][[1]][[1]]))
mav <- do.call(rbind, lapply(seq(length(av)), function(ii){
  lavi <- av[[ii]]
  mavi <- do.call(rbind, lapply(seq(length(lavi)), function(x){
    cgid <- names(lavi)[x]; sstat <- lavi[[1]][[1]]
    xssqv <- paste0(vvar, ":", format(sstat[,3], digits = 3), collapse = ";")
    pvalv <- paste0(vvar, ":", format(sstat[,5], digits = 3), collapse = ";")
    return(matrix(c(cgid, xssqv, pvalv), nrow = 1))}))
  colnames(mavi) <- c("cgid", "xssq", "pval"); return(mavi)}))

