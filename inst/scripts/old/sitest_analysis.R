#!/usr/bin/env R

# Author: Sean Maden
#
# Analyze the results of search index tests on blood subgroups.

library(ggplot2)
library(ggdist)
library(tidyverse)
library(tidyquant)
library(gghalves)

#----------
# load data
#----------
md.all.fpath <- "si1_all-md-2platforms.rda"
md.all <- get(load(md.all.fpath))

# sitest csv
sit.fpath <- "sitest_results_blood-groups-2platforms.csv"
sit <- data.table::fread(sit.fpath, sep = ",", header = T, data.table = F)
sit <- sit[,c(2:ncol(sit))]
colnames(sit) <- c("gsm", paste0("k=", c(10, 100, 500, 1000)), 
                   "subgroup", "platform")

lgsmi.fpath <- "lgsm-search-index-results_2platform_blood-groups.rda"

# plot info
lcol <- list("hm450k" = "blue", "epic" = "darkgoldenrod2")
lshape <- list("hm450k" = 16, "epic" = 15)

#-----------------
# helper functions
#-----------------
# results by sample GSM id
gsmi_siresults <- function(sit, md.all, age.int = 5){
  # assess si results across k nearest neighbors
  # sit: Data frame containing search index query results 
  # (colnames: gsm, k = 10, k = 100, k = 500, k = 1000)
  # md.all: Metadata for all samples in search indices
  # age.int: Absolute age difference magnitude for grouping samples 
  # (e.g. samples +/-age.int matche queried sample's age)
  lall <- lapply(sit[,1], function(gsmid){
    mi.ref <- md.all[md.all$gsm == gsmid,]
    sii <- sit[sit[,1] == gsmid,]; kv <- sii[2:5]
    message("Processing sample ", gsmid, "...")
    lgsmi <- lapply(kv, function(x){
      ki <- unlist(strsplit(x, ";")); mi <- md.all[ki[!ki == rownames(mi.ref)],]
      # get query sample terms
      tv.ref<-unlist(strsplit(mi.ref$tissue,";"))
      dx.ref<-unlist(strsplit(mi.ref$disease,";"))
      predsex.ref <- mi.ref$predsex; predage.ref <- as.numeric(mi.ref$predage)
      predage.minrange <- ifelse(predage.ref > age.int, predage.ref - age.int, 0)
      predage.maxrange <- predage.ref + age.int
      stype.ref <- gsub("^msraptype:|;msrapconf.*", "", mi.ref$sampletype);kdenom <- nrow(mi)
      tv.fract <- sapply(tv.ref, function(x){
        length(which(grepl(paste0("(^|;)", x, "($|;)"), mi$tissue)))/kdenom})
      dx.fract <- sapply(dx.ref, function(x){
        length(which(grepl(paste0("(^|;)", x, "($|;)"), mi$disease)))/kdenom})
      predsex.fract <- length(which(mi$predsex == predsex.ref))/kdenom
      mi.predage.num <- as.numeric(mi$predage)
      predage.fract <- length(which(mi.predage.num >= predage.minrange & 
                                      mi.predage.num <= predage.maxrange))/kdenom
      mi.stype.form <- gsub("^msraptype:|;msrapconf.*", "", mi$sampletype)
      stype.fract <- length(which(mi.stype.form == stype.ref))/kdenom
      lgsmi.ki <- list(gsm = gsmid, k = kdenom+1, tv.fract = tv.fract, 
                       dx.fract = dx.fract, predsex.fract = predsex.fract, 
                       predage.fract = predage.fract, stype.fract = stype.fract)
      return(lgsmi.ki)}); message("Finished sample ", gsmid, ".")
    return(lgsmi)}); names(lall) <- sit[,1]
    message("Finished all samples in `sit`.")
    return(lall)
}

# get flat results matrix by term
ktermdat <- function(lgsmi, termlvl1 = "", termlvl2 = NA){
  # Get the k-wise identity frequencies for a term
  # lgsmi: List of gsm-wise term frequency results for k queries
  # termlvl1: Name of the first term level in a k query results list.
  # termlvl2: Optional name of the second term level in a k query results list,
  #           accessible as kresults[[termlvl1]][[termlvl2]].
  # Returns mv matrix of term frequencies by kvalue.
  if(is.na(termlvl2)){
    mv <- do.call(rbind, lapply(seq(length(lgsmi)), function(x){
      gsmi <- names(lgsmi)[x]; ldati <- lgsmi[[x]]
      message("Beginning gsm ", gsmi)
      kv <- lapply(ldati, function(y){
        return(y[[termlvl1]])})
      mv <- matrix(c(gsmi, kv), nrow = 1)
      message("Finished gsm ", gsmi)
      return(mv)}))} else{
        mv <- do.call(rbind, lapply(seq(length(lgsmi)), function(x){
          gsmi <- names(lgsmi)[x]; ldati <- lgsmi[[x]]
          message("Beginning gsm ", gsmi)
          kv <- lapply(ldati, function(y){
            kval <- "NA"
            if(termlvl1 %in% names(y)){
              datlvl1 <- y[[termlvl1]]
              if(termlvl2 %in% names(datlvl1)){
                kval <- datlvl1[[termlvl2]]}}
            return(kval)});mv <- matrix(c(gsmi, kv), nrow = 1)
          message("Finished gsm ", gsmi);return(mv)}))
      };colnames(mv) <- c("gsm", c(10, 100, 500, 1000)); return(mv)
}

#-----------------------
# subgroup label results
#-----------------------
lgsm.all <- gsmi_siresults(sit, md.all, 5)
lgsmi.fpath <- "lgsm-search-index-results_2platform_blood-groups.rda"
save(lgsm.all, file = lgsmi.fpath)

#---------------------------
# results matrix by subgroup
#---------------------------
lgsm.all <- get(load(lgsmi.fpath))
groupv <- unique(sit[[6]])
platformv <- unique(sit$platform)
lmv <- list()
for(groupi in groupv){
  message("Beginning group ", groupi)
  for(platformi in platformv){
    message("Beginning platform ", platformi)
    gsmv <- sit[[1]][which(sit[[6]] == groupi & sit[[7]] == platformi)]
    cond <- names(lgsm.all) %in% gsmv & !duplicated(names(lgsm.all))
    lgsmi <- lgsm.all[cond];kv <- c("10", "100", "500", "1000")
    label <- ifelse(groupi == "all", "blood", 
                    ifelse(groupi == "cord_blood", 
                           "umbilical_cord", groupi))
    mv <- do.call(rbind, lapply(seq(length(lgsmi)), function(x){
      gsmi <- names(lgsmi)[x]; xv <- lgsmi[[gsmi]]
      resv <- lapply(xv, function(y){
        tvi <- y[["tv.fract"]]; tvii <- tvi[names(tvi) == label]
        return(tvii)});mv <- matrix(c(gsmi, resv), nrow = 1)
      return(mv)}));colnames(mv) <- c("gsm", kv)
    lmv[[paste0(label, "_", platformi)]] <- mv
    message("Finished with group ", groupi)}}

#-------------
# plot results
#-------------
lgsmi <- get(load(lgsmi.fpath))
ldfp <- list() # list the dfp objects for plots
for(groupi in c("umbilical_cord", "whole_blood", "blood",
                "peripheral_blood_mononuclear_cells")){
  ldfp[[groupi]] <- do.call(rbind, lapply(c("hm450k", "epic"), function(platformi){
    dfi <- lmv[[paste0(groupi,"_",platformi)]]
    dfpi <- do.call(rbind, lapply(c(10, 100, 500, 1000), function(x){
      dati <- dfi[,colnames(dfi) == as.character(x), drop = F]
      mvi.dat <- c(dati[,1], rep(x, nrow(dati)), rep(platformi, nrow(dati)))
      mvi <- as.data.frame(matrix(mvi.dat, ncol = 3));return(mvi)}))
    colnames(dfpi) <- c("id.freq", "kval", "platform")
    dfpi[,1] <- as.numeric(dfpi[,1])
    dfpi[,2] <- as.factor(as.numeric(dfpi[,2]))
    dfpi[,3] <- as.factor(as.character(dfpi[,3]))
    return(dfpi)}))}

# make the individual plots
lplot <- list()
for(ii in seq(length(ldfp))){ # ii = 3
  dfp <- ldfp[[ii]]; term <- names(ldfp)[ii]
  term <- ifelse(term == "peripheral_blood_mononuclear_cells", "PBMC", term)
  for(platform in c("hm450k", "epic")){
    pdf.fpath <- paste0("ggstorm_si-test_platform-",
                        platform,"-term-",term,".pdf")
    color <- ifelse(platform == "hm450k", "blue", "darkgoldenrod2")
    ptshape <- ifelse(platform == "hm450k", 16, 15)
    plot.title <- paste0(term, ", ", platform)
    lplot[[plot.title]] <- ggplot(dfp[dfp$platform == platform,], aes(x = kval, y = id.freq)) + 
      ggdist::stat_halfeye(adjust = .3, width = 0.5, .width = 0, 
                           color = color, fill = color,
                           justification = -.3, point_colour = color) + 
      geom_boxplot(width = .18, outlier.shape = NA, color = color) +
      gghalves::geom_half_point(side = "l", range_scale = .21, 
                                alpha = .2, color = color, shape = ptshape) +
      theme_bw() + xlab("K nearest neighbors") + ylab("Fraction with term") +
      ggtitle(plot.title)
    pdf(pdf.fpath, 3.6, 2); print(lplot[[plot.title]]); dev.off()}}

# make a grid of plots
lplot <- list()
for(ii in seq(length(ldfp))){ # ii = 3
  dfp <- ldfp[[ii]]; term <- names(ldfp)[ii]
  term <- ifelse(term == "peripheral_blood_mononuclear_cells", "PBMC", term)
  for(platform in c("hm450k", "epic")){
    ylab.txt <- ifelse(ii == 1, "HM450k", ifelse(ii == 2, "EPIC", ""))
    pdf.fpath <- paste0("ggstorm_si-test_platform-",
                        platform,"-term-",term,".pdf")
    color <- ifelse(platform == "hm450k", "blue", "darkgoldenrod2")
    ptshape <- ifelse(platform == "hm450k", 16, 15)
    plot.title <- paste0(term, ", ", platform)
    lplot[[plot.title]] <- ggplot(dfp[dfp$platform == platform,], aes(x = kval, y = id.freq)) + 
      ggdist::stat_halfeye(adjust = .3, width = 0.5, .width = 0, 
                           color = color, fill = color,
                           justification = -.3, point_colour = color) + 
      geom_boxplot(width = .18, outlier.shape = NA, color = color) +
      gghalves::geom_half_point(side = "l", range_scale = .21, 
                                alpha = .2, color = color, shape = ptshape) +
      theme_bw() + xlab("K nearest neighbors") + ylab(ylab.txt) +
      ggtitle(plot.title)
  }}

library(gridExtra)

plot.txt <- paste0("lplot[[", seq(length(lplot)), "]]", collapse = ",")
plot.txt <- paste0("grid.arrange(", plot.txt,
                   ", layout_matrix = matrix(nrow = 2, ncol = 4),",
                   "left = 'Fraction with term', bottom = 'K')")
eval(parse(text = paste0()))

grid.arrange()
pdf(pdf.fpath, 3.6, 2); print(lplot[[plot.title]]); dev.off()

#-------------------------
# results for common terms
#-------------------------
var.name <- "predsex.fract"
lgsmi <- lgsm.all[!duplicated(names(lgsm.all))]
mv <- do.call(rbind, lapply(seq(length(lgsmi)), function(x){
  gsmi <- names(lgsmi)[x]; ldati <- lgsmi[[x]]
  kv <- lapply(ldati, function(y){
    return(y[[var.name]])})
  mv <- matrix(c(gsmi, kv), nrow = 1)
  return(mv)}))
colnames(mv) <- c("gsm", c(10, 100, 500, 1000))

mv <- ktermdat(lgsmi[!duplicated(names(lgsmi))], termlvl1 = "predsex.fract")
mv <- ktermdat(lgsmi[!duplicated(names(lgsmi))], termlvl1 = "predage.fract")
mv <- ktermdat(lgsmi[!duplicated(names(lgsmi))], termlvl1 = "dx.fract", termlvl2 = "normal")
mv <- ktermdat(lgsmi[!duplicated(names(lgsmi))], termlvl1 = "dx.fract", termlvl2 = "control")
