#!/usr/bin/env R

# Author: Sean Maden
#
# Analyze predcell fraction similarity among blood subgroups.
#
#

library(ggplot)
library(ggdist)
library(gridExtra)
library(cowplot)

#----------
# load data
#----------
md <- get(load("si1_all-md-2platforms.rda"))
md2 <- get(load("si2_blood-md-2platforms.rda"))

#-----------------------
# get predcell summaries
#-----------------------
cnv <- colnames(md2)
cnv.predcell <- cnv[grepl("^predcell.*", cnv)]
catv <- gsub("predcell\\.", "", cnv.predcell)

groupv <- c("all", "cord_blood", "whole_blood",
            "peripheral_blood_mononuclear_cells")

# get all data
dfp <- do.call(rbind, lapply(catv, function(cati){
  do.call(rbind, lapply(groupv, function(groupi){
    if(groupi == "all"){which.samp <- seq(nrow(md2))} else{
      which.samp <- which(md2$blood_subgroup == groupi)}
    dati <- md2[which.samp, paste0("predcell.", cati)]
    matrix(c(rep(cati, length(dati)), 
             rep(groupi, length(dati)), 
             dati), ncol = 3)}))
})); colnames(dfp) <- c("predcell", "subgroup", "fraction")
dfp <- as.data.frame(dfp, stringsAsFactors = F)
dfp[,3] <- as.numeric(dfp[,3]);dfp$platform <- "2platforms"

# get data by platform
dfp.plat <- do.call(rbind, lapply(catv, function(cati){
  do.call(rbind, lapply(groupv, function(groupi){
    if(groupi == "all"){which.samp <- seq(nrow(md2))} else{
      which.samp <- which(md2$blood_subgroup == groupi)}
    do.call(rbind, lapply(c("hm450k", "epic"), function(platformi){
      dati <- md2[which.samp & md2$platform == platformi, 
                  paste0("predcell.", cati)]
      matrix(c(rep(cati, length(dati)), 
               rep(groupi, length(dati)), 
               dati, rep(platformi, length(dati))), 
             ncol = 4)}))}))}))
colnames(dfp.plat) <- c("predcell", "subgroup", "fraction", "platform")
dfp.plat <- as.data.frame(dfp.plat, stringsAsFactors = F)
dfp.plat <- rbind(dfp.plat, dfp)

# format vars
dfp.plat[,3] <- as.numeric(dfp.plat[,3])
dfp.plat[grepl("^peripheral.*", dfp.plat[,2]), 2] <- "PBMC"

# cell fraction plots

# define color blind palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

groupv <- c("all", "cord_blood", "whole_blood","PBMC")
lgg <- lapply(unique(dfp.plat$predcell), function(celli){
  dfpi <- dfp.plat[dfp.plat$predcell == celli,]
  # order dfpi on medians
  mediv <- sapply(groupv, function(x){
    which.dfpi <- dfpi[,2] == x & dfpi$platform == "2platforms"
    median(dfpi[which.dfpi,3])})
  dfpi[,2] <- factor(dfpi[,2], levels = names(mediv)[rev(order(mediv))])
  # get plot object
  ggplot(dfpi, aes(x = subgroup, y = fraction, fill = platform)) +
    geom_violin(draw_quantiles = 0.5) + coord_flip() + theme_bw() + 
    ggtitle(celli) + scale_fill_manual(values = cbp1) +
    xlab("") + ylab("") + theme(legend.position = "none")
}); names(lgg) <- unique(dfp.plat$predcell)

pleg <- ggplot(dfp.plat, aes(x = subgroup, y = fraction, fill = platform)) +
  geom_violin(draw_quantiles = 0.5) + scale_fill_manual(values = cbp1)
lgg[["legend"]] <- get_legend(pleg)

# make the composite plot
lm <- c(seq(1:6), rep(7,3))
plot.str <- paste0("lgg[[", seq(length(lgg)), "]]", collapse = ",")
plot.str <- paste0("grid.arrange(",plot.str,
                   ", layout_matrix = matrix(c(",
                   paste0(lm, collapse = ","),
                   "),ncol=3), bottom = 'Fraction', ",
                   "left = 'Blood group')")

pdf("violin-comp_predcell-bloodgroups_2platforms.pdf", 10, 10)
eval(parse(text = plot.str))
dev.off()

#-------------------
# search index corrs
#-------------------
# note: show blood cell fraction in relation to search index position
# do corr tests for each blood cell type fraction
# determine the cell types most predictive of search index position

cnv <- colnames(md)
cnv.predcell <- cnv[grepl("^predcell.*", cnv)]

si <- read.csv("sitest_results_blood-groups-2platforms.csv")
dim(si) # [1] 18671     8

#kv <- c(10,100,500,1000)
#ri <- 1; ki = 10
#dat <- si[ri,,drop = F]

si.cellcorr <- do.call(rbind, lapply(seq(nrow(si)), function(ri){
  message(ri); dat <- si[ri,,drop = F]
  do.call(rbind, lapply(kv, function(ki){
    dati <- as.character(dat[names(dat) == paste0("X",ki)])
    mdfii <- md[unlist(strsplit(dati, ";")),]
    samp.rank <- as.numeric(seq(ki))
    do.call(rbind, lapply(cnv.predcell, function(celli){
      samp.fracti <- as.numeric(mdfii[,celli])
      cti <- cor.test(samp.rank, samp.fracti, method = "spearman")
      matrix(c(ki, celli, cti$estimate, cti$p.value), nrow = 1)}))
  }))
}))
colnames(si.cellcorr) <- c("k", "predcell", "corr.rho", "corr.pnom")
si.cellcorr <- as.data.frame(si.cellcorr, stringsAsFactors = F)
for(c in c(1,3,4)){si.cellcorr[,c] <- as.numeric(si.cellcorr[,c])}
si.cellcorr$gsm.query <- rep(si[,2], each = length(kv)*6)
si.cellcorr$subgroup <- rep(si[,7], each = length(kv)*6)
si.cellcorr$platform <- rep(si[,8], each = length(kv)*6)

si.cellcorr.fname <- "si_cellpred_2platforms.rda"
save(si.cellcorr, file = si.cellcorr.fname)

# plot results
si.cellcorr <- get(load(si.cellcorr.fname))
si.cellcorr$abs.rho <- abs(si.cellcorr$corr.rho)

# plots: subgroups across k, for each cell type
# platformi <- "hm450k"

groupi <- "peripheral_blood_mononuclear_cells"
sif <- si.cellcorr[si.cellcorr$subgroup == groupi,]
sif[,1] <- factor(sif[,1])

lgg <- lapply(unique(sif$predcell), function(celli){
  sif.predcell <- sif[sif$predcell == celli,]
  ggplot(sif.predcell, aes(x = k, y = corr.rho, fill = platform)) +
    geom_violin(draw_quantiles = 0.5) + theme_bw() + ggtitle(celli)
}); names(lgg) <- unique(sif$predcell)

#-------------------------
# get bval data for new si
#-------------------------
library(minfi); library(HDF5Array); library(data.table)

# get paths
save.dpath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_results")
md.fpath <- file.path(save.dpath, "si2_blood-md-2platforms.rda")
gr.hm450k.fpath <- file.path("eternity", "recount-methylation", "recount-methylation-hm450k", "rmi",
                             "remethdb_hm450k_h5se_gr_merged_1619736651-1607018051_0-0-3")
gr.epic.fpath <- file.path("eternity", "recount-methylation", "recount-methylation-epic", "rmi",
                           "remethdb_h5se_gr_epic-hm850k_merged_1621537799-1607018051_0-0-3")
# define save path
save.fpath <- file.path(save.dpath, "bval-blood-merged-auto_row-gsm_2platforms.csv")

# make the bval matrix
md <- get(load(md.fpath))
gr.hm450k <- loadHDF5SummarizedExperiment(gr.hm450k.fpath)
gr.epic <- loadHDF5SummarizedExperiment(gr.epic.fpath)

# merge arrays
grf.hm450k <- gr.hm450k[,colnames(gr.hm450k) %in% rownames(md)]
grf.epic <- gr.epic[,colnames(gr.epic) %in% rownames(md)]
arr.type.str <- "IlluminaHumanMethylation450k"
gr.merge <- combineArrays(grf.hm450k, grf.epic, 
                          outType = arr.type.str)
dim(gr.merge) # [1] 453093  12858

# filter sex chr
which.cgsex <- getAnnotation(gr.merge)$chr %in% c("chrX", "chrY")
grf.merge <- gr.merge[!which.cgsex,]
dim(grf.merge) # [1] 442474  12858
table(getAnnotation(grf.merge)$chr)

# get bvals
bval <- t(getBeta(grf.merge)); dim(bval) # [1]  12858 442474
# write labels
samp.labels.fpath <- file.path(save.dpath, 
                               "bval-samplabels_blood-groups_2platforms.txt")
write(rownames(bval), sep = " ", file = samp.labels.fpath)
# write dnam bvals
num.samp <- 1000; indexv <- seq(1, nrow(bval), num.samp)
data.table::fwrite(as.data.frame(bval[1:num.samp,]), 
                   file = save.fpath, append = F,
                   col.names = F, row.names = F)
t1 <- Sys.time()
for(ri in indexv[2:length(indexv)]){
  start.index <- ri; end.index <- start.index + num.samp - 1
  end.index <- ifelse(end.index > nrow(bval), nrow(bval), end.index)
  # write.csv(dati, file = , append = T, row.names = F, col.names = F)
  data.table::fwrite(as.data.frame(bval[start.index:end.index,]), 
                     file = save.fpath, append = T, 
                     col.names = F, row.names = F)
  message("Finished writing row ", ri, ", time: ", Sys.time()-t1)
}
message("done")

# make si, query si, return results table
# see script: `si-blood-2platforms.py`

#------------------------------
# results, new si, with sex chr
#------------------------------
# si results for 1k nearest neighbors
md.fname <- "si2_blood-md-2platforms.rda"
md <- get(load(md.fname))
csv.fname <- "sitest_results_blood-groups-2platforms.csv"
tsi <- read.csv(csv.fname)
dim(tsi) # [1] 12857     4

# get stats for enrichment
sd.yr <- 5 # predage sd for interval check
# predcell colnames
predcell.cnv <- colnames(md)[grepl("predcell", colnames(md))]
#predcell stats labels
names.rvpc <- paste0(rep(predcell.cnv, each = 2), c("_rho", "_pnom"))
mr.knn <- do.call(rbind, lapply(seq(nrow(tsi)), function(ii){
  message(ii); gsmi <- tsi[ii,2]; which.gsm.query <- which(tsi[,2] == gsmi)
  knn.gsm.labelv <- unlist(strsplit(tsi[which.gsm.query,3], ";"))
  mdi <- md[knn.gsm.labelv,]
  which.gsm.md <- which(md$gsm == tsi[ii,2])
  md.gsmi <- md[which.gsm.md,]
  # study id
  gsei <- md.gsmi$gse
  gsei.ngsm.total <- nrow(md[md$gse == gsei,]) # get total study samples
  is.gei <- length(which(mdi$gse == gsei))
  fract.gsei <- is.gei/gsei.ngsm.total # get fraction knn in study
  # label
  gsm.group <- md.gsmi$blood_subgroup
  knn.group <- nrow(mdi[mdi$blood_subgroup == gsm.group,])
  # sex
  gsm.predsex <- md.gsmi$predsex
  knn.predsex <- nrow(mdi[mdi$predsex == gsm.predsex,])
  # age
  gsm.predage <- as.numeric(md.gsmi$predage)
  knn.predage <- as.numeric(mdi$predage)
  which.ageint <- knn.predage >= gsm.predage - sd.yr
  which.ageint <- which.ageint & knn.predage <= gsm.predage + sd.yr
  knn.ageint <- length(knn.predage[which.ageint])
  # corr tests with predcell
  rvpc <- unlist(lapply(predcell.cnv, function(cn){
    pc.gsm <- as.numeric(md.gsmi[,cn])
    pc.diff <- pc.gsm - as.numeric(mdi[,cn])
    cti <- suppressWarnings(cor.test(pc.diff, seq(1000), 
                                     method = "spearman"))
    return(c(cti$estimate, cti$p.value))}))
  names(rvpc) <- names.rvpc
  mdati <- c(gsmi, gsei.ngsm.total, fract.gsei, knn.group,
             knn.predsex, knn.ageint, rvpc)
  matrix(mdati, nrow = 1)
}))
colnames(mr.knn) <- c("gsm", "ngse_gsm", "fract_gse", "knn_group", "knn_predsex",
                      paste0("knn_ageint_", sd.yr), names.rvpc)
mr.knn <- as.data.frame(mr.knn, stringsAsFactors = F)
for(c in c(2:ncol(mr.knn))){mr.knn[,c] <- as.numeric(mr.knn[,c])}

mdf <- md[md$gsm %in% tsi[,2],]
mdf <- mdf[order(match(mdf$gsm, tsi[,2])),]
identical(mdf$gsm, tsi[,2])
mr.knn$blood_subgroup <- mdf$blood_subgroup
mr.knn$platform <- mdf$platform

library(ggplot2)

ggplot(mr.knn, aes(x = blood_subgroup, y = knn_predsex)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw()

ggplot(mr.knn, aes(x = blood_subgroup, y = predcell.CD4T_rho)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

ggplot(mr.knn, aes(x = blood_subgroup, y = predcell.CD8T_rho)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90))

md$predcell.CD8T <- as.numeric(md$predcell.CD8T)

#---------------------------
# results, new si, autosomal
#---------------------------
# si auto results for 1k nearest neighbors
csv.fname <- "sitest_results-auto_blood-groups-2platforms.csv"
tsi <- read.csv(csv.fname)
dim(tsi) # [1] 12857     4
md.fname <- "si2_blood-md-2platforms.rda"
md <- get(load(md.fname))

# get stats for enrichment
sd.yr <- 5 # predage sd for interval check
# predcell colnames
predcell.cnv <- colnames(md)[grepl("predcell", colnames(md))]
#predcell stats labels
names.rvpc <- paste0(rep(predcell.cnv, each = 2), c("_rho", "_pnom"))
mr.knn <- do.call(rbind, lapply(seq(nrow(tsi)), function(ii){
  message(ii); gsmi <- tsi[ii,2]; which.gsm.query <- which(tsi[,2] == gsmi)
  knn.gsm.labelv <- unlist(strsplit(tsi[which.gsm.query,3], ";"))
  mdi <- md[knn.gsm.labelv,]
  which.gsm.md <- which(md$gsm == tsi[ii,2])
  md.gsmi <- md[which.gsm.md,]
  # study id
  gsei <- md.gsmi$gse
  gsei.ngsm.total <- nrow(md[md$gse == gsei,]) # get total study samples
  is.gei <- length(which(mdi$gse == gsei))
  fract.gsei <- is.gei/gsei.ngsm.total # get fraction knn in study
  # label
  gsm.group <- md.gsmi$blood_subgroup
  knn.group <- nrow(mdi[mdi$blood_subgroup == gsm.group,])
  # sex
  gsm.predsex <- md.gsmi$predsex
  knn.predsex <- nrow(mdi[mdi$predsex == gsm.predsex,])
  # age
  gsm.predage <- as.numeric(md.gsmi$predage)
  knn.predage <- as.numeric(mdi$predage)
  which.ageint <- knn.predage >= gsm.predage - sd.yr
  which.ageint <- which.ageint & knn.predage <= gsm.predage + sd.yr
  knn.ageint <- length(knn.predage[which.ageint])
  # corr tests with predcell
  rvpc <- unlist(lapply(predcell.cnv, function(cn){
    pc.gsm <- as.numeric(md.gsmi[,cn])
    pc.diff <- pc.gsm - as.numeric(mdi[,cn])
    cti <- suppressWarnings(cor.test(pc.diff, seq(1000), 
                                     method = "spearman"))
    return(c(cti$estimate, cti$p.value))}))
  names(rvpc) <- names.rvpc
  mdati <- c(gsmi, gsei.ngsm.total, fract.gsei, knn.group,
             knn.predsex, knn.ageint, rvpc)
  matrix(mdati, nrow = 1)
}))
colnames(mr.knn) <- c("gsm", "ngse_gsm", "fract_gse", "knn_group", "knn_predsex",
                      paste0("knn_ageint_", sd.yr), names.rvpc)
mr.knn <- as.data.frame(mr.knn, stringsAsFactors = F)
for(c in c(2:ncol(mr.knn))){mr.knn[,c] <- as.numeric(mr.knn[,c])}

# append group and platform info
mdf <- md[md$gsm %in% tsi[,2],]
mdf <- mdf[order(match(mdf$gsm, tsi[,2])),]
identical(mdf$gsm, tsi[,2])
mr.knn$blood_subgroup <- mdf$blood_subgroup
mr.knn$platform <- mdf$platform

# analyze results -- predcell

get_dfp <- function(cnv, md){
  cnv <- c("platform", colnames(md)[grepl("predcell", colnames(md))])
  groupv <- c("whole_blood", "peripheral_blood_mononuclear_cells", 
              "cord_blood", "all")
  dfp <- do.call(rbind, lapply(groupv, function(groupi){
    if(groupi == "all"){
      dfp1 <- md[,cnv]
      dfp2 <- md[md$platform == "hm450k", cnv]
      dfp3 <- md[md$platform == "epic", cnv]
    } else{
      dfp1 <- md[md$blood_subgroup == groupi, cnv]
      dfp2 <- md[md$blood_subgroup == groupi & md$platform == "hm450k", cnv]
      dfp3 <- md[md$blood_subgroup == groupi & md$platform == "epic", cnv]
    }
    dfp2$platform <- "HM450k"; dfp3$platform <- "EPIC"; dfp1$platform <- "Combined"
    dfpi <- rbind(dfp1, rbind(dfp2, dfp3)); dfpi$group <- groupi
    return(dfpi)}))
}

make_lgg <- function(dfp, cnv, ylim.min = -1, ylim.max = 1){
  # convert pbmc labels
  dfp[grepl("peripheral", dfp$group),]$group <- "PBMC"
  # get predcell plots
  lgg <- lapply(cnv, function(cn){
    dfpi <- dfp;
    if(!cn %in% colnames(dfpi)){
      message("Couldn't find column ",cn, " in dfp...")
    } else{
      colnames(dfpi)[which(colnames(dfpi) == cn)] <- "predcell"
      dfpi[,"predcell"] <- as.numeric(dfpi[,"predcell"])
      lgg[[cn]] <- ggplot(dfpi, aes(x = group, y = predcell, color = platform)) + 
        geom_violin(draw_quantiles = 0.5) + theme_bw() +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90, 
                                         hjust = 0.95, vjust = 0.2)) + 
        ggtitle(cn) + ylim(ylim.min, ylim.max)
    }
  })
  # get the legend plot
  dfpi <- dfp;colnames(dfpi)[which(colnames(dfpi) == cnv[1])] <- "predcell"
  dfpi[,"predcell"] <- as.numeric(dfpi[,"predcell"])
  ggleg <- ggplot(dfpi, aes(x = group, y = predcell, color = platform)) + 
    geom_violin(draw_quantiles = 0.5) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(cn)
  lgg[["legend"]] <- cowplot::get_legend(ggleg)
  return(lgg)
}

plot_lgg <- function(pdf.fname, lgg, pdf.width = 10, pdf.height = 5,
                     xaxis_title = "Group", yaxis_title = "Predicted cell fraction"){
  # get plot str
  layout_str <- paste0("c(seq(6), 7, 7)")
  grid.str <- paste0("lgg[[", seq(length(lgg)),"]]", collapse = ",")
  grid.str <- paste0("grid.arrange(", grid.str, 
                     ", bottom = '",xaxis_title,
                     "', left = '",yaxis_title,"'",
                     ", layout_matrix = matrix(", layout_str, 
                     ", nrow = 2))")
  # save plot
  pdf(pdf.fname, 10, 5)
  eval(parse(text = grid.str)); dev.off()
  return(grid.str)
}

# plot rhos
md <- mr.knn; cnv <- colnames(md)[grepl("rho", colnames(md))]
dfp <- get_dfp(cnv, md)
lgg <- make_lgg(dfp, cnv, ylim.min = -0.75, ylim.max = 0.75)
pdf.fname <- "geom-violin-rho_col-plat_predcell-diff-knn_blood-groups-2plat.pdf"
plot_lgg(pdf.fname, lgg, yaxis_title = "Rho, rank vs. diff (query-k)")

# plot abs rhos
md <- mr.knn; cnv <- colnames(md)[grepl("rho", colnames(md))]
for(cn in cnv){
  md[,ncol(md) + 1] <- abs(md[,cn])
  colnames(md)[ncol(md)] <- paste0(cn, "_abs")}
cnv <- colnames(md)[grepl("rho_abs", colnames(md))]
dfp <- get_dfp(cnv, md)
lgg <- make_lgg(dfp, cnv, ylim.min = 0, ylim.max = 0.75)
pdf.fname <- "geom-violin-rho-abs_col-plat_predcell-diff-knn_blood-groups-2plat.pdf"
plot_lgg(pdf.fname, lgg, yaxis_title = "Absolute Rho, rank vs. diff (query-k)")


#---------------------------
# cell fractions by subgroup
#---------------------------
cnv <- c("platform", colnames(md)[grepl("predcell", colnames(md))])
groupv <- c("whole_blood", "peripheral_blood_mononuclear_cells", 
            "cord_blood", "all")

dfp <- do.call(rbind, lapply(groupv, function(groupi){
  if(groupi == "all"){
    dfp1 <- md[,cnv]
    dfp2 <- md[md$platform == "hm450k", cnv]
    dfp3 <- md[md$platform == "epic", cnv]
  } else{
    dfp1 <- md[md$blood_subgroup == groupi, cnv]
    dfp2 <- md[md$blood_subgroup == groupi & md$platform == "hm450k", cnv]
    dfp3 <- md[md$blood_subgroup == groupi & md$platform == "epic", cnv]
  }
  dfp2$platform <- "HM450k"; dfp3$platform <- "EPIC"; dfp1$platform <- "Combined"
  dfpi <- rbind(dfp1, rbind(dfp2, dfp3)); dfpi$group <- groupi
  return(dfpi)}))

dfp[grepl("peripheral", dfp$group),]$group <- "PBMC"

# get predcell plots
lgg <- lapply(colnames(dfp)[2:7], function(cn){
  dfpi <- dfp; which.cn <- which(colnames(dfpi) == cn)
  colnames(dfpi)[which.cn] <- "predcell"
  dfpi[,"predcell"] <- as.numeric(dfpi[,"predcell"])
  lgg[[cn]] <- ggplot(dfpi, aes(x = group, y = predcell, color = platform)) + 
    geom_violin(draw_quantiles = 0.5) + theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) + 
    ggtitle(cn) + ylim(0, 1)
})
# get the legend plot
ggleg <- ggplot(dfp, aes(x = group, y = predcell.CD8T, color = platform)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle(cn)
lgg[["legend"]] <- cowplot::get_legend(ggleg)
# get plot str
layout_str <- paste0("c(seq(6), 7, 7)")
grid.str <- paste0("lgg[[", seq(length(lgg)),"]]", collapse = ",")
grid.str <- paste0("grid.arrange(", grid.str, 
                   ", bottom = 'Group', left = 'Predicted cell fraction'",
                   ", layout_matrix = matrix(", layout_str, 
                   ", nrow = 2))")
# save plot
pdf("ggviolin-predcell_bloodgroups-2platforms.pdf", 10, 5)
eval(parse(text = grid.str)); dev.off()
