#!/usr/bin/env R

# Author: Sean Maden
#
# Simulate the result of GSE bias corrections on variances, inc.
# a general GSE bias correction across available studies (adj1) and
# a precise GSE bias correction on a subset of studies evaluated (adj2).
# For comparison, the variances from unadjusted DNAm is also calculated.
# 
# Main script steps: 1. get gr.all; 2. filter gr.all on blood samples, 
# autosomal probes, non-crossreactive probes; 3. get the test group sample
# subset; 4. filter probes on percentage detp > 0.01 == 0 for the group;
# 5. run the simulations and evaluations across iterations of probes and studies.
# The same probe set is tested in randomized study groups of varying sizes.

library(HDF5Array); library(minfi)
library(methyPre); library(limma); library(sva)

#----------
# load data
#----------
# get dir paths
save.dpath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_results")
load.hm450k.dpath <- file.path("eternity", "recount-methylation", 
                               "recount-methylation-hm450k", "rmi")
load.epic.dpath <- file.path("eternity", "recount-methylation", 
                             "recount-methylation-epic", "rmi")

# load metadata
md.fname <- "si2_blood-md-2platforms.rda"
md <- get(load(file.path(save.dpath, md.fname)))

# load se objects
gr.hm450k.fname <- "remethdb_hm450k_h5se_gr_merged_1619736651-1607018051_0-0-3"
gr.epic.fname <- "remethdb_h5se_gr_epic-hm850k_merged_1621537799-1607018051_0-0-3"
gr.hm450k <- loadHDF5SummarizedExperiment(file.path(load.hm450k.dpath, gr.hm450k.fname))
gr.epic <- loadHDF5SummarizedExperiment(file.path(load.epic.dpath, gr.epic.fname))

# load the detp tables
ptable1.path <- file.path(save.dpath, "detp-sstats_hm450k-blood-groups.rda")
ptable2.path <- file.path(save.dpath, "detp-sstats_epic-blood-groups.rda")
ptable1 <- get(load(ptable1.path))
ptable2 <- get(load(ptable2.path))

#----------------------------------
# helper functions, parallel method
#----------------------------------
# get anova data from model string
get_aovdat <- function(aov.str, labstr = "", pheno_subset, 
                       orderv = c("gse", "predage", "predcell.CD8T", "predcell.CD4T", 
                                  "predcell.NK", "predcell.Bcell", "predcell.Mono", 
                                  "predcell.Gran", "Residuals", "platform", "predsex")){
  xaov <- eval(parse(text = aov.str))
  xdat <- xaov[[1]]
  namev <- gsub(" ", "", rownames(xdat[1]))
  xdat.iter <- xdat[,c(3:5)]
  typev <- c("mssq", "fval", "pval")
  vdat <- unlist(lapply(1:3, function(ii){
    datii <- xdat.iter[,ii]; names(datii) <- namev
    for(vname in c("predsex", "platform")){
      if(!vname %in% namev){
        datii <- c(datii, "NA"); names(datii)[length(datii)] <- vname}}
    datii <- datii[order(match(names(datii), orderv))]
    names(datii) <- paste0(names(datii), "_", typev[ii])
    return(datii)}))
  names(vdat) <- paste0(names(vdat), "_", labstr)
  return(vdat)}

# parallel function for gse reps
par_gserep <- function(gse_rep, num.studies.subset, gsev, pheno, cgidv){
  # get study subset
  message("Working on rep ", gse_rep, " with ", num.studies.subset," studies...")
  gsev_filt <- sample(gsev, num.studies.subset)
  message("Using studies ", paste0(gsev_filt, collapse = ", "), "...")
  pheno_subset <- pheno[pheno$gse %in% gsev_filt,]
  pheno_subset$gse <- droplevels(pheno_subset$gse)
  # get the adj2 data
  which.bunadj <- which(grepl("_unadj", colnames(pheno_subset)))
  bunadj_subset <- t(as.matrix(pheno_subset[,which.bunadj]))
  rownames(bunadj_subset) <- paste0(gsub("_unadj", "", 
                                         rownames(bunadj_subset)), "_adj2")
  bval_adj2 <- removeBatchEffect(bunadj_subset, batch = pheno_subset$gse)
  pheno_subset <- cbind(pheno_subset, t(bval_adj2))
  # get the vector of model strings
  cnv_numeric <- colnames(pheno_subset)[grepl("predcell", colnames(pheno_subset))]
  cnv_numeric <- c("predage", cnv_numeric)
  lm.str <- paste0("~ ", paste0(cnv_numeric, collapse = " + "))
  num.lvl.predsex <- length(unique(pheno_subset$predsex))
  num.lvl.platform <- length(unique(pheno_subset$platform))
  if(num.lvl.predsex > 1){lm.str <- paste0(lm.str, " + predsex")}
  if(num.lvl.platform > 1){lm.str <- paste0(lm.str, " + platform")}
  # get results for 3 adj levels as a matrix
  cnv_all <- colnames(pheno_subset)
  mrep <- do.call(cbind, lapply(c("unadj", "adj", "adj2"), function(labstri){
    message("Working on models of type '", labstri, "'...")
    cg_cnv <- cnv_all[grepl(paste0(".*_", labstri, "$"), cnv_all)]
    lm.str.vect <- paste0(cg_cnv, " ", lm.str, ", data = pheno_subset")
    aov.str.vect <- paste0("summary(aov(", lm.str.vect, "))")
    mrep <- do.call(rbind, lapply(aov.str.vect, get_aovdat, 
                                  pheno_subset = pheno_subset, 
                                  labstr = labstri))
    mrep <- as.data.frame(mrep, stringsAsFactors = F)
    colnames(mrep) <- paste0(colnames(mrep), "_", labstri)
    return(mrep)}))
  mrep <- as.data.frame(mrep, stringsAsFactors = F)
  mrep$gse_rep <- gse_rep
  mrep$ngse <- num.studies.subset
  mrep$gsev <- paste0(unique(pheno_subset$gse), collapse = ";")
  mrep$cgid <- cgidv
  return(mrep)
}

# get all results for a series of study reps
get_mgserep <- function(ngse_rep, num.studies.subset, cgidv, pheno){
  gsev <- unique(pheno$gse) # get studies
  # format pheno colnames
  cnv_numeric <- c("predage", colnames(pheno)[grepl("predcell", colnames(pheno))])
  cnv_factor <- c("predsex", "platform", "gse")
  for(c in cnv_numeric){pheno[,c] <- as.numeric(pheno[,c])}
  for(c in cnv_factor){pheno[,c] <- as.factor(pheno[,c])}
  # get results matrix
  # process reps in parallel
  message("Processing studies using ", ngse_rep, " cores...")
  mgserep <- do.call(rbind, mclapply(seq(ngse_rep), par_gserep, 
                                     num.studies.subset = num.studies.subset,
                                     gsev = gsev, pheno = pheno, cgidv = cgidv, 
                                     mc.cores = ngse_rep))
  return(mgserep)
}

# get the results of all reps of random study selections
get_mgse_all <- function(grf, ngsev, num.gse = 5, ngse_rep = 10, num.probes = 500){
  message("Getting ",num.probes," random probes and ", 
          num.gse, " random studies...")
  cgidv <- sample(rownames(grf), num.probes)
  gsev <- sample(unique(grf$gse), num.gse)
  grff <- grf[cgidv, colnames(grf[,grf$gse %in% gsev])]
  bval_unadj <- getBeta(grff)
  # bval_unadj <- na.omit(bval_unadj) # handle NAs
  mval_unadj <- logit2(bval_unadj)
  message("Getting adjustment 1 data...")
  mval_adj <- removeBatchEffect(mval_unadj, 
                                batch = grff$gse)
  bval_adj <- ilogit2(mval_adj)
  message("Appending DNAm to pheno...")
  pheno <- colData(grff)
  rownames(bval_unadj) <- paste0(rownames(bval_unadj), "_unadj")
  rownames(bval_adj) <- paste0(rownames(bval_adj), "_adj")
  pheno <- cbind(pheno, 
                 cbind(t(as.matrix(bval_unadj)), 
                       t(as.matrix(bval_adj))))
  message("Beginning iterations of random study selection and ANVOAs...")
  mgse.all <- do.call(rbind, lapply(ngsev, function(num.studies.subset){
    get_mgserep(ngse_rep = ngse_rep, num.studies.subset = num.studies.subset, 
                cgidv = cgidv, pheno = pheno)}))
  message("Finished all ANOVA iterations; returning...")
  return(mgse.all)
}

# parallel wrapper for `get_mgse_all()`
par_mgse_all <- function(cgrep, grf, ngsev, num.gse = 5, 
                         ngse_rep = 10, num.probes = 500){
  t1 <- Sys.time(); message("Beginning cgrep ", cgrep, "...")
  mgse_all <- get_mgse_all(grf = grf, ngsev = ngsev, num.gse = num.gse, 
                           ngse_rep = ngse_rep, num.probes = num.probes)
  mgse_all <- as.data.frame(mgse_all, stringsAsFactors = F)
  message("Finished cgrep ", cgrep, ", time: ", Sys.time() - t1)
  mgse_all$cgrep <- cgrep; return(mgse_all)
}

#--------------------------------
# filter gr.all, samples & probes
#--------------------------------
# combine arrays
gr.all <- combineArrays(gr.epic, gr.hm450k, 
                        outType = "IlluminaHumanMethylation450k")
dim(gr.all) # [1] 453093  65260

# filter samples
gr.all <- gr.all[,colnames(gr.all) %in% rownames(md)]
dim(gr.all) # [1] 453093  12858
# match md
md <- md[order(match(rownames(md), colnames(gr.all))),]
cond <- identical(rownames(md), colnames(gr.all))
cond # TRUE
colData(gr.all) <- DataFrame(md)

# filter probes
anno <- getAnnotation(gr.all)
cg.keep <- rownames(anno[!anno$chr %in% c("chrY", "chrX"),])
data(chen_crxcg); cg.keep <- cg.keep[!cg.keep %in% chen.crxcg]
gr.all <- gr.all[rownames(gr.all) %in% cg.keep,]
dim(gr.all) # [1] 416045  12858

#----------------
# PBMC bias tests  
#----------------
# pbmc
grf <- gr.all

# get file info
groupi <- "peripheral_blood_mononuclear_cells"; groupi.label <- "pbmc"
save.fname <- paste0("mgse-adj-test_blood-group-",
                     groupi.label, "_2platforms.rda")
save.fpath <- file.path(save.dpath, save.fname)

# get the groupi sample ids
gsm.labelv <- rownames(md[md$blood_subgroup == groupi,])
grf <- grf[,colnames(grf) %in% gsm.labelv]
dim(grf) # [1] 399968    642
table(grf$platform)
# epic hm450k
# 256    386

# get detp filtered cgs
cname <- "perc_above_01.peripheral_blood_mononuclear_cells"
cg.keep <- unique(ptable1[ptable1[,cname] == 0,1],
                  ptable2[ptable2[,cname] == 0,1])
grf <- grf[rownames(grf) %in% cg.keep,]
dim(grf) # [1] 399968  642

# count studies
length(unique(grf$gse)) # 9
ngsev <- c(2, 3, 4, 5) # gse quantities vector for sim's

# run successive cg reps in batches
num_cgrep <- 5
ngsev <- c(2, 3, 4)
num.gse <- 5
ngse_rep <- 2
num.probes <- 100
num.batch <- 30
t1 <- Sys.time(); set.seed(0)
mcg <- do.call(rbind, lapply(seq(num.batch), function(batch){
  message("Beginning batch ",batch, ", time: ", Sys.time() - t1)
  lmcg_batch <- mclapply(seq(num_cgrep), par_mgse_all, 
                         grf = grf, ngsev = ngsev, num.gse = num.gse, 
                         ngse_rep = ngse_rep, num.probes = num.probes,
                         mc.cores = num_cgrep)
  mcg_batch <- matrix(ncol = 95, nrow = 0)
  for(ii in seq(length(lmcg_batch))){
    mcgii <- lmcg_batch[[ii]]
    if(is(mcgii, "data.frame")){
      if(ncol(mcgii) == 95){
        mcg_batch <- rbind(mcg_batch, mcgii)}}}
  mcg_batch <- as.data.frame(mcg_batch, stringsAsFactors = F)
  message("Finished batch ", batch, ", time: ", Sys.time() - t1)
  return(mcg_batch)}))

# save results matrix
save(mcg, file = save.fpath)

#-----------------------
# whole blood bias tests
#-----------------------
# whole blood

# get file info
groupi <- "whole_blood"; groupi.label <- "whole-blood"
save.fname <- paste0("mgse-adj-test_blood-group-",
                     groupi.label, "_2platforms.rda")
save.fpath <- file.path(save.dpath, save.fname)

# get detp filtered cgs
cg.keep <- unique(ptable1[ptable1$perc_above_01.whole_blood == 0,1],
                  ptable2[ptable2$perc_above_01.whole_blood == 0,1])
grf <- gr.all[rownames(gr.all) %in% cg.keep,]
dim(grf) # [1] 416045  12858

# filter na's
bval <- as.numeric(getBeta(grf[,1]))
cgidv.na <- rownames(grf)[which(is.na(bval))]
length(cgidv.na) # [1] 545
grf <- grf[!rownames(grf) %in% cgidv.na,]
dim(grf) # 415500   5757

# get the groupi sample ids
gsm.labelv <- rownames(md[md$blood_subgroup == groupi,])
grf <- grf[,colnames(grf) %in% gsm.labelv]
dim(grf) # [1] 416045   5980
length(unique(grf$gse)) # 25

# run successive cg reps in batches
num_cgrep <- 3
ngsev <- c(2, 3, 4)
num.gse <- 5
ngse_rep <- 2
num.probes <- 100
num.batch <- 30
t1 <- Sys.time(); set.seed(0)
mcg <- do.call(rbind, lapply(seq(num.batch), function(batch){
  message("Beginning batch ",batch, ", time: ", Sys.time() - t1)
  mcg_batch <- do.call(rbind, mclapply(seq(num_cgrep), par_mgse_all, 
                                       grf = grf, ngsev = ngsev, num.gse = num.gse, 
                                       ngse_rep = ngse_rep, num.probes = num.probes,
                                       mc.cores = num_cgrep))
  message("Finished batch ",batch, ", time: ", Sys.time() - t1)
  return(mcg_batch)}))

# save results matrix
save(mcg, file = save.fpath)

#----------------------
# cord blood bias tests
#----------------------
# cord blood
# get file info
groupi <- "cord_blood"; groupi.label <- "cord-blood"
save.fname <- paste0("mgse-adj-test_blood-group-",
                     groupi.label, "_2platforms.rda")
save.fpath <- file.path(save.dpath, save.fname)

# get detp filtered cgs
cg.keep <- unique(ptable1[ptable1$perc_above_01.cord_blood == 0,1],
                  ptable2[ptable2$perc_above_01.cord_blood == 0,1])
grf <- gr.all[rownames(gr.all) %in% cg.keep,]

# get the groupi sample ids
gsm.labelv <- rownames(md[md$blood_subgroup == groupi,])
grf <- grf[,colnames(grf) %in% gsm.labelv]

# run successive cg reps in batches
num_cgrep <- 3; ngsev <- c(2, 3, 4); num.gse <- 5; ngse_rep <- 2
num.probes <- 100; num.batch <- 30; t1 <- Sys.time(); set.seed(0)
mcg <- do.call(rbind, lapply(seq(num.batch), function(batch){
  message("Beginning batch ",batch, ", time: ", Sys.time() - t1)
  mcg_batch <- do.call(rbind, mclapply(seq(num_cgrep), par_mgse_all, grf = grf, ngsev = ngsev, 
                                       num.gse = num.gse, ngse_rep = ngse_rep, num.probes = num.probes,
                                       mc.cores = num_cgrep))
  message("Finished batch ",batch, ", time: ", Sys.time() - t1); return(mcg_batch)}))
# save results matrix
save(mcg, file = save.fpath)

#---------------------
# all blood bias tests
#---------------------
# all blood
# get file info
groupi <- groupi.label <- "all"
save.fname <- paste0("mgse-adj-test_blood-group-",
                     groupi.label, "_2platforms.rda")
save.fpath <- file.path(save.dpath, save.fname)

# filter detp 
cg.keep <- unique(ptable1[ptable1$perc_above_01.all == 0,1],
  ptable2[ptable2$perc_above_01.all == 0,1])
grf <- gr.all[rownames(gr.all) %in% cg.keep,]
dim(grf) # [1] 286798  12858

# run successive cg reps in batches
num_cgrep <- 3; ngsev <- c(2, 3, 4); num.gse <- 5; ngse_rep <- 2
num.probes <- 100; num.batch <- 30; t1 <- Sys.time(); set.seed(0)
mcg <- do.call(rbind, lapply(seq(num.batch), function(batch){
  message("Beginning batch ",batch, ", time: ", Sys.time() - t1)
  mcg_batch <- do.call(rbind, mclapply(seq(num_cgrep), par_mgse_all, grf = grf, ngsev = ngsev, 
                                       num.gse = num.gse, ngse_rep = ngse_rep, 
                                       num.probes = num.probes, mc.cores = num_cgrep))
  message("Finished batch ",batch, ", time: ", Sys.time() - t1); return(mcg_batch)}))
save(mcg, file = save.fpath)

#---------------------
# get results matrices
#---------------------
fnv <- c("mgse-adj-test_blood-group-cord-blood_2platforms.rda",
         "mgse-adj-test_blood-group-whole-blood_2platforms.rda",
         "mgse-adj-test_blood-group-pbmc_2platforms.rda",
         "mgse-adj-test_blood-group-all_2platforms.rda")
lmcg <- lapply(fnv, function(fn){get(load(fn))})
names(lmcg) <- gsub(".*group-|_2platforms.*", "", fnv)

# variables and model types
varv <- unique(gsub("_.*", "", colnames(lmcg[[1]])))[1:10]
typev <- c("unadj", "adj1", "adj2")

# fract explained variances (FEV) by probe, group
cnames.mfsq <- paste0(varv, "_", rep(typev, each = length(varv)))
lindex <- list("unadj" = 1:10, "adj" = 31:40, "adj2" = 61:70) # model type var indices
# get the list of fract mssq by group
lfmsq <- lapply(lmcg, function(lmcgi){
  mfsq <- t(apply(lmcgi, 1, function(dati){
    fract.rsq <- do.call(cbind, lapply(lindex, function(indexv){
      datii <- dati[indexv]; tot.rsq <- sum(as.numeric(datii), na.rm = t)
      fract.rsq <- as.numeric(datii)/tot.rsq; names(fract.rsq) <- names(datii)
      return(datii)})); return(fract.rsq)})); colnames(mfsq) <- cnames.mfsq
  mfsq <- as.data.frame(mfsq, stringsAsFactors = F)
  mfsq$ngse <- lmcgi$ngse; return(mfsq)}); names(lfmsq) <- names(lmcg)

# get 3x FEV differences
cname.diff <- paste0(rep(varv, each = 3), "_", 
                     rep(c("diff1", "diff2", "diff3"), times = length(varv)))
ldif <- lapply(lfmsq, function(mfsq){
  mdiff.all <- do.call(cbind, lapply(varv, function(vari){
    fii <- mfsq[,grepl(vari, colnames(mfsq))];namev <- names(fii)
    diff1 <- as.numeric(fii[,1]) - as.numeric(fii[,2])
    diff2 <- as.numeric(fii[,1]) - as.numeric(fii[,3])
    diff3 <- as.numeric(fii[,2]) - as.numeric(fii[,3])
    mdiff <- cbind(diff1, cbind(diff2, diff3))
    return(mdiff)}));mdiff.all <- as.data.frame(mdiff.all, stringsAsFactors = F)
    colnames(mdiff.all) <- cname.diff; mdiff.all$ngse <- mfsq$ngse
    return(mdiff.all)}); names(ldif) <- names(lfmsq)

save(ldif, file = "list_gse-bias-diffs_all-tx.rda")

#-------------------------------
# heatmap and supp table -- FEVs
#-------------------------------
# bind data
mfev <- do.call(rbind, lapply(names(lfmsq), function(txi){
  li <- lfmsq[[txi]]; li$tx <- txi; return(li)}))
# save
save(mfev, file = "matrix-fev_gse-bias_all-tx.rda")

# fev medians, all combined
mfev <- get(load("matrix-fev_gse-bias_all-tx.rda"))
med.fev <- apply(mfev[,c(1:30)],2,function(ci){
  median(as.numeric(ci), na.rm = T)})
dfp <- data.frame(label = names(med.fev), 
                  median.fev = as.numeric(med.fev),
                  stringsAsFactors = F)
dfp$fev.type <- gsub(".*_", "", dfp$label)
dfp$variable <- gsub("_.*", "", dfp$label)

# heatmap, combined fev medians
# format vars
dfp$tile_label <- format(dfp$median.fev, scientific = T, digits = 3)
# order variables
medv <- unlist(lapply(unique(dfp$variable), function(vi){
  median(dfp[dfp$variable == vi,]$median.fev, na.rm = T)}))
dfp$variable <- factor(dfp$variable, 
                       levels = unique(dfp$variable)[order(medv)])
pdf("ggtile_median-fev_gse-bias.pdf", 8, 3)
ggplot(dfp, aes(x = variable, y = fev.type)) +
  geom_tile(aes(fill = median.fev), color = "white") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2)) +
  geom_text(aes(label = tile_label), color = "white", size = 3) +
  scale_fill_gradient(low = "blue", high = "red", 
                      name = "Median\nFEV", 
                      labels = scales::label_scientific()) +
  xlab("Variable") + ylab("Model type")
dev.off()

# fev medians by tx
df.fev <- do.call(rbind, lapply(unique(mfev$tx), function(txi){
  mi <- mfev[mfev$tx == txi,]
  med.fevi <- apply(mi[,c(1:30)], 2, function(ci){
    median(as.numeric(ci), na.rm = T)})
  dfpi <- data.frame(label = names(med.fevi), 
                    median.fev = as.numeric(med.fevi),
                    stringsAsFactors = F)
  dfpi$tx <- txi; return(dfpi)}))
df.fev$fev.type <- gsub(".*_", "", df.fev$label)
df.fev$variable <- gsub("_.*", "", df.fev$label)

# heatmap, median fevs by tx
library(gridExtra)
pdf.fname <- "sfig_fev-bytx_gse-bias-results.pdf"
dfp <- df.fev
ltile <- lapply(unique(dfp$tx), function(txi){
  dfpi <- dfp[dfp$tx == txi,]
  # order vars
  medv <- unlist(lapply(unique(dfpi$variable), 
                        function(vari){
                          median(dfpi[dfpi$variable == vari,]$median.fev)}))
  dfpi$variable <- factor(dfpi$variable, levels = unique(dfpi$variable)[order(medv)])
  # get tile labels
  dfpi$labels <- format(dfpi$median.fev, scientific = T, digits = 3)
  ggtile <- ggplot(dfpi, aes(x=variable, y=fev.type)) +
    geom_tile(aes(fill = median.fev), color = "white") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_fill_gradient(low = "blue", high = "red", name = "Median\nFEV", 
                        labels = scales::label_scientific()) +
    ggtitle(txi) + geom_text(aes(label = labels), color = "white", size = 3)
  return(ggtile)})
# make composite plot
pdf(pdf.fname, 8, 8.3)
grid.arrange(ltile[[1]], ltile[[2]], ltile[[3]], ltile[[4]], 
             nrow = 4, bottom = "Variable", 
             left = "Model type")
dev.off()

#------------------------------------
# analyze the median abs diffs (MADs)
#------------------------------------
ldif <- get(load("list_gse-bias-diffs_all-tx.rda"))

# get mads by tx
mmad <- do.call(rbind, lapply(names(ldif), function(grp){
  lgrp.dif <- ldif[[grp]]
  varv <- unique(gsub("_.*", "", colnames(lgrp.dif)[1:30]))
  madv <- unlist(lapply(varv, function(vari){
    apply(lgrp.dif[,grepl(vari, colnames(lgrp.dif))], 2, 
          function(ci){median(abs(ci), na.rm = T)})
  }))
  mad.dfv <- data.frame(var = gsub("_.*", "", names(madv)),
                        diff = gsub(".*_", "", names(madv)),
                        tx = rep(grp, length(madv)),
                        mad = as.numeric(madv), stringsAsFactors = F)
  return(mad.dfv)}))

# get overall mads
mdif.all <- do.call(rbind, ldif)
varv <- unique(gsub("_.*", "", colnames(mdif.all)[1:30]))
madv <- unlist(lapply(varv, function(vari){
  apply(lgrp.dif[,grepl(vari, colnames(lgrp.dif))], 2, 
        function(ci){median(abs(ci), na.rm = T)})
}))
mad.dfv <- data.frame(var = gsub("_.*", "", names(madv)),
                      diff = gsub(".*_", "", names(madv)),
                      tx = rep("combined", length(madv)),
                      mad = as.numeric(madv), stringsAsFactors = F)
mmad <- rbind(mmad, mad.dfv)

# save table
mmad.fname <- "mads_gse-bias-sims"
save(mmad, file = paste0(mmad.fname,".rda"))
write.csv(mmad, file = paste0(mmad.fname,".csv"))
write.table(mmad, file = paste0(mmad.fname,".tsv"))

# get the heatmap of mads -- fig 2b
library(ggplot2); library(scales)

mmad <- get(load("mads_gse-bias-sims.rda"))
plot.fname <- "fig2b_ggtile_mads-labeled_gse-bias.pdf"
dfp <- mmad[mmad$tx == "combined",]
# order vars
medv <- unlist(lapply(unique(dfp$var), 
               function(vari){median(dfp[dfp$var == vari,]$mad)}))
dfp$var <- factor(dfp$var, levels = unique(dfp$var)[order(medv)])
# get tile labels
dfp$labels <- format(dfp$mad, scientific = T, digits = 3)
ggtile <- ggplot(dfp, aes(x=var, y=diff)) +
  geom_tile(aes(fill = mad), color = "white") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2)) +
  scale_fill_gradient(low = "blue", high = "red", 
                      name = "MAD", 
                      labels = scales::label_scientific()) +
  ylab("Difference type") + xlab("Variable") +
  geom_text(aes(label = labels), color = "white", size = 3)
pdf(plot.fname, 8, 3); print(ggtile);dev.off()

# get composite heatmaps of each subgroup type -- supp fig
library(gridExtra)
pdf.fname <- "sfig_mad-bytx_gse-bias-results.pdf"
dfp <- mmad[!mmad$tx == "combined",]
ltile <- lapply(unique(dfp$tx), function(txi){
  dfpi <- dfp[dfp$tx == txi,]
  # order vars
  medv <- unlist(lapply(unique(dfpi$var), 
                        function(vari){median(dfpi[dfpi$var == vari,]$mad)}))
  dfpi$var <- factor(dfpi$var, levels = unique(dfpi$var)[order(medv)])
  # get tile labels
  dfpi$labels <- format(dfpi$mad, scientific = T, digits = 3)
  ggtile <- ggplot(dfpi, aes(x=var, y=diff)) +
    geom_tile(aes(fill = mad), color = "white") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_fill_gradient(low = "blue", high = "red", name = "MAD", 
                        labels = scales::label_scientific()) +
    ggtitle(txi) +
    geom_text(aes(label = labels), color = "white", size = 3)
  return(ggtile)})
# make composite plot
pdf(pdf.fname, 8, 8.3)
grid.arrange(ltile[[1]], ltile[[2]], ltile[[3]], ltile[[4]], 
             nrow = 4, bottom = "Variable", 
             left = "Difference type")
dev.off()


#--------------------------------------------
# get the fev, diff median and variance stats
#--------------------------------------------
# helper functions
get_medians <- function(mdat, groupv){
  mmed <- do.call(rbind, lapply(groupv, function(typei){
    # get all merged medians
    col.str <- paste0(c(typei, "group"), collapse = "|")
    mff <- mdat[,c(grepl(col.str, colnames(mdat)))]
    statv <- apply(mff[,1:10], 2, median, na.rm = T)
    statv.format <- format(statv, scientific = T, digits = 3)
    med.merge <- c(statv.format, "merged_groups")
    # medians by sample group
    med.groups <- do.call(rbind, lapply(unique(mff$group), function(groupi){
      statv <- apply(mff[mff$group == groupi,1:10], 2, median, na.rm = T)
      statv.format <- format(statv, scientific = T, digits = 3)
      c(statv.format, groupi)}))
    med.all <- rbind(med.merge, med.groups)
    groupv <- med.all[,11]; med.all <- med.all[,c(1:10)]
    variablev <- gsub("_.*", "", colnames(med.all))
    typev <- gsub(".*_", "", colnames(med.all))
    med.all.format <- as.data.frame(t(med.all[,c(1:10)]), 
                                    stringsAsFactors = F)
    colnames(med.all.format) <- groupv
    med.all.format$variable <- varv
    med.all.format$model.type <- typev
    return(med.all.format)
  })); return(mmed)
}

get_variances <- function(mdat, groupv){
  mvar <- do.call(rbind, lapply(groupv, function(typei){
    # get all merged medians
    col.str <- paste0(c(typei, "group"), collapse = "|")
    mff <- mdat[,c(grepl(col.str, colnames(mdat)))]
    statv <- apply(mff[,1:10], 2, var, na.rm = T)
    statv.format <- format(statv, scientific = T, digits = 3)
    var.merge <- c(statv.format, "merged_groups")
    # medians by sample group
    var.groups <- do.call(rbind, lapply(unique(mff$group), function(groupi){
      statv <- apply(mff[mff$group == groupi,1:10], 2, var, na.rm = T)
      statv.format <- format(statv, scientific = T, digits = 3)
      c(statv.format, groupi)}))
    var.all <- rbind(var.merge, var.groups)
    groupv <- var.all[,11]; var.all <- var.all[,c(1:10)]
    variablev <- gsub("_.*", "", colnames(var.all))
    typev <- gsub(".*_", "", colnames(var.all))
    var.all.format <- as.data.frame(t(var.all[,c(1:10)]), stringsAsFactors = F)
    colnames(var.all.format) <- groupv; var.all.format$variable <- varv
    var.all.format$model.type <- typev; return(var.all.format)
  })); return(mvar)
}

get_mstat <- function(mdat, groupv){
  mmed <- get_medians(mdat, groupv); mvar <- get_variances(mdat, groupv)
  mmed$stat <- "median"; mvar$stat <- "variance"; mstat <- rbind(mmed, mvar)
  return(mstat)
}

# get matrices from lists
# fev matrix
mfsq <- do.call(rbind, lfmsq)
for(c in seq(ncol(mfsq))){mfsq[,c] <- as.numeric(mfsq[,c])}
mfsq$group <- gsub("\\..*", "", rownames(mfsq))
# diff matrix
mdif <- do.call(rbind, ldif)
mdif$group <- gsub("\\..*", "", rownames(mdif))

# get the mstat stables
mstat.fev <- get_mstat(mfsq, c("unadj", "adj1", "adj2"))
mstat.diff <- get_mstat(mdif, c("diff1", "diff2", "diff3"))

# save mstat stables
write.csv(mstat.fev, file = "st_mstat-fev.csv", row.names = F)
write.csv(mstat.diff, file = "st_mstat-diff.csv", row.names = F)

#-----------------
# results analyses
#-----------------
# summarize results
for(ii in seq(length(lmcg))){
  rmii <- lmcg[[ii]]; groupi <- names(lmcg)[ii]
  message(groupi, ": num unique cgs = ", length(unique(rmii[,94])))
  message(groupi, ": num unique studies = ", 
          length(unique(unlist(strsplit(rmii[,93], ";")))))
}
# cord-blood: num unique cgs = 12666
# cord-blood: num unique studies = 13
# whole-blood: num unique cgs = 14707
# whole-blood: num unique studies = 25
# pbmc: num unique cgs = 13269
# pbmc: num unique studies = 12
# all: num unique cgs = 14035
# all: num unique studies = 67

#------------------
# plot scatterplots
#------------------
# scatterplot, medians and variances
library(ggplot2); library(gridExtra); library(ggrepel); library(ggpubr)
groupv <- colnames(mstat)[4:8]
dfp <- do.call(rbind, lapply(groupv, function(groupi){
  dfp <- data.frame(med = mstat[mstat$stat == "median", groupi],
                    var = mstat[mstat$stat == "variance", groupi],
                    variable = mstat[mstat$stat == "variance", "variable"],
                    model_type = mstat[mstat$stat == "variance", "model.type"],
                    group = rep(groupi, nrow(mstat)),
                    stringsAsFactors = F);return(dfp)}))
for(c in 1:2){dfp[,c] <- as.numeric(dfp[,c])}
# make the full plots
ggfull.fname <- "ggpt-full_med-by-var_groups-2platforms.pdf"
lgg <- lapply(groupv, function(groupi){
  dfpi <- dfp[dfp$group == as.character(groupi),]
  dfpi$label_str <- paste0(dfpi$variable, "_", dfpi$model_type)
  dfpi <- dfpi[!duplicated(dfpi$label_str),]
  ggplot(dfpi, aes(x = med, y = var)) +
    geom_point(aes(color = variable, shape = model_type), alpha = 0.5) + theme_bw() +
    geom_text_repel(data = dfpi, aes(label = label_str, x = med, y = var),
                    segment.color = 'grey50', size = 4, force = 100, segment.size = 0.2, 
                    box.padding = 1, direction = "both") + ggtitle(groupi) + 
    theme(legend.position = "none")
}); names(lgg) <- groupv
plot <- ggplot(dfp, aes(x = med, y = var)) +
  geom_point(aes(color = variable, shape = model_type), alpha = 0.5) + 
  theme_bw()
lgg[["legend"]] <- as_ggplot(get_legend(plot))
pdf(ggfull.fname, 10, 10)
grid.arrange(lgg[[1]], lgg[[2]], lgg[[3]], lgg[[4]], lgg[[5]], lgg[[6]], 
             layout_matrix = matrix(seq(6), nrow = 2))
dev.off()
# make the zoomed plots
ggzoom.fname <- "ggpt-zoom_med-by-var_groups-2platforms.pdf"
xmax = 0.001; ymax = 0.01
lgg <- lapply(groupv, function(groupi){
  dfpi <- dfp[dfp$group == as.character(groupi),]
  dfpi$label_str <- paste0(dfpi$variable, "_", dfpi$model_type)
  dfpi <- dfpi[!duplicated(dfpi$label_str),]
  ggplot(dfpi, aes(x = med, y = var)) +
    geom_point(aes(color = variable, shape = model_type, size = 2), alpha = 0.5) + theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
    xlim(0, xmax) + ylim(0, ymax) + xlab("") + ylab("") + ggtitle(groupi)
    #geom_text_repel(data = dfpi, aes(label = label_str, x = med, y = var),
    #                segment.color = 'grey50', size = 3, force = 100, segment.size = 0.2, 
    #                box.padding = 0.5, direction = "both")
}); names(lgg) <- groupv
plot <- ggplot(dfp, aes(x = med, y = var)) +
  geom_point(aes(color = variable, shape = model_type), alpha = 0.5) + 
  theme_bw()
lgg[["legend"]] <- as_ggplot(get_legend(plot))
pdf(ggzoom.fname, 10, 9)
grid.arrange(lgg[[1]], lgg[[2]], lgg[[3]], lgg[[4]], lgg[[5]], lgg[[6]], 
             layout_matrix = matrix(seq(6), nrow = 2),
             bottom = "Median", left = "Variance")
dev.off()

#--------------
# plot heatmaps
#--------------

# mad heatmaps, fig 2b
mstat.mad <- mstat.diff
for(ci in c(1:5)){
  mstat.mad[,ci] <- median()
}



# make matrix for plotting


dfp <- do.call(rbind, lapply(c("unadj", "adj1", "adj2"), function(groupi){
  dfpi <- mmed[,groupi]; variablev <- gsub("\\_.*", "", names(dfpi))
  valv <- as.numeric(dfpi); groupv <- rep(groupi, length(valv))
  varv <- as.numeric(mvar[,groupi])
  return(data.frame(var = variablev, var = varv, 
                    value = valv, group = groupv))
})); dfp <- as.data.frame(dfp, stringsAsFactors = F)
dfp$group <- factor(dfp$group, levels = c("unadj", "adj1", "adj2"))
medv <- unlist(lapply(varv, function(x){mean(dfp[dfp$var==x,2])}))
dfp$var <- factor(dfp$var, levels = varv[order(medv)])
dfp$value.label <- format(dfp$value, scientific = T, digits = 3)
dfp$`Median\nFEV` <- dfp$value

# make heatmap for fig1
library(ggplot2)
pdf("ggtile.pdf", 4.5, 3)
ggtile <- ggplot(dfp, aes(x=group, y=var)) +
  geom_tile(aes(fill = `Median\nFEV`)) + theme_bw() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 2.8e-3) +
  geom_text(aes(label = value.label), color = "gray20") +
  xlab("Adjustment type") + ylab("Variable")
ggtile
dev.off()


# format for summary/plots
dfp <- do.call(rbind, lapply(seq(length(ldif)), function(ii){
  mdiff <- ldif[[ii]]
  mdiff$group <- names(ldif)[ii]
  return(mdiff)
}))
dim(dfp) # [1] 63829    32
# make the plot matrix
dsp <- do.call(rbind, lapply(c("diff1", "diff2", "diff3"), function(diffi){
  dfpf <- dfp[,grepl(diffi, colnames(dfp))]
  ds <- do.call(rbind, lapply(varv, function(variablei){
    which.variable <- which(grepl(variablei, colnames(dfpf)))
    medi <- median(dfpf[,which.variable], na.rm = T)
    vari <- var(dfpf[,which.variable], na.rm = T)
    c(variablei, medi, vari)}))
  ds <- as.data.frame(ds, stringsAsFactors = F)
  colnames(ds) <- c("variable", "median", "variance")
  ds$diff <- diffi
  return(ds)
}))
# make heatmap for fig1
library(ggplot2)
medv <- unlist(lapply(varv, function(x){median(dsp[dsp$variable==x,2])}))
dsp$variable <- factor(dsp$variable, levels = varv[order(medv)])
dsp$MAD <- as.numeric(dsp$median)
dsp$value.label <- format(dsp$MAD, scientific = T, digits = 3)
pdf("ggtile2.pdf", 4.5, 3)
ggtile <- ggplot(dsp, aes(x=diff, y=variable)) +
  geom_tile(aes(fill = MAD)) + theme_bw() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 1e-3) +
  geom_text(aes(label = value.label), color = "gray20") +
  xlab("Difference type") + ylab("Variable")
ggtile
dev.off()


df.fname <- "df-percmssq_gse-bias_blood-groups_2platforms.rda"
save(dfp, file = df.fname)

dim(dfp)
# [1] 328048     32




dfvar <- matrix(nrow = 0, ncol = 6)
for(vari in varv){
  dfpf <- dfp[,grepl(vari, colnames(dfp))]
  diff1.vari <- dfpf[,grepl("diff1", colnames(dfpf))]
  diff2.vari <- dfpf[,grepl("diff2", colnames(dfpf))]
  diff3.vari <- dfpf[,grepl("diff3", colnames(dfpf))]
  med.abs.diff1 <- format(median(abs(diff1.vari), na.rm = T),
                          scientific = T, digits = 3)
  med.abs.diff2 <- format(median(abs(diff2.vari), na.rm = T),
                          scientific = T, digits = 3)
  med.abs.diff3 <- format(median(abs(diff3.vari), na.rm = T),
                          scientific = T, digits = 3)
  num.models <- length(diff1.vari[!is.na(diff1.vari)])
  type <- "merge_all"
  dfvari <- matrix(c(vari, length(diff1.vari[!is.na(diff1.vari)]), 
                     med.abs.diff1, med.abs.diff2, med.abs.diff3, type),
                   nrow = 1)
  dfvar <- rbind(dfvar, dfvari)
  for(groupi in unique(dfp$group)){
    dfpfi <- dfp[dfp$group == groupi, grepl(vari, colnames(dfp))]
    diff1.vari <- dfpfi[,grepl("diff1", colnames(dfpfi))]
    diff2.vari <- dfpfi[,grepl("diff2", colnames(dfpfi))]
    diff3.vari <- dfpfi[,grepl("diff3", colnames(dfpfi))]
    med.abs.diff1 <- format(median(abs(diff1.vari), na.rm = T),
                            scientific = T, digits = 3)
    med.abs.diff2 <- format(median(abs(diff2.vari), na.rm = T),
                            scientific = T, digits = 3)
    med.abs.diff3 <- format(median(abs(diff3.vari), na.rm = T),
                            scientific = T, digits = 3)
    num.models <- length(diff1.vari[!is.na(diff1.vari)])
    type <- groupi
    dfvari <- matrix(c(vari, length(diff1.vari[!is.na(diff1.vari)]), 
                       med.abs.diff1, med.abs.diff2, med.abs.diff3, type),
                     nrow = 1)
    dfvar <- rbind(dfvar, dfvari)
  }
}
colnames(dfvar) <- c("varname", "num.models", "med.abs.diff1", 
                     "med.abs.diff2", "med.abs.diff3", "type")
dfvar

# get the magnitude mssq perc

rm(lmcg)
msq <- do.call(rbind, lfmsq)
rm(lfmsq)

for(c in seq(ncol(msq))){msq[,c] <- as.numeric(msq[,c])}
summary(msq)

for(c in seq(ncol(msq))){
  cname <- colnames(msq)[c]
  medc <- format(msq[,c], scientific = T, digits = 3)
  message(cname, ": ", median)
}


min(as.numeric(unlist(dfvar[,c(3:4)])))
max(as.numeric(unlist(dfvar[,c(3:4)])))

# make plots
library(ggplot2); library(ggdist)

dfpi <- dfp[!is.na(dfp$predage_diff3),]
dfpi$ngse <- droplevels(dfpi$ngse)

ggplot(dfpi, aes(x = ngse, y = predage_diff3, colour = group)) +
  geom_boxplot() + ylim(-1e-4,1e-4)
  
  geom_smooth(method = "loess", se = FALSE)


dfpi <- rbind(data.frame(ngse = dfp$ngse, diff = dfp$predage_diff1,
                         type = rep("diff1", nrow(dfp))),
              data.frame(ngse = dfp$ngse, diff = dfp$predage_diff2,
                         type = rep("diff2", nrow(dfp))))

ggplot(dfpi, aes(x = ngse, y = diff, colour = type, fill = type)) + ylim(0, 0.1) +
  stat_halfeye(width = 0.5, justification = -0.16, point_colour = NA, .width = 0)




dfpi <- rbind(data.frame(ngse = dfp$ngse, diff = dfp$predage_diff1,
                         type = rep("diff1", nrow(dfp))),
              data.frame(ngse = dfp$ngse, diff = dfp$predage_diff2,
                         type = rep("diff2", nrow(dfp))))

dfpi <- dfpi[!is.na(dfpi[,2]),]
dfpi[,1] <- droplevels(dfpi[,1])
dfpi[,2] <- as.numeric(dfpi[,2])

ggplot(dfpi, aes(x = ngse, y = diff, fill = type)) + geom_boxplot() +
  ylim(0,1)

dfpi <- rbind(data.frame(ngse = dfp$ngse, diff = dfp$predsex_diff1,
                         type = rep("diff1", nrow(dfp))),
              data.frame(ngse = dfp$ngse, diff = dfp$predsex_diff2,
                         type = rep("diff2", nrow(dfp))))
dfpi <- dfpi[!is.na(dfpi[,2]),]
dfpi[,1] <- droplevels(dfpi[,1])
dfpi[,2] <- as.numeric(dfpi[,2])

ggplot(dfpi, aes(x = ngse, y = diff, fill = type)) + geom_boxplot() +
  ylim(-0.01, 0.01)


plot(dfp[,1], dfp[,2])

ggplot(dfp, aes())


pdf("ggpt-gsebias-predage_diff1-diff2_2platforms-4groups.pdf", 5, 5)
ggpt <- ggplot(dfp, aes(x = predage_diff1, y = predage_diff2, color = group)) +
  geom_point(alpha = 0.5) + theme_bw()
ggpt
dev.off()

pdf("ggpt-gsebias_ngse-vs-var-diff1_2platforms-4groups.pdf", 5, 5)
ggpt <- ggplot() +
  geom_point(data = dfp, aes(x = ngse, y = predage_diff1, color = group), 
             alpha = 0.5) + theme_bw()
ggpt
dev.off()

# get dfp for stat smooths of vars
varv <- gsub("_.*", "", colnames(dfp)[1:30])
dfps <- do.call(rbind, lapply(varv, function(vari){
  dfpf <- dfp[,grepl(vari, colnames(dfp))]
  colnames(dfpf) <- c("diff1", "diff2", "diff3")
  dfpf$var <- vari; dfpf$group <- dfp$group
  dfpf$ngse <- dfp$ngse
  return(dfpf)
}))

dfp$ngse <- as.factor(dfp$ngse)
dfp <- dfp[!is.na(dfp$ngse),]
pdf("ggbox-gsebias-bygse_diff1-diff2-group_2platforms-4groups.pdf", 5, 5)
ggbox <- ggplot(data = dfp, aes(x = ngse, y = predage_diff1, fill = group)) +
  geom_boxplot(outlier.size = 0) + ylim(0,0.05) + theme_bw()
ggbox
dev.off()

dfp$log_predage_diff1 <- log10(dfp$predage_diff1)
pdf("ggviolin-gsebias-bygse_diff1-diff2-group_2platforms-4groups.pdf", 5, 5)
ggbox <- ggplot(data = dfp, aes(x = ngse, y = log_predage_diff1, fill = group)) +
  geom_violin(draw_quantiles = 0.5) + theme_bw()
ggbox
dev.off()

pdf("ggsmooth-gsebias_ngse-vs-var-diff1_2platforms-4groups.pdf", 5, 5)
ggsmooth <- ggplot() +
  stat_smooth(data = dfps, method = "loess", se = FALSE, group = var,
              aes(x = ngse, y = diff1)) + 
  theme_bw()
ggsmooth
dev.off()

pdf("ggpt-facet-gsebias_diff1-diff2_2platforms-4groups.pdf", 15, 5)
ggpt <- ggplot(dfp, aes(x = predage_diff1, y = predage_diff2, color = group)) +
  geom_point(alpha = 0.5) + theme_bw() + geom_abline(slope = 1, intercept = 0)
ggpt + facet_grid(cols = vars(group))
dev.off()

pdf("ggpt-facetzoom-gsebias_diff1-diff2_2platforms-4groups.pdf", 8, 2)
ggpt <- ggplot(dfp, aes(x = predage_diff1, y = predage_diff2, color = group)) +
  geom_point(alpha = 0.5) + theme_bw() + geom_abline(slope = 1, intercept = 0) +
  xlim(0, 10) + ylim(0,10) + xlab("Unadj. - Adj. 1") + ylab("Unadj. - Adj. 2")
ggpt + facet_grid(cols = vars(group))
dev.off()

dfp$ngse <- as.numeric(dfp$ngse)
pdf("ggviolin-gsebias_ngse-diff1_2platforms-4groups.pdf", 8, 2)
ggviolin <- ggplot(dfp, aes(x = ngse, y = predage_diff1, fill = group)) +
  geom_violin(draw_quantiles = 0.5) + theme_bw() + 
  xlab("Num. Studies") + ylab("Unadj. - Adj. 1")
# ggpt + facet_grid(cols = vars(group))
ggviolin
dev.off()
