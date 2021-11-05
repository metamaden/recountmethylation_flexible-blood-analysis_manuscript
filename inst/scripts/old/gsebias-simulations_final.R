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
# combine arrays
gr.all <- combineArrays(gr.epic, gr.hm450k, 
                        outType = "IlluminaHumanMethylation450k")
dim(gr.all) # [1] 453093  65260

# load the detp tables
ptable1.path <- file.path(save.dpath, "detp-sstats_hm450k-blood-groups.rda")
ptable2.path <- file.path(save.dpath, "detp-sstats_epic-blood-groups.rda")
ptable1 <- get(load(ptable1.path))
ptable2 <- get(load(ptable2.path))

#--------------------------------
# filter gr.all, samples & probes
#--------------------------------
# filter samples
gr.all <- gr.all[,colnames(gr.all) %in% rownames(md)]
dim(gr.all) # [1] 453093  12858

# filter probes
anno <- getAnnotation(gr.all)
cg.keep <- rownames(anno[!anno$chr %in% c("chrY", "chrX"),])
data(chen_crxcg); cg.keep <- cg.keep[!cg.keep %in% chen.crxcg]
gr.all <- gr.all[rownames(gr.all) %in% cg.keep,]
dim(gr.all) # [1] 416045  12858

#-----------------
# helper functions
#-----------------
# test the impact of gse bias corrections in simulations
eval_adj <- function(grf, nrep_probes = 10, nrep_gse = 2, 
                     nprobes = 5, ngsev = c(5, 8, 10, 20, 30, 40)){
  set.seed(0); gsev <- unique(grf$gse)
  message("Defining pheno and format columns...")
  pheno.all <- colData(grf)
  cnv <- c("predage", colnames(pheno.all)[grepl("predcell", colnames(pheno.all))])
  for(cn in cnv){pheno.all[,cn] <- as.numeric(pheno.all[,cn])}
  message("Doing data reps...")
  mcgrep <- do.call(rbind, lapply(nrep_probes, function(cgrep){
    cgv.rep <- sample(rownames(grf), nprobes); grff <- grf[cgv.rep,]
    message("Getting DNAm across all probes, studies in rep...")
    bval.unadj <- getBeta(grff)
    # bval.unadj <- na.omit(bval.unadj) # handle na's
    mval.unadj <- logit2(bval.unadj)
    mval.adj <- removeBatchEffect(mval.unadj, batch = grff$gse)
    bval.adj <- ilogit2(mval.adj)
    message("Iterating on study counts...")
    mgserep <- do.call(rbind, lapply(ngsev, function(ngse){
      mgse <- do.call(rbind, lapply(seq(nrep_gse), function(nrep){
        message("Beginning rep ",nrep," with num gse ",ngse,"...")
        gsev.iter <- sample(gsev, ngse)
        pheno <- pheno.all[pheno.all$gse %in% gsev.iter,]
        message("Using studies: ", paste0(unique(pheno$gse), collapse = ","))
        which.gsm.iter <- which(colnames(bval.unadj) %in% rownames(pheno))
        bval.unadj.filt <- bval.unadj[,which.gsm.iter,drop = F]
        bval.adj.filt <- bval.adj[,which.gsm.iter,drop = F]
        cond <- identical(colnames(bval.unadj.filt), rownames(pheno))
        if(!cond){stop("Couldn't match pheno rownames to bval colnames.")}
        message("Getting the exact GSE adjusted Beta-values (adj2)...")
        mval.unadj.filt <- mval.unadj[,which.gsm.iter,drop=F]
        bval.adj2 <- removeBatchEffect(mval.unadj.filt, batch = pheno$gse)
        message("Appending Beta-values to pheno...")
        cgidv <- rownames(bval.unadj.filt)
        bval.unadj.filt.bind <- as.matrix(t(bval.unadj.filt))
        bval.adj.filt.bind <- as.matrix(t(bval.adj.filt))
        bval.adj2.filt.bind <- as.matrix(t(bval.adj2))
        colnames(bval.unadj.filt.bind) <- paste0(cgidv, "_unadj")
        colnames(bval.adj.filt.bind) <- paste0(cgidv, "_adj")
        colnames(bval.adj2.filt.bind) <- paste0(cgidv, "_adj2")
        pheno <- cbind(pheno, 
                       cbind(bval.unadj.filt.bind, 
                             cbind(bval.adj.filt.bind, 
                                   bval.adj2.filt.bind)))
        mres <- do.call(rbind, lapply(cgidv, function(cgid){
          message("Working on probe ", cgid, "...")
          lm.str <- paste0(" ~ gse + ", paste0(cnv, collapse = " + "))
          lm.str.unadj <- paste0(cgid, "_unadj", lm.str)
          lm.str.adj <- paste0(cgid, "_adj", lm.str)
          lm.str.adj2 <- paste0(cgid, "_adj2", lm.str)
          if(length(unique(pheno$predsex)) > 1){
            lm.str.unadj <- paste0(lm.str.unadj, " + predsex")
            lm.str.adj <- paste0(lm.str.adj, " + predsex")
            lm.str.adj2 <- paste0(lm.str.adj2, " + predsex")}
          aov.str.unadj <- paste0("summary(aov(", lm.str.unadj, ", data = pheno))")
          aov.str.adj <- paste0("summary(aov(", lm.str.adj, ", data = pheno))")
          aov.str.adj2 <- paste0("summary(aov(", lm.str.adj2, ", data = pheno))")
          message(cgid,": getting ANOVA results...")
          aov.unadj <- eval(parse(text = aov.str.unadj))
          aov.adj <- eval(parse(text = aov.str.adj))
          aov.adj2 <- eval(parse(text = aov.str.adj2))
          unadjv <- unlist(aov.unadj[[1]][,3:5])
          adjv <- unlist(aov.adj[[1]][,3:5])
          adj2v <- unlist(aov.adj2[[1]][,3:5])
          if(!length(unique(pheno$predsex)) > 1){
            unadjv <- c(unadjv[1:8], "NA", unadjv[9],
                        unadjv[10:18], "NA", unadjv[19],
                        unadjv[20:28], "NA", unadjv[29])
            adjv <- c(adjv[1:8], "NA", adjv[9],
                      adjv[10:18], "NA", adjv[19],
                      adjv[20:28], "NA", adjv[29])
            adj2v <- c(adj2v[1:8], "NA", adj2v[9],
                       adj2v[10:18], "NA", adj2v[19],
                       adj2v[20:28], "NA", adj2v[29])
          }
          dati <- c(unadjv, adjv, adj2v)
          mres <- matrix(c(cgrep, nrep, cgid, ngse, 
                           paste0(unique(pheno$gse), collapse = ","), 
                           dati), 
                         nrow = 1); return(mres)
        }))
        message("Formatting results...")
        vnames <- c("gse", "predage", "predcell.CD8T", 
                    "predcell.CD4T", "predcell.NK", "predcell.Bcell", 
                    "predcell.Mono", "predcell.Gran", "predsex", "Residuals")
        vtype <- c("mean-sq", "fval", "pval")
        label <- c(rep("bval-unadj", length(vnames)*length(vtype)),
                   rep("bval-adj1", length(vnames)*length(vtype)),
                   rep("bval-adj2", length(vnames)*length(vtype)))
        vartype <- rep(vtype, each = length(vnames))
        labelsv <- gsub(" ", "_", paste(vnames, vartype, label))
        cnames <- c("cgrep", "gserep", "cgid", "ngse", "gsev", labelsv)
        colnames(mres) <- cnames
        return(mres)
      }))
      message("Finished gse rep ", nrep, "!"); return(mgse)}))
    message("Finished mgserep ", ngse, "!")
    return(mgserep)
  }))
  message("Finished cgrep ", cgrep, "!")
  return(mcgrep)
}

#----------------------------------------------------
# run bias correction experiments across blood groups
#----------------------------------------------------
set.seed(0)

# all blood
# filter detp 
cg.keep <- unique(ptable1[ptable1$perc_above_01.all == 0,1],
  ptable2[ptable2$perc_above_01.all == 0,1])
grf <- gr.all[rownames(gr.all) %in% cg.keep,]
dim(grf) # [1] 286798  12858
groupi <- groupi.label <- "all"
save.fname <- paste0("mgse-adj-test_blood-group-",
                     groupi.label, "_2platforms.rda")
save.fpath <- file.path(save.dpath, save.fname)
which.sample <- which(md$blood_subgroup == groupi)
length(unique(grf$gse)) # 66
rm(gr.all); rm(gr.hm450k); rm(epic)
mgse <- eval_adj(grf, nrep_probes = 10, nrep_gse = 5,
                 nprobes = 1000, ngsev = c(5, 10, 30, 50))
save(mgse, file = save.fpath)

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
dim(grf) # [1] 316700  12858
# get the groupi sample ids
gsm.labelv <- rownames(md[md$blood_subgroup == groupi,])
grf <- grf[,colnames(grf) %in% gsm.labelv]
dim(grf) # [1] 316700   1417
length(unique(grf$gse)) # 9
mgse <- eval_adj(grf, nrep_probes = 10, nrep_gse = 5,
                 nprobes = 1000, ngsev = c(2, 4, 6, 8))
save(mgse, file = save.fpath)

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
set.seed(0)
bval <- as.numeric(getBeta(grf[,1]))
cgidv.na <- rownames(grf)[which(is.na(bval))]
length(cgidv.na) # [1] 545
grf <- grf[!rownames(grf) %in% cgidv.na,]
dim(grf)
# get the groupi sample ids
gsm.labelv <- rownames(md[md$blood_subgroup == groupi,])
grf <- grf[,colnames(grf) %in% gsm.labelv]
dim(grf) # [1] 416045   5980
length(unique(grf$gse)) # 30
mgse <- eval_adj(grf, nrep_probes = 10, nrep_gse = 5,
                 nprobes = 1000, ngsev = c(2, 5, 10, 22))
save(mgse, file = save.fpath)

# pbmc
# get file info
groupi <- "peripheral_blood_mononuclear_cells"; groupi.label <- "pbmc"
save.fname <- paste0("mgse-adj-test_blood-group-",
                     groupi.label, "_2platforms.rda")
save.fpath <- file.path(save.dpath, save.fname)
# get detp filtered cgs
cg.keep <- unique(ptable1[ptable1$perc_above_01.peripheral_blood_mononuclear_cells == 0,1],
                  ptable2[ptable2$perc_above_01.peripheral_blood_mononuclear_cells == 0,1])
grf <- gr.all[rownames(gr.all) %in% cg.keep,]
dim(grf) # [1] 399968  12858
# get the groupi sample ids
gsm.labelv <- rownames(md[md$blood_subgroup == groupi,])
grf <- grf[,colnames(grf) %in% gsm.labelv]
dim(grf) # [1] 399968    642
length(unique(grf$gse)) # 9
mgse <- eval_adj(grf, nrep_probes = 10, nrep_gse = 5,
                 nprobes = 1000, ngsev = c(2, 4, 6, 8))
save(mgse, file = save.fpath)
