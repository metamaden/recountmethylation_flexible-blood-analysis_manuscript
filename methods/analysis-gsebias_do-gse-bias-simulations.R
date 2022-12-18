#!/usr/bin/env R

# Author: Sean Maden
#
# Simulate the result of GSE bias corrections on variances, inc.
# a general GSE bias correction across available studies (adj1) and
# a precise GSE bias correction on a subset of studies evaluated (adj2).
# For comparison, the variances from unadjusted DNAm is also calculated.
# 
# Main script steps: 
# 1. Load the grset (autosomal probes noob-norm only, not GSE-adjusted).
# 3. Filter samples for a given blood subgroup.
# 4. Run the simulations, evaluations across iterations of probes and 
#     studies. Adjustments are performed in the same probes across 
#     randomized study sets.
#
#

library(HDF5Array); library(minfi)
library(methyPre); library(limma); library(sva)
library(data.table)

#--------------------
# load data -- server
#--------------------
save.dpath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_final_results")
# load noob-norm data (pre gse-adj)
gr.fname <- "gr-noob_h5se_hm450k-epic-merge_0-0-3"
gr <- loadHDF5SummarizedExperiment(file.path(save.dpath, gr.fname))

#----------------------------------
# helper functions, parallel method
#----------------------------------
# get anova data from model string
get_aovdat <- function(aov.str, labstr = "", pheno_subset, 
                       orderv = c("gse", "glint.epi.pc1", "glint.epi.pc2", "predage", 
                                  "predcell.CD8T", "predcell.CD4T", "predcell.NK", 
                                  "predcell.Bcell", "predcell.Mono", "predcell.Gran", 
                                  "Residuals", "predsex", "platform")){
  xaov <- eval(parse(text = aov.str)); xdat <- xaov[[1]]
  namev <- gsub(" ", "", rownames(xdat[1])); xdat.iter <- xdat[,c(2,4,5)]
  typev <- c("sumsq", "fval", "pval")
  vdat <- unlist(lapply(1:3, function(ii){
    datii <- xdat.iter[,ii]; names(datii) <- namev
    for(vname in c("predsex", "platform")){
      if(!vname %in% namev){
        datii <- c(datii, "NA"); names(datii)[length(datii)] <- vname}}
    datii <- datii[order(match(names(datii), orderv))]
    names(datii) <- paste0(names(datii), "_", typev[ii])
    return(datii)}))
  names(vdat) <- paste0(names(vdat), "_", labstr)
  return(vdat)
}

# parallel function for gse reps
par_gserep <- function(gse_rep, num.studies.subset, gsev, pheno, cgidv, table.fpath){
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
  cnv_numeric <- c("glint.epi.pc1", "glint.epi.pc2", "predage",
                   colnames(pheno_subset)[grepl("predcell", colnames(pheno_subset))])
  lm.str <- paste0("~ gse + ", paste0(cnv_numeric, collapse = " + "))
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
    return(mrep)}))
  mrep <- as.data.frame(mrep, stringsAsFactors = F)
  mrep$gse_rep <- gse_rep
  mrep$ngse <- num.studies.subset
  mrep$gsev <- paste0(unique(pheno_subset$gse), collapse = ";")
  mrep$cgid <- cgidv
  # write new results rows
  message("Writing new results rows to file ", table.fpath, "...")
  data.table::fwrite(mrep, file = table.fpath, sep = ",", row.names = F,
                     col.names = F, append = T)
  return(mrep)
}

# get all results for a series of study reps
get_mgserep <- function(ngse_rep, num.studies.subset, cgidv, pheno, table.fpath){
  gsev <- unique(pheno$gse) # get studies
  # format pheno colnames
  cnv_numeric <- c("predage", "glint.epi.pc1", "glint.epi.pc2", 
                   colnames(pheno)[grepl("predcell", colnames(pheno))])
  cnv_factor <- c("predsex", "platform", "gse")
  for(c in cnv_numeric){pheno[,c] <- as.numeric(pheno[,c])}
  for(c in cnv_factor){pheno[,c] <- as.factor(pheno[,c])}
  # get results matrix, process reps in parallel
  message("Processing studies using ", ngse_rep, " cores...")
  mgserep <- do.call(rbind, mclapply(seq(ngse_rep), par_gserep, 
                                     num.studies.subset = num.studies.subset,
                                     gsev = gsev, pheno = pheno, cgidv = cgidv, 
                                     mc.cores = ngse_rep, table.fpath = table.fpath))
  return(mgserep)
}

# get the results of all reps of random study selections
get_mgse_all <- function(grf, ngsev, table.fpath, num.gse = 5, 
                         ngse_rep = 10, num.probes = 500){
  message("Getting ",num.probes," random probes and ", 
          num.gse, " random studies...")
  cgidv <- sample(rownames(grf), num.probes)
  gsev <- sample(unique(grf$gse), num.gse)
  grff <- grf[cgidv, colnames(grf[,grf$gse %in% gsev])]
  bval_unadj <- getBeta(grff)
  mval_unadj <- logit2(bval_unadj)
  message("Getting adjustment 1 data...")
  mval_adj <- removeBatchEffect(mval_unadj, batch = grff$gse)
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
                cgidv = cgidv, pheno = pheno, table.fpath = table.fpath)}))
  message("Finished all ANOVA iterations; returning...")
  return(mgse.all)
}

# parallel wrapper for `get_mgse_all()`
par_mgse_all <- function(cgrep, grf, ngsev, table.fpath, num.gse = 5, 
                         ngse_rep = 10, num.probes = 500){
  t1 <- Sys.time(); message("Beginning cgrep ", cgrep, "...")
  mgse_all <- get_mgse_all(grf = grf, ngsev = ngsev, num.gse = num.gse, 
                           ngse_rep = ngse_rep, num.probes = num.probes,
                           table.fpath = table.fpath)
  mgse_all <- as.data.frame(mgse_all, stringsAsFactors = F)
  message("Finished cgrep ", cgrep, ", time: ", Sys.time() - t1)
  mgse_all$cgrep <- cgrep; return(mgse_all)
}

# script to process a subgroup
do_subgroup_tests <- function(grf, num.cgrep = 5, ngsev = c(2, 3, 4), num.gse = 6, sepsym = ",", 
                              ngse.rep = 5, num.probes = 100, num.batch = 30, seed = 1,
                              rda.fname = "mgse-adj-test_blood-group-all_2platforms.rda",
                              table.fname = "mgse-adj-test_blood-group-all_2platforms.csv",
                              save.dpath = file.path("home", "metamaden", "bioinfo_appnote", 
                                                     "manuscript_final_results")){
  set.seed(seed);t1 <- Sys.time()
  # write new table, then write new lines as processes finish
  cnv <- c("gse", "glint.epi.pc1", "glint.epi.pc2", "predage", "predcell.CD8T",
           "predcell.CD4T", "predcell.NK", "predcell.Bcell", "predcell.Mono",
           "predcell.Gran", "Residuals", "predsex", "platform")
  cnv <- c(paste0(cnv, "_sumsq"), paste0(cnv, "_fval"), paste0(cnv, "_pval"))
  cnv <- c(paste0(cnv, "_unadj"), paste0(cnv, "_adj"), paste0(cnv, "_adj2"))
  cnv <- c(cnv, "gse_rep", "ngse", "gsev", "cgid")
  table.fpath <- file.path(save.dpath, table.fname)
  mcnv <- matrix(nrow = 0, ncol = length(cnv));colnames(mcnv) <- cnv
  data.table::fwrite(mcnv, file = table.fpath, sep = sepsym, append = F, 
                     col.names = T, row.names = F)
  # do iterations, writing new lines as they complete
  for(batch in seq(num.batch)){
    message("Beginning batch ",batch, ", time: ", Sys.time() - t1)
    mclapply(seq(num.cgrep), par_mgse_all, grf = grf, ngsev = ngsev, 
             num.gse = num.gse, ngse_rep = ngse.rep, num.probes = num.probes,
             mc.cores = num.cgrep, table.fpath = table.fpath)
  }
  if(file.exists(table.fpath)){return(TRUE)
  } else{
    message("Couldn't find new table file at ", table.fpath)
    return(FALSE)
  }
  return(NULL)
}

#---------------------------------------
# process subgroups -- server, test runs
#---------------------------------------
# process all -- takes about 20hr to complete on remote server
new.fname <- "mcg-gsebias-final_blood-4stypes-2platforms"
rda.fname <- paste0(new.fname, ".rda"); table.fname <- paste0(new.fname, ".csv")
mcg.test <- do_subgroup_tests(grf = gr, num.cgrep = 3, ngsev = c(2,3,4), num.gse = 5,
                             ngse.rep = 3, num.probes = 500, num.batch = 20, seed = 1,
                             table.fname = table.fname, rda.fname = rda.fname)

# append colnames
csv.fname <- "mcg-gsebias-final_blood-4stypes-2platforms.csv"
csv.fpath <- file.path("home", "metamaden", "bioinfo_appnote", 
                       "manuscript_final_results")
csv.fpath <- file.path(csv.fpath, csv.fname)
csv <- fread(csv.fpath, sep = ",", header = F, data.table = F)
# get colnames
cnv <- c("gse", "glint.epi.pc1", "glint.epi.pc2", "predage", "predcell.CD8T",
         "predcell.CD4T", "predcell.NK", "predcell.Bcell", "predcell.Mono",
         "predcell.Gran", "Residuals", "predsex", "platform")
cnv <- c(paste0(cnv, "_sumsq"), paste0(cnv, "_fval"), paste0(cnv, "_pval"))
cnv <- c(paste0(cnv, "_unadj"), paste0(cnv, "_adj"), paste0(cnv, "_adj2"))
cnv <- c(cnv, "gse_rep", "ngse", "gsev", "cgid")
colnames(csv) <- cnv
# save new table
csv.fname <- "mcg-gsebias-final_blood-4stypes-2platforms.csv"
csv.fpath <- file.path(csv.fname)
fwrite(csv, file = csv.fpath, sep = ",", col.names = T, row.names = F)

#--------------------------------
# get variance fract by variables
#--------------------------------
save.dpath <- file.path("home", "metamaden", "bioinfo_appnote", 
                        "manuscript_final_results")
# append colnames
csv.fname <- "mcg-gsebias-final_blood-4stypes-2platforms.csv"
csv.fpath <- file.path(save.dpath, csv.fname)
csv <- fread(csv.fpath, sep = ",", header = T, data.table = F)

# get sum of sq vars by model type
typev <- c("unadj", "adj1", "adj2")
varv <- c("gse", "glint.epi.pc1", "glint.epi.pc2", "predage", "predcell.CD8T",
          "predcell.CD4T", "predcell.NK", "predcell.Bcell", "predcell.Mono",
          "predcell.Gran", "Residuals", "predsex", "platform")
cnames.mfsq <- paste0(varv, "_", rep(typev, each = length(varv))) # groups by probe, group
lindex <- list("unadj" = 1:13, "adj" = 40:52, "adj2" = 79:91) # model type var indices
# get the list of fract sumsq by group
mcg <- csv
mfsq <- t(apply(mcg, 1, function(dati){
  do.call(cbind, lapply(lindex, function(indexv){
    datii <- dati[indexv]
    tot.rsq <- sum(as.numeric(datii), na.rm = t)
    fract.indexv.rsq <- as.numeric(datii)/tot.rsq
    names(fract.indexv.rsq) <- names(datii)
    return(fract.indexv.rsq)
  }))
}))
colnames(mfsq) <- cnames.mfsq
mfsq <- as.data.frame(mfsq, stringsAsFactors = F)
mfsq$ngse <- mcg$ngse
mfsq$cgid <- mcg$cgid
mfsq$gsev <- mcg$gsev

# save mfsq matrix
msq.fname <- "msq-gse-bias_all-blood-2-platforms.rda"
msq.fpath <- file.path(save.dpath, msq.fname)
save(mfsq, file = msq.fpath)

#--------------------------------
# get variance diffs by variables
#--------------------------------
# get 3x FEV differences
cname.diff <- paste0(rep(varv, each = 3), "_", 
                     rep(c("diff1", "diff2", "diff3"), 
                         times = length(varv)))
mdiff.all <- do.call(cbind, lapply(varv, function(vari){
  fii <- mfsq[,grepl(vari, colnames(mfsq))];namev <- names(fii)
  diff1 <- as.numeric(fii[,1]) - as.numeric(fii[,2])
  diff2 <- as.numeric(fii[,1]) - as.numeric(fii[,3])
  diff3 <- as.numeric(fii[,2]) - as.numeric(fii[,3])
  mdiff <- cbind(diff1, cbind(diff2, diff3))
  return(mdiff)}))
mdiff.all <- as.data.frame(mdiff.all, stringsAsFactors = F)
colnames(mdiff.all) <- cname.diff
mdiff.all$ngse <- mfsq$ngse

# save mdiff.all matrix
mdiff.fname <- "mdiff-gse-bias_all-blood-2-platforms.rda"
mdiff.fpath <- file.path(save.dpath, mdiff.fname)
save(mdiff.all, file = mdiff.fpath)
message("done")