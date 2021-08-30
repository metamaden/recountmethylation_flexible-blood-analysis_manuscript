#!/usr/bin/env R

# Author: Sean Maden
# 
# Perform ANOVA-based filter for predicted sex, age, and blood cell type 
# fractions. Uses parallelization to expedite ANOVA calculations across sample 
# subgroups and platforms.
#
# This is a *fast* implementation of aov() for ANOVAs on blood subgroups. The
# gist is a single bind combines chunk probe Beta-values data with the metadata, 
# then a vector of ANOVA formula strings is passed to aov() with vectorization 
# using lapply(). The data.frame lookup is much faster than binding individual 
# columns inside of lapply().

library(HDF5Array); library(purrr)

#-----------------
# helper functions
#-----------------
# parallelize function
dopar_adf <- function(sii, mbval, md, t1, max.na.perc = 0.2){
    # Calculates a series of anova tests, to be run in parallel.
    # mbval : A subsetted matrix of beta-values.
    # md : Sample metadata, where rows are the same as the mbval columns.
    # max.na.perc : Max allowed percentage of samples with NA values.
    message("Beginning at index ", min(sii), " to ", max(sii))
    mbvali <- mbval[sii,]; cgidv <- rownames(mbvali)
    # filter max na values
    max.na <- round(max.na.perc * ncol(mbvali), digits = 0)
    num.na <- apply(mbvali,1,function(x){length(x[is.na(x)])})
    mbvali <- mbvali[num.na <= max.na,]
    message("NA filter removed ", length(sii) - nrow(mbvali)," probes.")
    # store NA results
    mna <- matrix(c(cgidv, num.na), ncol = 2)
    # impute medians for na values
    message("Imputing NAs with row medians...")
    mbvali <- t(apply(mbvali, 1, function(x){
        x[is.na(x)] <- median(x, na.rm = T); return(x)}))
    # bind data.all
    dat.all <- cbind(md, t(as.matrix(mbvali)))
    formula.str <- paste0("~ ", paste0(colnames(mdff), collapse = " + "))
    formulae <- lapply(cgidv, function(x) as.formula(paste0(x, formula.str)))
    # do ANOVAs
    res <- lapply(seq(length(formulae)), function(x){
        message("start index: ", min(sii)," probe number: ", x)
        fx <- formulae[[x]]
        aovi <- try(aov(fx, data = dat.all))
        if(class(aovi)[1] == "aov"){return(summary(aovi))}else{
            return(aovi)}})
    names(res) <- cgidv
    message("Finished at index ", min(sii), " to ", max(sii), 
        ", time elapsed: ", Sys.time() - t1)
    lavi <- list(mna = mna, tav = res)
    return(lavi)
}

# dopar blood subgroups
dopar_subgroups <- function(lfv, mdf, fnstem = "anova-df_hm450k_blood-",
    num.sessions = 100, num.cores = 20, save.fpath = ".", 
    groupv = c("cord_blood", "whole_blood", "all",
        "peripheral_blood_mononuclear_cells")){
    # Manage a series of parallel sessions, across sample subgroups, using 
    # dopar_adf() calls.
    # lfv : list of paths to h5se files containing the sample groups to analyze.
    # mdf : Filtered metadata
    # fnstem : Stem of the new results objects to save
    # num.sessions : Total number of sessions/cores to use
    # save.path : Path to save new results object
    # groupv : Vector of sample groups to analyze
    message("Beginning to work on subgroups...")
    for(group in groupv){
        message("Beginning group ", group, "...")
        # load the h5se data
        gr.fname <- lfv[grepl(paste0(".*group-", group), lfv)]
        gri <- HDF5Array::loadHDF5SummarizedExperiment(gr.fname)
        # filter metadata
        mdff <- mdf[rownames(mdf) %in% colnames(gri),]
        mdff <- mdff[order(match(rownames(mdff), colnames(gri))),]
        cond <- identical(rownames(mdff), colnames(gri)); cond
        if(!cond){stop("Couldn't match metadata to gri colnames.")}
        # remove sex if 1 level only
        if(length(unique(mdff$predsex)) == 1){
            mdff <- mdff[,!colnames(mdff)=="predsex"]}
        # retain only autosomal cpgs for anova tests
        cg.anno <- getAnnotation(gri)
        cgv.filt <- rownames(cg.anno[!cg.anno$chr %in% c("chrY", "chrX"),])
        mbval <- getBeta(gri[rownames(gri) %in% cgv.filt,])
        # get all anova results, with parallelization
        message("Beginning ",num.sessions," processes across ",num.cores," cores...")
        ladf <- parallel::mclapply(splitIndices(nrow(mbval), num.sessions), 
            FUN = dopar_adf, mbval = mbval, md = mdff, t1 = Sys.time(), 
            mc.cores = num.cores)
        save.fpath.anova <- file.path(save.fpath, paste0(fnstem, group, ".rda"))
        message("Saving result to ", save.fpath.anova)
        save(ladf, file = save.fpath.anova)
        message("Finished with group ", group, "...")
    }; message("Finished all subgroups. Returning..."); return(NULL)
}

#----------
# load data
#----------
lfv <- list.files()
lfv.hm450k <- lfv[grepl("^gr-adj_hm450k.*", lfv)]
lfv.epic <- lfv[grepl("^gr-adj_epic.*", lfv)]

# hm450k metadata
md.fname <- "md-blood-filt_hm450k-merged_0-0-3.rda"
md <- get(load(md.fname)); mdf <- md[,c(2, 10:17)]
for(c in c(1,3)){mdf[,c] <- as.factor(mdf[,c])}
for(c in c(2, 4:9)){mdf[,c] <- as.numeric(mdf[,c])}
md.hm450k <- mdf

# epic metadata
# md.fname <- "md-blood-filt_epic-merged_0-0-3.rda"
md.fname <- "si2_blood-md-2platforms.rda"
md <- get(load(md.fname)); mdf <- md[md$platform == "epic",]
mdf <- mdf[,c(2, 10:17)]
for(c in c(1,3)){mdf[,c] <- as.factor(mdf[,c])}
for(c in c(2, 4:9)){mdf[,c] <- as.numeric(mdf[,c])}
md.epic <- mdf

# file paths
lmav.fpath <- file.path("lanova-matrix-results_2platforms_blood-all-groups.rda")
lmavf.fpath <- file.path(paste0("lanova-matrix-results-filt_",
    "max-xssq-20perc-max-pval-1eNeg5_2platforms_blood-all-groups.rda"))
lcg.fpath <- file.path(paste0("lcg-anovaf_max-xssq-20perc-max-pval-1eNeg5",
    "_2platforms_blood-all-groups"))
lcgfinal.fpath <- paste0("list-cgv-final-filt_anova-and-detp001_",
    "2platforms_blood-all-groups.rda")
lgrp.stat.fpath <- "lstat-bv-cgoverlap_2platforms_blood-groups.rda"

#----------------------
# subgroups anova tests
#----------------------
# Use parallelization to caluculate ANOVAs for autosomal probes across
# sample subgroups and platforms.
platforms <- c("hm450k", "epic")
anova_fnstems <- c("anova-df_hm450k_blood-", "anova-df_epic_blood-")
llfv <- list("hm450k" = lfv.hm450k, "epic" = lfv.epic)
lmd <- list("hm450k" = md.hm450k, "epic" = md.epic)
for(pi in seq(length(platforms))){
    platform <- platforms[pi]; message("Beginning platform ", platform, "...")
    dopar_subgroups(lfv = llfv[[platform]], mdf = lmd[[platform]], 
        fnstem = anova_fnstems[pi], num.sessions = 1000, 
        num.cores = 20, save.fpath = ".")}

#------------------------------
# read results into flat matrix
#------------------------------
# get 1 matrix per subgroup
lfv <- list.files(); lfv <- lfv[grepl("anova-df", lfv)]; 
lmav <- list(); t1 <- Sys.time()
for(fn in lfv){
    message("Beginning with file ", fn); lavi <- get(load(fn))
    # unlist probe data by 2 levels
    lcgi <- unlist(unlist(lapply(lavi, function(x){x[[2]]}), 
        recursive = F), recursive = F)
    vvar <- gsub(" ", "", rownames(lcgi[[1]]))
    # map anova data
    mpval <- do.call(rbind, map(lcgi, function(x){
        res <- "NA"; if(length(x)>1){res <- x[[5]]}; return(res)}))
    colnames(mpval) <- paste0(vvar, ".pval")
    mxssq <- do.call(rbind, map(lcgi, function(x){
        res <- "NA"; if(length(x)>1){res <- x[[3]]}; return(res)}))
    colnames(mxssq) <- paste0(vvar, ".xssq")
    mav <- cbind(matrix(c(names(lcgi)), ncol = 1),
        cbind(mpval, mxssq))
    colnames(mav) <- c("cgid", colnames(mav)[2:ncol(mav)])
    mav <- as.data.frame(mav, stringsAsFactors = F)
    for(c in seq(2, ncol(mav), 1)){mav[,c] <- as.numeric(mav[,c])}
    lmav[[fn]] <- mav; message("Finished file ", lfi, 
        ", time elapsed: ", Sys.time() - t1)}

# save lmav
lmav.fpath <- file.path("lanova-matrix-results_2platforms_blood-all-groups.rda")
save(lmav, file = lmav.fpath)

#---------------------------
# evaluate the anova results
#---------------------------
lmav <- get(load(lmav.fpath))

# define filters
max.xssq <- 0.2; min.pval <- 1e-3

which.mavi <- 1
mavi <- lmav[[which.mavi]]

num.pval <- apply(mavi[,c(2:10)], 1, function(x){length(x[x <= 1e-3])})
summary(num.pval)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.000   3.000   5.000   4.315   6.000   9.000

num.pval <- apply(mavi[,c(2:10)], 1, function(x){length(x[x <= 1e-4])})
# summary(num.pval)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.000   2.000   4.000   3.847   5.000   9.000

num.xssq <- apply(mavi[,c(12:21)], 1, function(x){
    xdenom <- sum(x); xperc <- x/xdenom; length(xperc[xperc >= 0.1])})
summary(num.xssq)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   1.000   1.000   2.000   2.174   3.000  10.000

num.xssq <- apply(mavi[,c(12:21)], 1, function(x){
    xdenom <- sum(x); xperc <- x/xdenom; length(xperc[xperc >= 0.2])})
# summary(num.xssq)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.000   1.000   1.000   1.466   2.000  10.000

# get probes after filter
max.xssq <- 0.2; min.pval <- 1e-5; lmavf <- list()
for(i in seq(length(lmav))){
    mavi <- lmav[[i]]; namei <- names(lmav)[i]
    # get probes with high var contributions
    num.xssq <- apply(mavi[,c(12:20)], 1, function(x){
        xperc <- x/sum(x); length(xperc[xperc >= max.xssq])})
    which.cg.xssq <- which(num.xssq > 0)
    # get subset with significant contributions
    mavi.filt <- mavi[which.cg.xssq,]
    num.pval <- apply(mavi.filt[,c(2:10)], 1, function(x){
        length(x[x <= min.pval])})
    cgidv.rm <- mavi.filt[which(num.pval > 0),1]
    # get the filtered matrix
    mavi.final <- mavi.filt[!mavi.filt[,1] %in% cgidv.rm,]
    lmavf[[namei]] <- mavi.final
    message(namei); message(dim(mavi.final))
}
# anova-df_epic_blood-all.rda
# 37820 21
# anova-df_epic_blood-cord_blood.rda
# 226989 21
# anova-df_epic_blood-peripheral_blood_mononuclear_cells.rda
# 459458 21
# anova-df_epic_blood-whole_blood.rda
# 492483 21
# anova-df_hm450k_blood-all.rda
# 36778 21
# anova-df_hm450k_blood-cord_blood.rda
# 132219 21
# anova-df_hm450k_blood-peripheral_blood_mononuclear_cells.rda
# 181945 21
# anova-df_hm450k_blood-whole_blood.rda
# 38932 21

# anova-df_epic_blood-all.rda
# 37820 21
# anova-df_hm450k_blood-all.rda
# 36778 21

# anova-df_epic_blood-cord_blood.rda
# 226989 21
# anova-df_hm450k_blood-cord_blood.rda
# 132219 21

# anova-df_hm450k_blood-peripheral_blood_mononuclear_cells.rda
# 181945 21
# anova-df_epic_blood-peripheral_blood_mononuclear_cells.rda
# 459458 21

# anova-df_epic_blood-whole_blood.rda
# 492483 21
# anova-df_hm450k_blood-whole_blood.rda
# 38932 21


# save lmavf
lmavf.fpath <- file.path(paste0("lanova-matrix-results-filt_",
    "max-xssq-20perc-max-pval-1eNeg5_2platforms_blood-all-groups.rda"))
save(lmavf, file = lmavf.fpath)

# save filtered probe ids
lcg.fpath <- file.path(paste0("lcg-anovaf_max-xssq-20perc-max-pval-1eNeg5",
    "_2platforms_blood-all-groups"))
lcg.avf <- lapply(lmavf, function(x){x[,1]})
names(lcg.avf) <- names(lmavf); save(lcg.avf, file = lcg.fpath)

#-------------------
# evaluate NA values
#-------------------
for(fn in lfv){
    message(fn); lavi <- get(load(fn))
    mna <- do.call(rbind, lapply(lavi, function(x){x$mna}))
    print(summary(as.numeric(mna[,2])))
    print(nrow(mna[as.numeric(mna[,2]) > 10,]))
    print(nrow(mna[as.numeric(mna[,2]) > 50,]))
    print(nrow(mna[as.numeric(mna[,2]) > 100,]))
    print(nrow(mna[as.numeric(mna[,2]) > 500,]))
    print(nrow(mna[as.numeric(mna[,2]) > 1000,]))
}   
# anova-df_epic_blood-all.rda
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#    0.000    0.000    0.000    1.791    0.000 2117.000
# [1] 873
# [1] 873
# [1] 873
# [1] 657
# [1] 657
# anova-df_epic_blood-cord_blood.rda
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   0.0000   0.0000   0.0000   0.4519   0.0000 539.0000
# [1] 923
# [1] 696
# [1] 696
# [1] 696
# [1] 0
# anova-df_epic_blood-peripheral_blood_mononuclear_cells.rda
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   0.0000   0.0000   0.0000   0.1167   0.0000 142.0000
# [1] 696
# [1] 696
# [1] 696
# [1] 0
# [1] 0
# anova-df_epic_blood-whole_blood.rda
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#   0.0000   0.0000   0.0000   0.2027   0.0000 223.0000
# [1] 923
# [1] 923
# [1] 923
# [1] 0
# [1] 0
# anova-df_hm450k_blood-all.rda
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#       0       0       0       0       0       0
# [1] 0
# [1] 0
# [1] 0
# [1] 0
# [1] 0
# anova-df_hm450k_blood-cord_blood.rda
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      0       0       0       0       0       0
# [1] 0
# [1] 0
# [1] 0
# [1] 0
# [1] 0
# anova-df_hm450k_blood-peripheral_blood_mononuclear_cells.rda
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#       0       0       0       0       0       0
# [1] 0
# [1] 0
# [1] 0
# [1] 0
# [1] 0
# anova-df_hm450k_blood-whole_blood.rda
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#       0       0       0       0       0       0
# [1] 0
# [1] 0
# [1] 0
# [1] 0
# [1] 0

#------------------------------------------------------------------------
# apply detp filter to anova-filtered probes, summarize platform overlaps
#------------------------------------------------------------------------
# load anova-filtered probes
lcg.fpath <- file.path(paste0("lcg-anovaf_max-xssq-20perc-max-pval-1eNeg5",
    "_2platforms_blood-all-groups"))
lcg <- get(load(lcg.fpath))
# load the detp stats
detp.epic <- get(load("detp-sstats_epic-blood-groups.rda"))
detp.hm450k <- get(load("detp-sstats_hm450k-blood-groups.rda"))

# get the epic and hm450k detp stats for the same group
# note: filter on less than 10% of samples having detp >= 1e-3
fnv <- names(lcg); which.col.filt <- 7; perc.thresh <- 0.1; lcgfinal <- list()
for(groupi in c("all", "cord_blood", "peripheral_blood_mononuclear_cells", 
    "whole_blood")){ # groupi <- "all"
    cond.groupi<- grepl(paste0(".*blood-",groupi,"\\.rda$"), fnv)
    fnvi.epic <- fnv[cond.groupi & grepl(".*_epic_.*", fnv)]
    fnvi.hm450k <- fnv[cond.groupi & grepl(".*_hm450k_.*", fnv)]
    cgv.epic <- lcg[[fnvi.epic]]; cgv.hm450k <- lcg[[fnvi.hm450k]]
    dpf.cond.str <- paste0("\\.", groupi, "$")
    dpf.cond <- grepl(dpf.cond.str, colnames(detp.epic))
    dpf.epic <- detp.epic[detp.epic[,1] %in% cgv.epic, dpf.cond]
    dpf.cond <- grepl(dpf.cond.str, colnames(detp.hm450k))
    dpf.hm450k <- detp.hm450k[detp.hm450k[,1] %in% cgv.hm450k, dpf.cond]
    # apply the detp filter
    dpff.epic <- dpf.epic[dpf.epic[,which.col.filt]<perc.thresh,]
    dpff.hm450k <- dpf.hm450k[dpf.hm450k[,which.col.filt]<perc.thresh,]
    # get overlapping probes
    cgf.ol <- intersect(dpff.epic[,1], dpff.hm450k[,1])
    dpff.epic.ol <- dpff.epic[dpff.epic[,1] %in% cgf.ol,]
    dpff.hm450k.ol <- dpff.hm450k[dpff.hm450k[,1] %in% cgf.ol,]
    epic.reorder <- order(match(dpff.epic.ol[,1], dpff.hm450k.ol[,1]))
    dpff.epic.ol <- dpff.epic.ol[epic.reorder,]
    identical(dpff.epic.ol[,1], dpff.hm450k.ol[,1])
    ct.median <- cor.test(dpff.epic.ol[,2], dpff.hm450k.ol[,2], 
        method = "spearman")
    ct.mean <- cor.test(dpff.epic.ol[,3], dpff.hm450k.ol[,3], 
        method = "spearman")
    ct.sd <- cor.test(dpff.epic.ol[,4], dpff.hm450k.ol[,4], 
        method = "spearman")
    # print results
    print(groupi); print(dim(dpff.epic)); print(dim(dpff.hm450k))
    print(length(intersect(dpff.epic[,1], dpff.hm450k[,1])))
    message("spearman, median: ", ct.median$estimate)
    message("spearman, mean: ", ct.mean$estimate)
    message("spearman, sd: ", ct.sd$estimate)
    lcgfinal[[groupi]] <- list("epic" = dpff.epic[,1], 
        "hm450k" = dpff.hm450k[,1], "overlap" = cgf.ol)
}
# [1] "all"
# [1] 34576     7
# [1] 34790     7
# [1] 6425
# spearman, median: 0.299666349871948
# spearman, mean: 0.41405611228793
# spearman, sd: 0.364256275256728
# [1] "cord_blood"
# [1] 215578      7
# [1] 120273      7
# [1] 75661
# spearman, median: 0.33713181553039
# spearman, mean: 0.29466025568977
# spearman, sd: 0.236195837357384
# [1] "peripheral_blood_mononuclear_cells"
# [1] 447186      7
# [1] 175388      7
# [1] 144580
# spearman, median: 0.455258207919038
# spearman, mean: 0.560128498796311
# spearman, sd: 0.535737760488615
# [1] "whole_blood"
# [1] 477950      7
# [1] 37232     7
# [1] 33852
# spearman, median: 0.419865681375051
# spearman, mean: 0.366033238844855
# spearman, sd: 0.320251811328634

# save the final filtered probe sets
lcgfinal.fpath <- paste0("list-cgv-final-filt_anova-and-detp001_",
    "2platforms_blood-all-groups.rda")
save(lcgfinal, file = lcgfinal.fpath)

"list-cgv-final-filt_anova-and-detp001_2platforms_blood-all-groups.rda"

#----------------------------
# get filt probes means, vars
#----------------------------
lcgfinal <- get(load(lcgfinal.fpath)); lgrp.stat <- list()
for(groupi in c("cord_blood", "whole_blood", "all",
        "peripheral_blood_mononuclear_cells")){
    message("subgroup: ", groupi); lcgi <- lcgfinal[[groupi]]
    gri.hm450k.fname <- paste0("gr-adj_hm450k_blood-group-",groupi,".rda")
    gri.epic.fname <- paste0("gr-adj_epic_blood-group-",groupi,".rda")
    gri.hm450k <- HDF5Array::loadHDF5SummarizedExperiment(gri.hm450k.fname)
    gri.epic <- HDF5Array::loadHDF5SummarizedExperiment(gri.epic.fname)
    message("Getting the Beta-values as matrices...")
    grf.hm450k <- gri.hm450k[rownames(gri.hm450k) %in% lcgi$hm450k,]
    grf.epic <- gri.epic[rownames(gri.epic) %in% lcgi$epic,]
    bv.hm450k <- as.matrix(getBeta(grf.hm450k))
    bv.epic <- as.matrix(getBeta(grf.epic))
    message("Getting the overlapping probe Beta-values as a matrix...")
    bvf.hm450k <- bv.hm450k[lcgi$overlap,];bvf.epic <- bv.epic[lcgi$overlap,]
    order.cgv <- order(match(rownames(bvf.hm450k), rownames(bv.epic)))
    bvf.hm450k <- bvf.hm450k[order.cgv,] 
    cond <- identical(rownames(bvf.hm450k), rownames(bvf.epic))
    if(!cond){stop("Couldn't match bvf.hm450k and bvf.epic ",
        "on overlapping probes.")}; bv.shared <- cbind(bvf.hm450k, bvf.epic)
    lbv <- list("hm450k" = bv.hm450k, "epic" = bv.epic, "shared" = bv.shared)
    message("Calculating summary statistics...")
    lstat <- lapply(lbv, function(x){list("mean" = rowMeans(x, na.rm =T),
        "var" = matrixStats::rowVars(x, na.rm = T))})
    names(lstat) <- c("hm450k", "epic", "shared")
    lgrp.stat[[groupi]] <- lstat; message("Finished group ", groupi)}
lgrp.stat.fpath <- "lstat-cgfilt_bval-mean-var_2platforms_blood-groups.rda"
save(lgrp.stat, file = lgrp.stat.fpath)

#--------------------------------------------
# compare overlapping probes across platforms
#--------------------------------------------
lgrp.stat <- get(load(lgrp.stat.fpath))
for(groupi in c("cord_blood", "whole_blood", "all",
        "peripheral_blood_mononuclear_cells")){
    message("subgroup: ", groupi); cgv <- unique(unlist(lcgfinal))
    
    # 
    cond <- identical(rownames(bv.hm450k), rownames(bv.epic))
    if(!cond){stop("Couldn't match cgids in gr objects.")}
    # get bval stats
    lbv <- list("hm450k" = bv.hm450k, "epic" = bv.epic, "combined" = bv.combined)
    lstat <- lapply(lbv, function(x){list("mean" = rowMeans(x, na.rm =T),
        "var" = matrixStats::rowVars(x, na.rm = T))})
    # do corr tests and report rhos
    # ct.median <- cor.test(lstat[[1]]$median,lstat[[2]]$median,method="spearman")
    # ct.mean <- cor.test(lstat[[1]]$mean, lstat[[2]]$mean, method = "spearman")
    # ct.var <- cor.test(lstat[[1]]$var, lstat[[2]]$var, method = "spearman")
    # message("c.t. median, rho: ", ct.median$estimate)
    # message("c.t. mean, rho: ", ct.mean$estimate)
    # message("c.t. var, rho: ", ct.var$estimate)
    names(lstat) <- c("hm450k", "epic", "combined")
    lgrp.stat[[groupi]] <- lstat}
# subgroup: cord_blood
# c.t. mean, rho: 0.902686503786354
# c.t. var, rho: 0.830665796016492
# subgroup: whole_blood
# c.t. mean, rho: 0.939576416515208
# c.t. var, rho: 0.869592747439162
# subgroup: all
# c.t. mean, rho: 0.948025800858743
# c.t. var, rho: 0.892401464920395
# subgroup: peripheral_blood_mononuclear_cells
# c.t. mean, rho: 0.974320877005387
# c.t. var, rho: 0.917640295774072

# make figure -- storm plots for mean, var cg diffs

# save means and vars for power analyses
lgrp.stat.fpath <- "lstat-bv-cgoverlap_2platforms_blood-groups.rda"
save(lgrp.stat, file = lgrp.stat.fpath)
