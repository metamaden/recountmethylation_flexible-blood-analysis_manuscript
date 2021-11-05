#!/usr/bin/env R

# Author: Sean Maden
#
# Power analysis using pwrEWAS. Note, this uses a modified pwrEWAS() function, 
# loaded from the script `pwrEWAS_revised.R`.
#
# The supplemental table, main figure 2, and supplemental figure are generated 
# showing the power analysis results.
#
#

library(ggplot2); library(cowplot); library(gridExtra)
library(HDF5Array); library(minfi)
library(methyPre); library(limma); library(sva)

#----------
# load data
#----------
# load data
source("pwrEWAS_revised.R")

# make save dir
# pdir.name <- "power_results"
# if(!dir.exists(pdir.name)){dir.create(pdir.name)}

#-----------------
# helper functions
#-----------------
# 
lpwer_groupi_composite <- function(dfpl, groupi = "cord_blood", 
    xlab.str = "", ylab.str = "", grouplab = NULL){
    # make a composite plot for each delta, with legend
    lp <- list(); dfpi <- dfpl[dfpl$group == groupi,]
    message("Getting plot legend...")
    lplot <- ggplot(data=dfpi[dfpi$delta == "05",], aes(x=num.samples, y=mean, 
            color=platform, group = platform)) + geom_line() + geom_point() + 
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = platform), 
            width=0.5, position=position_dodge(.9)) + theme_bw()
    legend <- get_legend(lplot)
    message("Getting list of composite plots...")
    for(di in c("05", "1", "2")){
        if(is.null(grouplab)){grouplab = groupi}
        ggtitle.str <- ifelse(di == "05", 
            paste0("Group = ",grouplab, "\n"), " \n")
        ylab.str <- ifelse(di == "05", ylab.str, "")
        di.form <- ifelse(di == "05", "0.05",
            ifelse(di == "1", "0.10", "0.20"))
        dfpi.plot <- dfpi[dfpi$delta == di,]
        lp[[di]] <- ggplot(data=dfpi.plot, aes(x=num.samples, y=mean, 
            color=platform, group = platform)) + geom_line() + geom_point() + 
        geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = platform), 
            width=0.5, position=position_dodge(.9)) + theme_bw() +
        ggtitle(paste0(ggtitle.str, "Delta = ", di.form)) + 
        theme(legend.position = "none") + xlab(xlab.str) + ylab(ylab.str) +
        geom_hline(yintercept = 0.8, color = "blue")}
    lp[['legend']] <- legend
    pdf.fpath <- paste0("lineplot-pwr_delta-all_blood-group-",groupi,
        "_platform-all.pdf"); message("Saving plot '", pdf.fpath, "'...")
    pdf(pdf.fpath, 10, 3)
    grid.arrange(lp[[1]], lp[[2]], lp[[3]], legend,
        layout_matrix=matrix(seq(4),nrow=1)); dev.off(); return(lp)
}

lpwer_80 <- function(dfpl, power.perc = 80){
    message("Getting estimated samples at given power of ",
        power.perc, "% power...")
    pwr.dec <- power.perc/100
    dpp <- as.data.frame(matrix(nrow = 0, ncol = 5))
    colnames(dpp) <- c("group", "platform", "delta", 
        "pwr_mean", "pwr_selo", "pwr_sehi")
    for(groupi in unique(dfpl$group)){
        for(platformi in unique(dfpl$platform)){
            for(deltai in unique(dfpl$delta)){
                dfpl.filt <- dfpl$group == "cord_blood"
                dfpl.filt <- dfpl.filt & dfpl$platform == "hm450k"
                dfpl.filt <- dfpl.filt & dfpl$delta == "05"
                dii <- dfpl[dfpl.filt,]
                dii$sdlo <- dii$mean - dii$sd
                dii$sdhi <- dii$mean + dii$sd
                mean.ns = predict.glm(glm(num.samples ~ mean, data = dii),
                    newdata = data.frame(mean = pwr.dec))
                # at est samples, get sdlo, sdhi
                sdlo = predict.glm(glm(sdlo ~ num.samples, data = dii),
                    newdata = data.frame(num.samples = mean.ns))
                sdhi = predict.glm(glm(sdhi ~ num.samples, data = dii),
                    newdata = data.frame(num.samples = mean.ns))
                dppi <- data.frame(group = groupi, platform = platformi,
                    delta = deltai, mean_ns = mean.ns, sdlo = selo, 
                    sdhi = sehi); dpp <- rbind(dpp, dppi)
            }
        }
    }
}

#-------------------------------
# get the summary stats by group
#-------------------------------
lgrp.stat <- list()
# load data

# get dir paths
save.dpath <- ""
gr.fnv <- list.files()
gr.fnv <- gr.fnv[grepl("^gr-adj.*", gr.fnv)]
gr.fnv
# [1] "gr-adj_epic_blood-group-all.rda"
# [2] "gr-adj_epic_blood-group-cord_blood.rda"
# [3] "gr-adj_epic_blood-group-peripheral_blood_mononuclear_cells.rda"
# [4] "gr-adj_epic_blood-group-whole_blood.rda"
# [5] "gr-adj_hm450k_blood-group-all.rda"
# [6] "gr-adj_hm450k_blood-group-cord_blood.rda"
# [7] "gr-adj_hm450k_blood-group-peripheral_blood_mononuclear_cells.rda"
# [8] "gr-adj_hm450k_blood-group-whole_blood.rda"

# load metadata
md.fname <- "si2_blood-md-2platforms.rda"
# md <- get(load(file.path(save.dpath, md.fname)))
md <- get(load(md.fname))

# load the detp tables
# ptable1.path <- file.path(save.dpath, "detp-sstats_hm450k-blood-groups.rda")
# ptable2.path <- file.path(save.dpath, "detp-sstats_epic-blood-groups.rda")
ptable1.path <- "detp-sstats_hm450k-blood-groups.rda"
ptable2.path <- "detp-sstats_epic-blood-groups.rda"
ptable1 <- get(load(ptable1.path)); ptable2 <- get(load(ptable2.path))

#-----------------------
# get lgrf for each group
#------------------------

# manage groups
# append md to colData, filter probes
library(matrixStats); library(methyPre); data(chen_crxcg); lgrf <- list()
liter <- list(); t1 <- Sys.time()
save.fname <- "lstat-pwr_blood-groups_2platforms.rda"
save.fpath <- file.path(save.fname)

# process probes in batches, parallelized by gr object type
parfun_dfstat <- function(typei, num.probes = 50000){
  # grff <- grf.combined
  if(typei == "hm450k"){grff <- grf.hm450k}
  if(typei == "epic"){grff <- grf.epic}
  if(typei == "combined"){grff <- grf.combined}
  bval <- getBeta(grff); 
  indexv <- seq(0, nrow(bval), num.probes); t1_stat <- Sys.time()
  dfstat <- do.call(rbind, lapply(indexv, function(starti){
    endi <- starti + num.probes - 1
    endi <- ifelse(endi > nrow(bval), nrow(bval), endi)
    mbval <- as.matrix(bval[starti:endi,])
    rmi <- rowMeans(mbval); rvi <- rowVars(mbval)
    message(typei, ":", starti, 
            ", time stat: ", Sys.time() - t1_stat)
    dfi <- data.frame(mean = rmi, var = rvi, 
                      stringsAsFactors = F)
    return(dfi)}))
  return(dfstat)
}

# process group probes with batching, parallelization
for(groupi in c("all", "whole_blood", "cord_blood", 
                "peripheral_blood_mononuclear_cells")){
  message("Beginning group ", groupi, ", time: ", Sys.time() - t1)
  fnvi <- gr.fnv[grepl(groupi, gr.fnv)]
  fnvi.hm450k <- fnvi[grepl("hm450k", fnvi)]
  fnvi.epic <- fnvi[grepl("epic", fnvi)]
  gr.hm450k <- HDF5Array::loadHDF5SummarizedExperiment(fnvi.hm450k)
  gr.epic <- HDF5Array::loadHDF5SummarizedExperiment(fnvi.epic)
  message("Appending colData...")
  # hm450k
  mdf <- md[md$platform == "hm450k",]
  mdf <- mdf[rownames(mdf) %in% colnames(gr.hm450k),]
  grf.hm450k <- gr.hm450k[,colnames(gr.hm450k) %in% rownames(mdf)]
  mdf <- mdf[order(match(rownames(mdf), colnames(grf.hm450k))),]
  cond <- identical(rownames(mdf), colnames(grf.hm450k))
  if(cond){colData(grf.hm450k) <- DataFrame(mdf)}else{
    stop("Failed to match md for gr.hm450k.")}
  # epic
  mdf <- md[md$platform == "epic",]
  mdf <- mdf[rownames(mdf) %in% colnames(gr.epic),]
  grf.epic <- gr.epic[,colnames(gr.epic) %in% rownames(mdf)]
  mdf <- mdf[order(match(rownames(mdf), colnames(grf.epic))),]
  cond <- identical(rownames(mdf), colnames(grf.epic))
  if(cond){colData(grf.epic) <- DataFrame(mdf)}else{
    stop("Failed to match md for gr.epic.")}
  message("Filtering on probes...")
  cname <- paste0("perc_above_01.", groupi)
  cg.keep <- unique(c(ptable1[ptable1[,cname] == 0,1], 
                      ptable2[ptable2[,cname] == 0,1]))
  anno <- getAnnotation(gr.epic)
  cg.sexchr <- rownames(anno[anno$chr %in% c("chrY", "chrX"),])
  cg.keep <- cg.keep[!cg.keep %in% cg.sexchr]
  cg.keep <- cg.keep[!cg.keep %in% chen.crxcg]
  length(cg.keep)
  # get the filtered grf objects
  grf.hm450k <- grf.hm450k[rownames(grf.hm450k) %in% cg.keep,]
  #lgrf[[groupi]][["hm450k"]] <- parfun_dfstat("hm450k")
  grf.epic <- grf.epic[rownames(grf.epic) %in% cg.keep,]
  #lgrf[[groupi]][["epic"]] <- parfun_dfstat("epic")
  grf.combined <- combineArrays(grf.hm450k, grf.epic, 
                                outType = "IlluminaHumanMethylation450k")
  rm(grf.hm450k); rm(grf.epic)
  lgrf[[groupi]][["combined"]] <- parfun_dfstat("combined")
  rm(grf.combined)
  #message("Final dim grf.hm450k: ", dim(grf.hm450k))
  #message("Final dim grf.epic: ", dim(grf.epic))
  #message("Final dim grf.combined: ", dim(grf.combined))
  #message("Beginning probe stats calculations...")
  #typev <- c("hm450k", "epic", "combined")
  # names(lgrf[[groupi]]) <- typev
  #message("dfstat hm450k dim: ", dim(lgrf[["hm450k"]]))
  #message("dfstat epic dim: ", dim(lgrf[["epic"]]))
  message("dfstat combined: ", dim(lgrf[["combined"]]))
  message("Saving new data..."); save(lgrf, file = save.fpath)
  message("Finished with group ", groupi, ", time: ", Sys.time() - t1)
}

# save results
save(lgrf, file = save.fpath)

#---------------------------
# power analyses at 2 scales
#---------------------------
# get bval stats, means and vars
lgrf.fname <- "lstat-pwr_blood-groups_2platforms.rda"
lgrf.fpath <- file.path(lgrf.fname)
lgrf <- get(load(lgrf.fpath))

# helper function
parfun_pwr <- function(level, dfi.stat, num.core = 1, num.sim = 100){
  if(level == "zoom"){
    message("Computing power tests at zoomed scale...")
    lpwr <- pwrEWAS_itable(tissueType = dfi.stat, minTotSampleSize = 10, 
                                maxTotSampleSize = 80, SampleSizeSteps = 5, NcntPer = 0.5, 
                                targetDelta = c(0.05, 0.1, 0.2), J = 10000, targetDmCpGs = 500,  
                                detectionLimit = 0.01, DMmethod = "limma", FDRcritVal = 0.05, 
                                core = num.core, sims = num.sim, maxCnt.tau = 100)
  } else{
    message("Computing power tests at large scale...")
    lpwr <- pwrEWAS_itable(tissueType = dfi.stat, minTotSampleSize = 20, 
                                 maxTotSampleSize = 400, SampleSizeSteps = 50, 
                                 NcntPer = 0.5, targetDelta = c(0.05, 0.1, 0.2), 
                                 J = 10000, targetDmCpGs = 500, detectionLimit = 0.01, 
                                 DMmethod = "limma", FDRcritVal = 0.05, core = num.core, 
                                 sims = num.sim, maxCnt.tau = 100)
  }; return(lpwr)
}

# get group/type labels
groupv <- c("cord_blood", "whole_blood", "all",
            "peripheral_blood_mononuclear_cells")
typev <- c("combined")
catv <- paste0(rep(groupv, each = 1), ";", 
               rep(typev, times = length(groupv)))
# power analyses 1: large scale
# do standard power analyses
t1 <- Sys.time()
lapply(catv, function(cati){
  message("Working on category ", cati, "...")
  groupi <- gsub(";.*", "", cati); typei <- gsub(".*;", "", cati)
  results.fname <- paste0("lpwr-results_type-", typei, "_group-", groupi, ".rda")
  results.fpath <- file.path(results.fname)
  message("Getting the stats df..."); dfi.stat <- lgrf[[groupi]][[typei]]
  colnames(dfi.stat) <- c("mu", "var"); levelv <- c("zoom", "large")
  lpwr <- lapply(levelv, function(leveli){parfun_pwr(leveli, dfi.stat)})
  names(lpwr) <- c("zoom", "large"); 
  message("Saving results..."); save(lpwr, file = results.fpath)
  message("Finished with ",groupi,", ", typei, ", time: ", Sys.time() - t1)
})

# load results and save in single list
fnv <- paste0("lpwr-results_type-", gsub(".*;", "", catv), "_group-", 
              gsub(";.*", "", catv), ".rda")
lpwr <- lapply(fnv, function(fn){get(load(fn))})
names(lpwr) <- gsub(";.*", "", catv)
lpwr.fname <- "lpwr-results_gseadj-combined-2platforms_4groups.rda"
save(lpwr, file = lpwr.fname)

#------------------------------
# get the composite power array
#------------------------------
# get results list
lpwr.fname <- "lpwr-results_gseadj-combined-2platforms_4groups.rda"
lpwr <- get(load(lpwr.fname))

ppdf <- do.call(rbind, lapply(names(lpwr), function(groupi){
  lpi <- lpwr[[groupi]]
  ppdf <- do.call(rbind, lapply(lpi, function(da){
    pa <- da$powerArray
    paf <- as.data.frame(as.matrix(pa), stringsAsFactors = F)
    colnames(paf) <- "alpha"; repv <- dimnames(pa)[[1]]
    ssizev <- dimnames(pa)[[2]]; deltav <- dimnames(pa)[[3]]
    paf$rep <- rep(repv, times = length(ssizev)*length(deltav))
    paf$sample.size <- rep(rep(ssizev, each = length(repv)), times = length(deltav))
    paf$delta <- rep(deltav, each = length(repv)*length(ssizev))
    return(paf)})); ppdf$group <- groupi
  return(ppdf)
}))
for(c in 1:4){ppdf[,c] <- as.numeric(ppdf[,c])}
ppdf[grepl("^peripheral.*", ppdf$group),]$group <- "PBMC"
ppdf[,5] <- as.factor(ppdf[,5])

#--------------------
# get composite plots
#--------------------
library(ggpubr); library(gridExtra); library(patchwork)
# library(cowplot)

ymax.zoomv <- c(0.6, 0.6, 0.7)
ymin.zoomv <- c(0, 0, 0.2)
xmax.zoomv <- c(75, 75, 40)
xmin.zoomv <- c(0, 0, 10)
dv <- c(0.05, 0.1, 0.2)
inset.leftv <- c(0.45, 0.35, 0.3)
inset.topv <- c(0.63, 0.7, 0.7)
xlabv <- c("Sample size", "Sample size", "Sample size")
ylabv <- c("Alpha", "Alpha", "Alpha")
ppdf$Subgroup <- ppdf$group

# get the magnifying glass image
# library(png); library(grid)

library(cowplot)
library(magick)
library(ggsci)

# define the palette
pal <- c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#868686FF")
# coordinates for png images (magnifying glass, right arrpw)
coordv <- list(c(xmin.mg = 60, xmax.mg = 100, ymin.mg = 0, ymax.mg = 0.15,
                 xmin.ar = 95, xmax.ar = 140, ymin.ar = 0, ymax.ar = 0.20),
               c(xmin.mg = 50, xmax.mg = 92, ymin.mg = 0, ymax.mg = 0.15,
                 xmin.ar = 85, xmax.ar = 115, ymin.ar = 0, ymax.ar = 0.20),
               c(xmin.mg = 30, xmax.mg = 70, ymin.mg = 0.20, ymax.mg = 0.35,
                 xmin.ar = 60, xmax.ar = 95, ymin.ar = 0.20, ymax.ar = 0.40))
# get list of plot objets
lgg <- lapply(seq(3), function(ii){
  deltai <- dv[ii]
  ymax.zoom <- ymax.zoomv[ii]
  xmax.zoom <- xmax.zoomv[ii]
  ymin.zoom <- ymin.zoomv[ii]
  xmin.zoom <- xmin.zoomv[ii]
  coordvi <- coordv[[ii]]
  ppdfi <- ppdf[ppdf$delta == deltai,]
  plot.main <- ggplot(ppdfi, aes(x = sample.size, y = alpha, colours = Subgroup)) + 
    annotate("rect", xmin = xmin.zoom, xmax = xmax.zoom, 
             ymin = ymin.zoom, ymax = ymax.zoom, alpha = 0.34, fill = "grey") +
    stat_smooth(aes(color = group), method="loess", se = F, alpha = 0.5) +
    geom_hline(yintercept = 0.8, color = "black", linetype = "dashed") +
    theme_bw() + xlim(0, 320) + ylim(0,0.95) + ggtitle(paste0("Delta = ", deltai)) +
    theme(legend.position = "none") + xlab(xlabv[ii]) + ylab(ylabv[ii]) +
    scale_colour_manual(values = pal) +
    annotation_raster(readPNG("magnifying_glass_bgtransparent.png"), 
                      xmin = coordvi[["xmin.mg"]], xmax = coordvi[["xmax.mg"]], 
                      ymin = coordvi[["ymin.mg"]], ymax = coordvi[["ymax.mg"]]) +
    annotation_raster(readPNG("rightarrow_bgtransparent.png"), 
                      xmin = coordvi[["xmin.ar"]], xmax = coordvi[["xmax.ar"]], 
                      ymin = coordvi[["ymin.ar"]], ymax = coordvi[["ymax.ar"]])
    
  plot.inset <- ggplot(ppdfi, aes(x = sample.size, y = alpha, colours = group)) + 
    stat_smooth(aes(color = group), method="loess", se = F) +
    theme_bw() + ylim(ymin.zoom, ymax.zoom) + xlim(xmin.zoom, xmax.zoom) + 
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_colour_manual(values = pal)
  
  plot.composite <- plot.main #+ 
    #inset_element(plot.inset, left = inset.leftv[ii], bottom = 0.02, 
    #              right = 0.98, top = inset.topv[ii])
  return(plot.composite)
})
# get legend
ppdf$Subgroup <- ppdf$group
lplot <- ggplot(ppdf, aes(x = sample.size, y = alpha)) + 
  stat_smooth(aes(color = Subgroup), method="loess", se = F) + 
  scale_colour_manual(values = pal) + 
  theme_bw()
lgg[["legend"]] <- as_ggplot(get_legend(lplot))
# make composite plot
pdf.fname <- "ggsmooth-pwr_composite-3deltas_4groups-2platforms.pdf"
pdf(pdf.fname, 6, 5)
ggarrange(lgg[[1]], lgg[[2]], lgg[[3]], lgg[[4]], ncol = 2, nrow = 2)
dev.off()

# get individual plots
fnv <- paste0("ggsmooth-pwr_delta-",c("05", "1", "2"),
              "_4groups-2platforms.pdf")
lapply(seq(length(fnv)), function(ii){
  pdf.fname <- fnv[ii]; pdf(pdf.fname, 4, 3); print(lgg[[ii]]); dev.off()})
