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

#----------
# load data
#----------
# load data
source("pwrEWAS_revised.R")
# get bval stats, means and vars
lgrp.stat.fpath <- "lstat-bv-cgoverlap_2platforms_blood-groups.rda"
lgrp.stat <- get(load(lgrp.stat.fpath))
# make save dir
pdir.name <- "power_results"
if(!dir.exists(pdir.name)){dir.create(pdir.name)}

#-----------------
# helper functions
#-----------------
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

#------------------
# do power analyses
#------------------
# do standard power analyses
num.core <- 10; num.sim <- 100
for(groupi in names(lgrp.stat)){
    lgroupi <- lgrp.stat[[groupi]]
    for(platformi in names(lgroupi)){
        results.fname <- paste0("lpwr_platform-", platformi, "_blood-group-",
            groupi, ".rda");results.fpath <- file.path(pdir.name, results.fname)
        if(file.exists(results.fpath)){
            message("Found results file ", results.fpath, ", skipping tests...")
        } else{
            message("Computing power tests for new file ", results.fpath)
            lplatformi <- lgroupi[[platformi]]
            dfi <- data.frame(mu = lplatformi$mean, 
                var = lplatformi$var, stringsAsFactors = F)
            lpwr <- pwrEWAS_itable(tissueType = dfi, minTotSampleSize = 50, 
                maxTotSampleSize = 1100, SampleSizeSteps = 150, NcntPer = 0.5, 
                targetDelta = c(0.05, 0.1, 0.2), J = 10000, targetDmCpGs = 500,  
                detectionLimit = 0.01, DMmethod = "limma", FDRcritVal = 0.05, 
                core = num.core, sims = num.sim, maxCnt.tau = 100)
            save(lpwr, file = results.fpath)
        }; message("Finished with results file ", results.fpath, "...")
    }
}

# test lpwr 
lpwr.fpath <- file.path(pdir.name, "lpwr_platform-combined_blood-group-all.rda")
lpwr <- get(load(lpwr.fpath))

# get the lpwr results as a df
groupv <- c("all", "cord_blood", "peripheral_blood_mononuclear_cells", 
    "whole_blood")
platformv <- c("hm450k", "epic", "combined")

df.pa <- do.call(rbind, lapply(groupv, function(groupi){
    dfi <- do.call(rbind, lapply(platformv, function(platformi){
        lpwr.fname <- paste0("lpwr_platform-",platformi,
            "_blood-group-",groupi,".rda")
        lpwr.fpath <- file.path(pdir.name, lpwr.fname)
        lpwr <- get(load(lpwr.fpath))
        dfi.pa <- as.data.frame(lpwr$powerArray)
        dfi.pa$group <- groupi; dfi.pa$platform <- platformi 
        return(dfi.pa)})); return(dfi)}))

df.mean <- do.call(rbind, lapply(groupv, function(groupi){
    dfi <- do.call(rbind, lapply(platformv, function(platformi){
        lpwr.fname <- paste0("lpwr_platform-",platformi,
            "_blood-group-",groupi,".rda")
        lpwr.fpath <- file.path(pdir.name, lpwr.fname)
        lpwr <- get(load(lpwr.fpath))
        dfi.mean <- as.data.frame(lpwr$meanPower)
        dfi.mean$group <- groupi; dfi.mean$platform <- platformi 
        return(dfi.mean)})); return(dfi)}))

lpwr <- list(pa = df.pa, mean = df.mean)
lpwr.fpath <- "lpwr_pwrmd-100sims-3delta-500dmp_2platforms_blood-groups.rda"
save(lpwr, file = lpwr.fpath)

# transfer to local
# scp -P 21747 metamaden@69.168.53.28:/home/metamaden/bioinfo_appnote/manuscript_results/lpwr_pwrmd-100sims-3delta-500dmp_2platforms_blood-groups.rda ./

#--------------------------
# get power composite fig 2
#--------------------------
lpwr <- get(load(lpwr.fpath))
dfi <- lpwr[["pa"]]; dfi <- dfi[dfi$platform == "combined",]

# save supp table
st.fname <- "st_power-results-std_4groups_2platforms"
save(dfi, file = paste0(st.fname, ".rda"))
write.csv(dfi, file = paste0(st.fname, ".csv"))

# make dft for plotting fig 2
dft <- as.data.frame(matrix(nrow = 0, ncol = 5))
colnames(dft) <- c("num.samples", "power", "rep", "delta", "group")
for(groupi in unique(dfi$group)){
    for(di in c("05", "1", "2")){
        cnv.format <- gsub(".*\\.", "", colnames(dfi))
        col.filt <- cnv.format == di; row.filt <- dfi$group == groupi
        dffi <- dfi[row.filt, col.filt]
        num.samplesv <- gsub("\\..*", "", colnames(dffi))
        num.samplesv <- rep(num.samplesv, each = 100)
        repv <- rep(seq(100), 8)
        dati <- unlist(dffi)
        dfbind <- data.frame(num.samples = num.samplesv,
            power = dati, rep = repv, delta = di, group = groupi,
            stringsAsFactors = F)
        dft <- rbind(dft, dfbind)}}
dft$num.samples <- as.numeric(dft$num.samples)
dft[grepl("^peripheral.*", dft$group),]$group <- "PBMC"

lgg <- list()
for(di in c("05", "1", "2")){
    di.labstr <- ifelse(di == "05", "0.05",
        ifelse(di == "1", "0.10", "0.20"))
    xlab.str <- ifelse(di == "1", "Number of samples", "")
    ylab.str <- ifelse(di == "05", "Power", "")
    lgg[[di]] <- ggplot(dft[dft$delta == di,], 
        aes(x=num.samples, y=power, color=group, group=group)) +
        scale_fill_manual(values = cbp1.map) +
        geom_hline(yintercept=0.8, color = "blue") +
        geom_smooth(method = "loess") + theme_bw() +
        ylim(0.5, 1) + xlim(0, 1200) +
        xlab(xlab.str) + ylab(ylab.str) + 
        ggtitle(paste0("Delta = ", di.labstr)) +
        theme(legend.position = "none",
            axis.text.x = element_text(angle = 90, hjust = 0.99, vjust = 0.2))}

dft$`Blood group` <- dft$group
ploti <- ggplot(dft, 
    aes(x=num.samples, y=power, color=`Blood group`, group=group)) + 
    theme_bw() + geom_smooth(method = "loess")
lgg[["legend"]] <- get_legend(ploti)

lm <- matrix(c(rep(c(1,2,3), each = 2), 4), nrow = 1)
pdf.fname <- "ggsmooth-composite_platform-combined_4groups.pdf"
pdf(pdf.fname, 8.2, 2.5); 
grid.arrange(lgg[[1]], lgg[[2]], lgg[[3]], lgg[[4]],
    layout_matrix = lm)
dev.off()

#-----------------------------
# get power composite supp fig
#-----------------------------
# make dfp
lpwr <- get(load(lpwr.fpath))
dfp <- lpwr[["mean"]]; colnames(dfp)[1:3] <- paste0(c("05", "1", "2"), ".mean")
dfp$num.samples <- rep(rownames(dfp)[1:8], 12)
dfi <- lpwr[["pa"]]; deltav <- c("05", "1", "2")
dfq <- do.call(cbind, lapply(deltav, function(deltai){
    cnv <- gsub(".*\\.", "", colnames(dfi));which.deltai <- grepl(deltai, cnv)
    dfii <- dfi[,which.deltai];colnames(dfii) <- gsub("\\..*","",colnames(dfii))
    dfq.all <- as.data.frame(matrix(nrow = 0, ncol = 5))
    colnames(dfq.all) <- c(paste0(deltai, ".", c("q25", "q75", "sd")), 
        "group", "platform")
    for(groupi in unique(dfi$group)){
        for(platformi in unique(dfi$platform)){
            samples.filt <- dfi$group == groupi & dfi$platform == platformi
            dfii.filt <- dfii[which(samples.filt),]
            q25v <- apply(dfii.filt, 2, function(x){quantile(x)[2]})
            q75v <- apply(dfii.filt, 2, function(x){quantile(x)[4]})
            sev <- apply(dfii.filt, 2, function(x){sd(x)})
            dfq <- data.frame(q25v, q75v, sev, stringsAsFactors = F)            
            colnames(dfq) <- paste0(deltai, ".", c("q25", "q75", "sd"))
            rownames(dfq) <- colnames(dfii.filt)
            dfq$group <- groupi; dfq$platform <- platformi
            dfq.all <- rbind(dfq.all, dfq)}}
    return(dfq.all)}))
dfq <- dfq[,c(1:3, 6:8,11:15)]
cond <- identical(dfp$platform, dfq$platform) & 
    identical(dfp$group, dfq$group)
if(!cond){stop("Couldn't match groups and platforms for dfp and dfq.")}
dfp.final <- cbind(dfp, dfq)
# save
dfp.final.fpath <- "table_pwr-results-std_all-groups-2platforms"
save(dfp.final, file = paste0(dfp.final.fpath, ".rda"))
write.csv(dfp.final, file = paste0(dfp.final.fpath, ".csv"))

# make plots
dfp <- get(load(paste0(dfp.final.fpath, ".rda")))
# dfp <- get(load(paste0(dfp.final.fpath, ".rda")))
# make dfp long
dfp.pre <- get(load(paste0(dfp.final.fpath, ".rda")))
dfpl <- do.call(rbind, lapply(c("05", "1", "2"), function(di){
    di.str.pattern <- paste0("^",di,"\\..*")
    dfpi <- dfp.pre[,grepl(di.str.pattern, colnames(dfp.pre))]
    dfpi$delta <- di; colnames(dfpi) <- gsub(".*\\.", "", colnames(dfpi))
    dfpi$group <- dfp.pre$group; dfpi$platform <- dfp.pre$platform
    dfpi$num.samples <- dfp.pre$num.samples
    return(dfpi)}))
dfpl$num.samples <- as.numeric(dfpl$num.samples)

# single plot, group all, delta 0.05, all platforms
groupi <- "all"; platformi <- "any"; deltai <- "05"
which.dfpl <- dfpl$group=="all" & dfpl$delta=="05";dfp.plot<-dfpl[which.dfpl,]

ploti <- ggplot(data=dfp.plot, aes(x=num.samples, y=mean, 
    color=platform, group = platform)) +
  geom_line() + geom_point() + 
  geom_errorbar(data = dfp.plot, aes(ymin=mean-sd, ymax=mean+sd, color = platform), 
    width=0.5, position=position_dodge(.9)) + 
  theme_bw()

plot.fname <- paste0("lineplot-pwr_delta-",deltai,"_blood-group-",groupi,
    "_platform-",platformi,".pdf")
pdf(plot.fname, 5, 4); ploti; dev.off()

# composite plot, group cord_blood, all platforms, all deltas
lp <- list(); groupi <- "cord_blood"; dfpi <- dfpl[dfpl$group == groupi,]
# get the legend plot
lplot <- ggplot(data=dfpi[dfpi$delta == "05",], aes(x=num.samples, y=mean, 
        color=platform, group = platform)) + geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = platform), 
        width=0.5, position=position_dodge(.9)) + theme_bw()
legend <- get_legend(lplot)
# get the composite plots
for(di in c("05", "1", "2")){
    di.form <- ifelse(di == "05", "0.05",
        ifelse(di == "1", "0.10", "0.20"))
    dfpi.plot <- dfpi[dfpi$delta == di,]
    lp[[di]] <- ggplot(data=dfpi.plot, aes(x=num.samples, y=mean, 
        color=platform, group = platform)) + geom_line() + geom_point() + 
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = platform), 
        width=0.5, position=position_dodge(.9)) + theme_bw() +
    ggtitle(paste0("Delta = ", di.form)) + theme(legend.position = "none") +
    geom_hline(yintercept = 0.8, color = "blue")}
pdf.fpath <- paste0("lineplot-pwr_delta-all_blood-group-",groupi,
    "_platform-all.pdf")
pdf(pdf.fpath, 10, 3)
grid.arrange(lp[[1]], lp[[2]], lp[[3]], legend,
    layout_matrix = matrix(c(1,2,3,4), nrow = 1))
dev.off()

# make composite of the group-wise plots
llpwr <- list()
for(ii in seq(4)){
    groupi <- unique(dfpl$group)[ii]
    grouplab.str <- ifelse(grepl("^peripheral.*", groupi), "PBMC", groupi)
    llpwr[[groupi]] <- lpwer_groupi_composite(dfpl, groupi = groupi,
        grouplab = grouplab.str)}
lm <- matrix(c(1,2,3,14,4,5,6,7,8,9,10,7,11,12,13,14), 
    byrow = T, nrow = 4)
pdf.fpath <- "lplot-pwr-std_by-deltas_4groups_2platforms.pdf"
pdf(pdf.fpath, 10, 10)
grid.arrange(
    llpwr[[1]][[1]], llpwr[[1]][[2]], llpwr[[1]][[3]], 
    llpwr[[2]][[1]], llpwr[[2]][[2]], llpwr[[2]][[3]], llpwr[[2]][[4]], 
    llpwr[[3]][[1]], llpwr[[3]][[2]], llpwr[[3]][[3]], 
    llpwr[[4]][[1]], llpwr[[4]][[2]], llpwr[[4]][[3]], 
    layout_matrix = lm, bottom = "Number of samples",
    left = "Power")
dev.off()

# get the smooths for combined probes at each delta
# define colorblind friendly palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
platformi <- "combined"; di <- "05"
which.dfpl <- dfpl$platform==platformi & dfpl$delta==di;
dfp.plot <- dfpl[which.dfpl,]
dfp.plot[grepl("^periph.*", dfp.plot$group),]$group <- "PBMC"

ploti <- ggplot(data=dfp.plot, 
    aes(x = num.samples, y = mean, color = group, group = group)) + 
    geom_smooth(method = loess) + theme_bw() +
    geom_hline(yintercept = 0.8, color = "blue") +
    scale_fill_manual(values = cbp1)

pdf("smooth-glm_exe.pdf"); ploti; dev.off()


  geom_line() + geom_point() + 
  geom_errorbar(data = dfp.plot, aes(ymin=mean-sd, ymax=mean+sd, color = platform), 
    width=0.5, position=position_dodge(.9)) + 
  theme_bw()


# get prediction
dfpl.filt <- dfpl$group == "cord_blood"
dfpl.filt <- dfpl.filt & dfpl$platform == "hm450k"
dfpl.filt <- dfpl.filt & dfpl$delta == "05"

dii <- dfpl[dfpl.filt,]
lm <- lm(num.samples ~ mean, data = dii)
predict(lm, newdata = data.frame(mean = 0.80))

predict.glm(glm(mean ~ num.samples, data = dii))

