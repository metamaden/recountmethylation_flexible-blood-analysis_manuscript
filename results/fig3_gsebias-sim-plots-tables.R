#!/usr/bin/env R

# Author: Sean Maden
#
# Make heatmaps of GSE bias simulations results.

library(ggplot2); library(data.table)
library(scales); library(gridExtra)
library(ggridges)

#----------
# load data
#----------
# load data tables
# sum of squared variances table
msq.fname <- "msq-gse-bias_all-blood-2-platforms.rda"
msq <- get(load(msq.fname))
# fev differences table
mdif.fname <- "mdiff-gse-bias_all-blood-2-platforms.rda"
mdif <- get(load(mdif.fname))

#-----------------
# helper functions
#-----------------
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

#---------------------------------------
# fraction explained variances plot data
#---------------------------------------
dfp.fev <- apply(msq[,c(1:39)], 2, function(ci){
  median(as.numeric(ci), na.rm=T)})
# format heatmap data
lvlv <- c("gse", "predsex", "predcell.Mono", "predcell.NK",
          "predcell.CD4T", "predage", "predcell.Bcell", "predcell.CD8T", 
          "predcell.Gran", "platform", "glint.epi.pc2", "glint.epi.pc1", 
          "Residuals")
dfp.fev <- data.frame(var = names(dfp.fev), value = as.numeric(dfp.fev))
dfp.fev$model <- gsub(".*_", "", dfp.fev$var)
dfp.fev$var <- gsub("_.*", "", dfp.fev$var)
dfp.fev$`Median\nFEV` <- as.numeric(dfp.fev$value)
dfp.fev$var <- factor(dfp.fev$var, levels = lvlv)
dfp.fev$model <- factor(dfp.fev$model, 
                        levels = c("unadj", "adj1", "adj2"))
dfp.fev$value.label <- round(100*dfp.fev$value, digits = 2)

#--------------------------------------
# get plot dfp data for fev plots, fast
#--------------------------------------
# get fev binned on type
typev <- c("unadj", "adj1", "adj2")
lvarv <- list(technical = c("platform"),
              demographic = c("predage", "predsex", "glint.epi.pc2", 
                              "glint.epi.pc1"),
              biological = c("predcell.CD8T", "predcell.CD4T", "predcell.NK", 
                             "predcell.Bcell", "predcell.Mono", "predcell.Gran"))

msqf <- msq[,!grepl("^Residuals.*", colnames(msq))]
ltot.fev <- lapply(typev, function(ti){
  apply(msqf[,grepl(ti, colnames(msqf))], 1, sum, na.rm = T)
})
names(ltot.fev) <- typev

# get plot data
dfp <- do.call(rbind, lapply(names(lvarv), function(vari){
  varvii <- lvarv[[vari]]
  do.call(rbind, lapply(names(ltot.fev), function(ti){
    fev.fract.denom <- ltot.fev[[ti]]
    msqff <- msqf[,grepl(ti, colnames(msqf))]
    # get vector of category ssq var
    which.cnamev <- grepl(paste0(varvii, collapse = "|"), colnames(msqff))
    msqff <- msqff[, which.cnamev, drop = F]
    ssqv <- apply(msqff, 1, sum, na.rm = T)
    fev.cat.fractv <- ssqv/fev.fract.denom # get fraction fev by cat
    dfi <- data.frame(fev.fract = fev.cat.fractv)
    dfi$vartype <- vari
    dfi$modeltype <- ti
    return(dfi)
  }))
}))

#--------------------
# make fev dist plots
#--------------------
# load data
#dfp.fname <- "dfp_fev-bycat_gse-bias_blood-4stypes.csv"
#dfp <- data.table::fread(dfp.fname, sep = ",", header = T, data.table = F)

# format vars
dfp$`Model type` <- ifelse(dfp$modeltype == "unadj", "unadjusted",
                           ifelse(dfp$modeltype == "adj1", "adjustment 1", "adjustment 2"))
dfp$FEV <- dfp$fev.fract

# get violin plot objects
ggviolin <- ggplot(dfp, aes(y = FEV, x = `Model type`, fill = `Model type`)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())
ggviolin <- ggviolin + facet_wrap(~vartype, ncol = 4)

# save new plot
plot.fname <- "ggviolin_fev-byvarcat_gsebias.pdf"
pdf(paste0(plot.fname, ".pdf"), 5.6, 1.6)
print(ggviolin); dev.off()

# get summary stats
grpv <- c("unadj", "adj1", "adj2")
filtv <- c("all", "biological", "demographic", "technical")
# all vars
dfpi <- dfp
for(filti in filtv){
  dfpi <- dfp[dfp$vartype==filti,]
  for(grpi in grpv){
    dfpii <- dfpi[dfpi$modeltype == grpi,]
    message(filti, ", ", grpi, ": ", median(dfpii$fev, na.rm = T))
  }
}

# biological, unadj: 0.153656211806912
# biological, adj1: 0.365648631715492
# biological, adj2: 0.373661371181732

# demographic, unadj: 0.254445061411055
# demographic, adj1: 0.611092139972859
# demographic, adj2: 0.615988921116518

# technical, unadj: 0.0369246897171495
# technical, adj1: 0.0857439298076734
# technical, adj2: 0.0885693329339579





#----------------------------
# get data for fev dist plots
#----------------------------
# get fev binned on type
varv.technical <- c("platform")
varv.dem <- c("predage", "predsex", "glint.epi.pc2", "glint.epi.pc1")
varv.bio <- c("predcell.CD8T", "predcell.CD4T", 
              "predcell.NK", "predcell.Bcell",
              "predcell.Mono", "predcell.Gran")
lvarv <- list(technical = varv.technical,
              demographic = varv.dem,
              biological = varv.bio)
typev <- c("unadj", "adj1", "adj2")

# write to new results table
dfp.fname <- "dfp_fev-bycat_gse-bias_blood-4stypes.csv"
mcname <- matrix(c("vartype", "fev", "modeltype"), nrow = 1)
data.table::fwrite(mcname, file = dfp.fname, sep = ",", 
                   row.names = F, col.names = F, append = F)
# iterate on sims
dfp <- do.call(rbind, lapply(seq(nrow(msq)), function(ri){
  #message(ri); 
  ridat <- msq[ri,]
  dfi <- do.call(rbind, lapply(typev, function(ti){
    rii <- ridat[grepl(ti, names(ridat))]
    rii <- rii[!grepl('Residuals', names(rii))]
    total.var <- sum(as.numeric(rii), na.rm = T)
    do.call(rbind, lapply(names(lvarv), function(vi){
      rii.vi <- rii[paste0(lvarv[[vi]], "_", ti)]
      if(!vi == "technical"){rii.vi <- sum(rii.vi, na.rm = T)}
      rfract <- as.numeric(rii.vi)/total.var
      data.frame(vartype = vi, fev = rfract, modeltype = ti)
    }))
  }))
  data.table::fwrite(dfi, file = dfp.fname, sep = ",", 
                     row.names = F, col.names = F, append = T)
}))

# save plot data
dfp.fname <- "dfp_fev-bycat_gse-bias_blood-4stypes.rda"
save(dfp, file = dfp.fname)

#--------------------------------------
# table s2 -- median fevs by model, var
#--------------------------------------
# get var categories
grpv <- c("unadj", "adj1", "adj2")
filtv <- c("biological", "demographic", "technical")
tfev <- do.call(cbind, lapply(grpv, function(grpi){
  dfpi <- dfp[dfp$modeltype==grpi,]
  unlist(lapply(filtv, function(filti){
    dfpii <- dfpi[dfpi$vartype==filti,]
    median(dfpii$fev, na.rm = T)
  }))
}))
colnames(tfev) <- grpv
rownames(tfev) <- filtv
# get variable-wise fev
dim(msq)
msq.filt <- msq[,!grepl("^Residuals.*", colnames(msq))]
dim(msq.filt)
# total vars by sim
ltot.fev <- lapply(grpv, function(grpi){
  apply(msq.filt[,grepl(grpi, colnames(msq.filt))], 1, 
        function(ri){sum(ri, na.rm = T)})
})
names(ltot.fev) <- grpv
# parse fevs by var
filtv <- unique(gsub("_.*", "", colnames(msq.filt)[1:36]))
tfev.bind <- do.call(cbind, lapply(grpv, function(grpi){
  msqi <- msq.filt[,grepl(grpi, colnames(msq.filt))]
  tot.fevi <- ltot.fev[[grpi]]
  unlist(lapply(filtv, function(filti){
    fractvi <- msqi[,grepl(filti, colnames(msqi))]/tot.fevi
    median(fractvi, na.rm = T)
  }))
}))
colnames(tfev.bind) <- grpv
rownames(tfev.bind) <- filtv
# bind all results
st2 <- rbind(tfev, tfev.bind)
t(round(st2, 3))
