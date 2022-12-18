#!/usr/bin/env R

# Author: Sean Maden
#
# Make heatmaps of GSE bias simulations results.

library(ggplot2); library(data.table)
library(scales); library(gridExtra)
library(ggpubr)
# library(magick)

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

#-----------------------------------------------
# main fig -- compare var cat dist, violin plots
#-----------------------------------------------
# get plot data
# get fev binned on type
typev <- c("unadj", "adj1", "adj2")
lvarv <- list(technical = c("platform", "gse"),
              demographic = c("predage", "predsex", "glint.epi.pc2", 
                              "glint.epi.pc1"),
              biological = c("predcell.CD8T", "predcell.CD4T", "predcell.NK", 
                             "predcell.Bcell", "predcell.Mono", "predcell.Gran"))

msqf <- msq[,!grepl("^Residuals.*", colnames(msq))]
ltot.fev <- lapply(typev, function(ti){apply(msqf[,grepl(ti, colnames(msqf))], 1, sum, na.rm = T)})
names(ltot.fev) <- typev
# get plot data object
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

# get plot objects
dfp$`Model type` <- ifelse(dfp$modeltype=="unadj", "unadjusted",
                           ifelse(dfp$modeltype=="adj1", "adjustment 1", "adjustment 2"))
lvlv <- c("unadjusted", "adjustment 1", "adjustment 2")
dfp$`Model type` <- factor(dfp$`Model type`, levels = lvlv)
dfp$FEV <- dfp$fev.fract
# format plot vars
catv <- c("technical", "demographic", "biological")
tech.str <- paste0(paste0(rep(" ", 13), collapse = ""), "Technical", collapse = "")
biol.str <- paste0(paste0(rep(" ", 5), collapse = ""), "Biological", collapse = "")
demo.str <- paste0(paste0(rep(" ", 2), collapse = ""), "Demographic", collapse = "")
# get list of plot objects
text.size <- 10; title.size <- 12
lgg <- lapply(catv, function(cati){
  dfpi <- dfp[dfp$vartype == cati,]
  ggvp <- ggplot(dfpi, aes(y = FEV, x = `Model type`, fill = `Model type`)) + 
    geom_violin(draw_quantiles = 0.5) + theme_bw() +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          legend.position = "none", plot.title = element_text(size = title.size),
          axis.text.y = element_text(size = text.size))
  if(cati == "biological"){
    ggvp <- ggvp + ggtitle(biol.str) + 
      theme(axis.title.y = element_text(size = title.size))
  }
  if(cati == "demographic"){
    ggvp <- ggvp + ggtitle(demo.str) +
      theme(axis.title.y = element_blank())
  }
  if(cati == "technical"){
    ggvp <- ggvp + ggtitle(tech.str) +
      theme(axis.title.y = element_blank())
  }
  return(ggvp)
})
names(lgg) <- catv
# get zoom panel for technical
# get plot legend
pl <- ggplot(dfp, aes(y = FEV, x = `Model type`, fill = `Model type`)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  theme(legend.title = element_text(size = title.size),
        legend.text = element_text(size = text.size))
lgg[["legend"]] <- get_legend(pl)

# save new plot
plot.fname <- "ggvp_fev-byvarcat_gsebias"
mg.pgn.fname <- "magnifying_glass_bgtransparent.png"
# get plot params
lm <- matrix(c(1,1,2,2,3,3,4), nrow = 1)
# save new pdf
pdf(paste0(plot.fname, ".pdf"), 7.8, 1.8)
grid.arrange(lgg[["biological"]], lgg[["demographic"]], lgg[["technical"]], 
             lgg[["legend"]], layout_matrix = lm)
dev.off()


#--------------------
# sfigs, compare fevs
#--------------------
# get plot data
dfp1 <- data.frame(unadj = ltot.fev$unadj, adj.val = ltot.fev$adj1)
dfp2 <- data.frame(unadj = ltot.fev$unadj, adj.val = ltot.fev$adj2)
dfp1$adj.type <- "adj. 1"; dfp2$adj.type <- "adj. 2"
dfp <- rbind(dfp1, dfp2)
# fract fev
dfp$fract.fev <- dfp$adj.val/dfp$unadj

# plot scatterplot fev
ggpt <- ggplot(dfp, aes(x = unadj, y = adj.val)) +
  geom_point(draw_quantiles = 0.1) + theme_bw()
ggpt <- ggpt + facet_wrap(~adj.type, ncol = 2)
pdf("ggpt_fev-adj-unadj_gsebias.pdf", 5.5, 3.5)
print(ggpt); dev.off()

# 2d density plot
ggpt <- ggplot(dfp, aes(x = unadj, y = adj.val)) + geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") + theme_bw() +
  xlab("Unadjusted FEV") + ylab("Adjusted FEV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggpt <- ggpt + facet_wrap(~adj.type, ncol = 2)
pdf("ggdensity_fev-adj-unadj_gsebias.pdf", 3.5, 1.8)
print(ggpt); dev.off()

# plot fraction fev
# get plot object
ggvp <- ggplot(dfp, aes(x = adj.type, y = fract.fev, group = adj.type)) +
  geom_violin(show_quantiles = 0.5) + theme_bw() + 
  ylab("FEV fraction\n(Adj./Unadj.)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
# save new plot
pdf('ggvp_fev-fract_gsebias.pdf', 2.5, 1.5)
print(ggvp);dev.off()

# summary stats for reporting
median(dfp[dfp$adj.type=="adj. 1",]$fract.fev) # 0.6882031
median(dfp[dfp$adj.type=="adj. 2",]$fract.fev) # 0.6841613
var(dfp[dfp$adj.type=="adj. 1",]$fract.fev) # 0.07646976
var(dfp[dfp$adj.type=="adj. 2",]$fract.fev) # 0.07648388
sd(dfp[dfp$adj.type=="adj. 1",]$fract.fev) # 0.2765317
sd(dfp[dfp$adj.type=="adj. 2",]$fract.fev) # 0.2765572

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


#---------------------------------
# violin plots with technical zoom -- OLD
#---------------------------------
source("facet_zoom2.R")
# library(png)

catv <- c("technical", "demographic", "biological")
tech.str <- paste0(paste0(rep(" ", 13), collapse = ""),
                   "Technical", collapse = "")
biol.str <- paste0(paste0(rep(" ", 5), collapse = ""),
                   "Biological", collapse = "")
demo.str <- paste0(paste0(rep(" ", 2), collapse = ""),
                   "Demographic", collapse = "")
text.size <- 10
title.size <- 12
lgg <- lapply(catv, function(cati){
  dfpi <- dfp[dfp$vartype == cati,]
  ggvp <- ggplot(dfpi, aes(y = FEV, x = `Model type`, fill = `Model type`)) + 
    geom_violin(draw_quantiles = 0.5) + theme_bw() +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          legend.position = "none", plot.title = element_text(size = title.size),
          axis.text.y = element_text(size = text.size))
  if(cati == "biological"){
    ggvp <- ggvp + ggtitle(biol.str) + 
      theme(axis.title.y = element_text(size = title.size))
  }
  if(cati == "demographic"){
    ggvp <- ggvp + ggtitle(demo.str) +
      theme(axis.title.y = element_blank())
  }
  if(cati == "technical"){
    ggvp <- ggvp + ggtitle(tech.str) +
      theme(axis.title.y = element_blank()) +
      facet_zoom2(ylim = c(0, 0.01))
  }
  return(ggvp)
})
names(lgg) <- catv
# get zoom panel for technical
# get plot legend
pl <- ggplot(dfp, aes(y = FEV, x = `Model type`, fill = `Model type`)) + 
  geom_violin(draw_quantiles = 0.5) + theme_bw() +
  theme(legend.title = element_text(size = title.size),
        legend.text = element_text(size = text.size))
lgg[["legend"]] <- get_legend(pl)

# save new plot
plot.fname <- "ggviolin_fev-byvarcat_gsebias"
mg.pgn.fname <- "magnifying_glass_bgtransparent.png"
# get plot params
lm <- matrix(c(1,1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3,3,3,3,3,3,
               4,4,4,4,4,4), nrow = 1)
# save new pdf
pdf(paste0(plot.fname, ".pdf"), 7.8, 1.8)
grid.arrange(lgg[["biological"]], lgg[["demographic"]],
             lgg[["technical"]], lgg[["legend"]], 
             layout_matrix = lm)
dev.off()

# Produce image using graphics device
# fig <- image_graph(width = 800, height = 200, res = 110)
# ggplot2::qplot(mpg, wt, data = mtcars, colour = cyl)
#grid.arrange(lgg[["biological"]], lgg[["demographic"]], lgg[["technical"]], 
#             lgg[["legend"]], layout_matrix = lm)
#dev.off()
#mg.image <- image_scale(image_read(mg.pgn.fname), "x22")
#out <- image_composite(fig, mg.image, 
#                       offset = geometry_point(475, -175))
#print(out)