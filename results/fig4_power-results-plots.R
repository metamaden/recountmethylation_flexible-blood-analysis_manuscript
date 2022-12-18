#!/usr/bin/env R

# Author: Sean Maden
#
# Plot pwrEWAS power analysis results across 4 blood sample types.
# 

library(ggplot2)
library(gridExtra)
library(ggpubr)
library(png)
library(patchwork)

#----------
# load data
#----------
# get results list
lpwr.fname <- "lpwr-results_gseadj-combined-2platforms_4groups.rda"
lpwr <- get(load(lpwr.fname))
# define the palette
pal <- c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#868686FF")

#---------------
# get plot data
#---------------
dfp <- do.call(rbind, lapply(names(lpwr), function(groupi){
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

# format vars
for(c in 1:4){
  dfp[,c] <- as.numeric(dfp[,c])}
dfp[,5] <- as.factor(dfp[,5])

#--------------------
# get composite plots
#--------------------
# coordinates for plot insets
ymax.zoomv <- c(0.6, 0.6, 0.7)
ymin.zoomv <- c(0, 0, 0.2)
xmax.zoomv <- c(75, 75, 40)
xmin.zoomv <- c(0, 0, 10)
dv <- c(0.05, 0.1, 0.2)
inset.leftv <- c(0.45, 0.45, 0.3)
inset.topv <- c(0.63, 0.65, 0.7)
xlabv <- c("Sample size", "Sample size", "Sample size")
ylabv <- c("Alpha", "Alpha", "Alpha")
ppdf$Subgroup <- ppdf$group

# coordinates for png images (magnifying glass, right arrpw)
coordv <- list(c(xmin.mg = 60, xmax.mg = 100, ymin.mg = 0, ymax.mg = 0.15,
                 xmin.ar = 95, xmax.ar = 140, ymin.ar = 0, ymax.ar = 0.20),
               c(xmin.mg = 40, xmax.mg = 80, ymin.mg = 0.10, ymax.mg = 0.25,
                 xmin.ar = 75, xmax.ar = 110, ymin.ar = 0.10, ymax.ar = 0.30),
               c(xmin.mg = 30, xmax.mg = 70, ymin.mg = 0.20, ymax.mg = 0.35,
                 xmin.ar = 60, xmax.ar = 95, ymin.ar = 0.20, ymax.ar = 0.40))

# get list of plot objets
lgg <- lapply(seq(3), function(ii){
  deltai <- dv[ii]
  ymax.zoom <- ymax.zoomv[ii]; xmax.zoom <- xmax.zoomv[ii]
  ymin.zoom <- ymin.zoomv[ii]; xmin.zoom <- xmin.zoomv[ii]
  coordvi <- coordv[[ii]]; ppdfi <- ppdf[ppdf$delta == deltai,]
  plot.main <- ggplot(ppdfi, aes(x = sample.size, y = alpha, colours = Subgroup)) + 
    annotate("rect", xmin = xmin.zoom, xmax = xmax.zoom, 
             ymin = ymin.zoom, ymax = ymax.zoom, alpha = 0.34, fill = "grey") +
    stat_smooth(aes(color = group), method="loess", se = F, alpha = 0.5) +
    geom_hline(yintercept = 0.8, color = "black", linetype = "dashed") +
    theme_bw() + xlim(0, 320) + ylim(0,0.95) + ggtitle(paste0("Delta = ", deltai)) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    scale_colour_manual(values = pal) + ylab("Power") +
    annotation_raster(readPNG("magnifying_glass_bgtransparent.png"), 
                      xmin = coordvi[["xmin.mg"]], xmax = coordvi[["xmax.mg"]], 
                      ymin = coordvi[["ymin.mg"]], ymax = coordvi[["ymax.mg"]]) +
    annotation_raster(readPNG("rightarrow_bgtransparent.png"), 
                      xmin = coordvi[["xmin.ar"]], xmax = coordvi[["xmax.ar"]], 
                      ymin = coordvi[["ymin.ar"]], ymax = coordvi[["ymax.ar"]])
  plot.inset <- ggplot(ppdfi, aes(x = sample.size, y = alpha, colours = group)) + 
    stat_smooth(aes(color = group), method="loess", se = F) +
    theme_bw() + ylim(ymin.zoom, ymax.zoom) + xlim(xmin.zoom, xmax.zoom) + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_colour_manual(values = pal)
  plot.composite <- plot.main + 
    inset_element(plot.inset, left = inset.leftv[ii], bottom = 0.02, 
                  right = 0.98, top = inset.topv[ii])
  return(plot.composite)
})

# get legend
ppdf$Subgroup <- ppdf$group
lplot <- ggplot(ppdf, aes(x = sample.size, y = alpha)) + 
  stat_smooth(aes(color = Subgroup), method="loess", se = F) + 
  scale_colour_manual(values = pal) + theme_bw()
lgg[["legend"]] <- as_ggplot(get_legend(lplot))

#--------------
# save new plot
#--------------
# plot vars
ymax.zoomv <- c(0.6, 0.6, 0.7)
ymin.zoomv <- c(0, 0.16, 0.2)
xmax.zoomv <- c(75, 55, 40)
xmin.zoomv <- c(0, 0, 10)
dv <- c(0.05, 0.1, 0.2)
inset.leftv <- c(0.45, 0.35, 0.3)
inset.topv <- c(0.63, 0.7, 0.7)
xlabv <- c("Sample size", "Sample size", "Sample size")
ylabv <- c("Alpha", "Alpha", "Alpha")
ppdf$Subgroup <- ppdf$group

# save new composite pdf
pdf.fname <- "fig3c_pwrewas-results-3deltas-4types.pdf"
# axis labs
xlab.str <- paste0("Total samples (N)", 
                   paste0(rep(" ", 20), collapse = ""), 
                   collapse = "")
ylab.str <- "Power"
# layout matrix
lm <- matrix(c(rep(seq(3), each = 2), 4), nrow = 1)
pdf(pdf.fname, 9, 2.5)
grid.arrange(ggarrange(lgg[[1]]), ggarrange(lgg[[2]]), ggarrange(lgg[[3]]), 
             ggarrange(lgg[[4]]), layout_matrix = lm, 
             bottom = xlab.str, left = ylab.str)
dev.off()