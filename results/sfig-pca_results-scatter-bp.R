#!/usr/bin/env R

# Author: Sean Maden
#
# Make PCA scatter plots. Makes scatter plots and 95% confidence ellipses
# for top variables explaining PCA variances.
#


library(ggplot2)
library(gridExtra)
library(gtools)
library(ggpubr)
library(grid)

#----------
# load data
#----------
# get plot lists
lpca1.fpath <- "lpca-mpc10-sd-av_fh10k-2plat-autodnam.rda"
lpca <- get(load(lpca1.fpath)); dfp <- as.data.frame(t(lpca$mpc10))
md.anova <- lpca$md # sample metadata

#--------------
# get plot data
#--------------
# get eigenvalues percent vector
eigv <- lpca$pcsd^2
eigv.perc <- round(100*eigv/sum(eigv),2)

#-----------------------------
# format categorical variables
#-----------------------------
# format subgroup var
dfp$subgroup <- md.anova$group
dfp$Subgroup <- ifelse(dfp$subgroup == "other/NOS", "Other/NOS",
                       ifelse(dfp$subgroup == "whole_blood", "Whole blood",
                              ifelse(dfp$subgroup == "cord_blood", "Cord blood", "PBMC")))

# format platform var
dfp$Platform <- md.anova$platform
dfp$Platform <- ifelse(dfp$Platform == "epic", "EPIC", "HM450K")

#----------------------------
# format continuous variables
#----------------------------
# quantile cutoffs for cell fractions
num.quant <- 5
if(min(md.anova$predcell.CD4T) < 0){
  md.anova[md.anova$predcell.CD4T < 0,]$predcell.CD4T <- 0}
if(min(md.anova$predcell.CD8T < 0)){
  md.anova[md.anova$predcell.CD8T < 0,]$predcell.CD8T <- 0}
if(min(md.anova$predcell.Bcell < 0)){
  md.anova[md.anova$predcell.Bcell < 0,]$predcell.Bcell <- 0}
md.anova$predcell.CD4T <- round(md.anova$predcell.CD4T, 2)
md.anova$predcell.CD8T <- round(md.anova$predcell.CD8T, 2)
md.anova$predcell.Bcell <- round(md.anova$predcell.Bcell, 2)
dfp$`CD4+ T-cell quantile` <- gtools::quantcut(md.anova$predcell.CD4T, q = num.quant)
dfp$`CD8+ T-cell quantile` <- gtools::quantcut(md.anova$predcell.CD8T, q = num.quant)
dfp$`B-cell quantile` <- gtools::quantcut(md.anova$predcell.Bcell, q = num.quant)
dfp$gse <- md.anova$gse

# quantile cutoffs for glint.epi.pc1
num.quant <- 5
md.anova$glint.epi.pc1 <- round(md.anova$glint.epi.pc1, 2)
dfp$glint.epi.pc1.quantile <- gtools::quantcut(md.anova$glint.epi.pc1, 
                                               q = num.quant)
md.anova$glint.epi.pc2 <- round(md.anova$glint.epi.pc2, 2)
dfp$glint.epi.pc1.quantile <- gtools::quantcut(md.anova$glint.epi.pc2, 
                                               q = num.quant)

#-----------------
# helper functions
#-----------------
get_ptplots <- function(dfp, varname, eigv.perc = eigv.perc, 
                        xmin = -25, xmax = 105, ymin = -45, ymax = 50, 
                        plot.width = 5, plot.height = 4.5, fname = "newplot.pdf", 
                        pal = NULL, axis.text.size = 16, axis.title.size = 20,
                        legend.text.size = 16, legend.title.size = 20){
  lpt <- list()
  message("Getting axis labels...")
  xlab.str <- paste0("PC1 (",eigv.perc[1],"%)", paste0(rep(" ", 25), collapse = ""), collapse = "")
  ylab.str <- paste0(paste0(rep(" ", 10),collapse =""), "PC2 (",eigv.perc[2],"%)", collapse ="")
  message("Getting plot objects...")
  ggp.str <- paste0("ggplot(dfp, aes(x = PC1, y = PC2, color = ",varname,"))")
  ggplot <- eval(parse(text = ggp.str))
  if(is.null(pal)){
    lpt[["pt.scatter"]] <- ggplot + geom_point(alpha = 0.5) + 
      ylim(ymin, ymax) + xlim(xmin, xmax) + theme_bw() +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = axis.text.size),
            axis.text.y = element_text(size = axis.text.size),
            legend.position = "none")
    lpt[["pt.ellipse"]] <- ggplot + stat_ellipse(level = 0.95) + 
      ylim(ymin, ymax) + xlim(xmin, xmax) + theme_bw() + 
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = axis.text.size),
            axis.text.y = element_text(size = axis.text.size),
            legend.position = "none")
    plot.legend <- ggplot + geom_point(alpha = 0.5) + 
      ylim(ymin, ymax) + xlim(xmin, xmax) + theme_bw() +
      theme(legend.text=element_text(size=legend.text.size),
            legend.title=element_text(size=legend.title.size))
    lpt[["pt.legend"]] <- ggpubr::get_legend(plot.legend)
  }else{
    lpt[["pt.scatter"]] <- ggplot + geom_point(alpha = 0.5) + 
      ylim(ymin, ymax) + xlim(xmin, xmax) +
      theme_bw() + scale_color_manual(values = pal) +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = axis.text.size),
            axis.text.y = element_text(size = axis.text.size),
            legend.position = "none")
    lpt[["pt.ellipse"]] <- ggplot + stat_ellipse(level = 0.95) + 
      ylim(ymin, ymax) + xlim(xmin, xmax) +
      theme_bw() + scale_color_manual(values = pal) +
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = axis.text.size),
            axis.text.y = element_text(size = axis.text.size),
            legend.position = "none")
    plot.legend <- ggplot(dfp, aes(x = PC1, y = PC2, color = Subgroup)) +
      geom_point(alpha = 0.5) + ylim(ymin, ymax) + xlim(xmin, xmax) + theme_bw() + 
      scale_color_manual(values = pal) +
      theme(legend.text=element_text(size=legend.text.size),
            legend.title=element_text(size=legend.title.size))
    lpt[["pt.legend"]] <- ggpubr::get_legend(plot.legend)
  }
  message("Saving pdf...")
  pdf(fname, plot.width, plot.height)
  grid.arrange(lpt[[1]], lpt[[2]], as_ggplot(lpt[[3]]),
               layout_matrix = matrix(c(1,2,3,3), nrow = 2),
               bottom = textGrob(xlab.str, gp=gpar(fontsize=axis.title.size, font=8)), 
               left = textGrob(ylab.str, rot=90, gp=gpar(fontsize=axis.title.size, font=8)))
  dev.off(); return(lpt)
}

#----------------------
# subgroup scatterplots
#----------------------
# get scatterplot list
lpt.pc1v2 <- list()
pal <- c("Other/NOS" = "#3B3B3BFF", "Whole blood" = "#868686FF", 
         "Cord blood" = "#EFC000FF", "PBMC" = "#CD534CFF")

lpt.group <- get_ptplots(dfp, "Subgroup", eigv.perc = eigv.perc, 
                         pal = pal, 
                         fname = "fig2a_pca-pt-subgroup.pdf")

#-------------------------------
# epistructure pc1 scatter plots
#-------------------------------
dfp$`EPISTRUCTURE\nPC1 eigenvalue\nquantile` <- dfp$glint.epi.pc1.quantile
lpt.group <- get_ptplots(dfp, "`EPISTRUCTURE\nPC1 eigenvalue\nquantile`", 
                         eigv.perc = eigv.perc, 
                         fname = "fig2c_pca-pt-epipc1.pdf")

#---------------------------
# platform pc1 scatter plots
#---------------------------
lpt.platform <- get_ptplots(dfp, "Platform",eigv.perc = eigv.perc, 
                            fname = "fig2b_pca-pt-platform.pdf")

#-----------------------------------
# cd4t quantile pc1vs2 scatter plots
#-----------------------------------
dfp$`CD4+ T-cell\nquantile` <- dfp$`CD4+ T-cell quantile`
lpt.cd4t <- get_ptplots(dfp, "`CD4+ T-cell\nquantile`", eigv.perc = eigv.perc, 
                        fname = "fig2e_pca-pt-cd4t.pdf")

#-----------------------------------
# cd8t quantile pc1vs2 scatter plots
#-----------------------------------
dfp$`CD8+ T-cell\nquantile` <- dfp$`CD8+ T-cell quantile`
lpt.cd4t <- get_ptplots(dfp, "`CD8+ T-cell\nquantile`", eigv.perc = eigv.perc, 
                        fname = "fig2d_pca-pt-cd8t.pdf")

#------------------------------------
# bcell quantile pc1vs2 scatter plots
#------------------------------------
dfp$`B-cell\nquantile` <- dfp$`B-cell quantile`
lpt.bcell <- get_ptplots(dfp, "`B-cell\nquantile`", eigv.perc = eigv.perc, 
                         fname = "fig2f_pca-pt-bcell.pdf")

