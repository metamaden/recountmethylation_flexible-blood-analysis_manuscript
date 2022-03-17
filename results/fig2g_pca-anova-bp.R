#!/usr/bin/env R

# Author: Sean Maden
#
# Get stacked barplots showing results of performing ANOVAs on autosomal DNAm PCAs.
# Stacked barplots and axis labels show the percent of variances explained for each
# of the top 10 components in PCA.
#

library(ggplot2)
library(gridExtra)
library(gtools)
library(ggpubr)
library(grid)

#-----------------
# helper functions
#-----------------
format_varlab <- function(avdf, pcsdv, num.round = 2,
                          varv = c("group", "platform", "predage", "predcell.Bcell",
                                   "predcell.CD4T", "predcell.CD8T", "predcell.Gran",
                                   "predcell.Mono", "predcell.NK", "predsex", "gse", 
                                   "glint.epi.pc1", "glint.epi.pc2", "Residuals"),
                          varlabv = c("Blood subgroup", "Array platform", "Predicted age",
                                      "Predicted B-cell fraction", "Predicted CD4+ T-cell fraction",
                                      "Predicted CD8+ T-cell fraction", "Predicted Granylocyte fraction",
                                      "Predicted Moncyte fraction", "Predicted NK cell fraction",
                                      "Predicted sex", "Study/GSE ID", "EPISTRUCTURE PC1", 
                                      "EPISTRUCTURE PC2", "Residuals")){
  # format var labels
  avdf$var.label <- "NA"
  for(ii in seq(length(varv))){
    if(varv[ii] %in% avdf$var){
      avdf[avdf$var == varv[ii],]$var.label <- varlabv[ii]
    }
  }
  avdf$pca.label <- "NA"; pcav <- c() # format pca labels
  # get eigenvalue magnitudes, percentages
  eigv <- pcsdv^2; eigv.perc <- round(100*eigv/sum(eigv),num.round)
  for(pci in seq(max(avdf$pc))){
    labi <- paste0("PC",pci," (",eigv.perc[pci],"%)"); pcav <- c(pcav, labi)
    avdf[avdf$pc == pci,]$pca.label <- labi}
  avdf$pca.label <- factor(avdf$pca.label, levels = pcav[rev(order(pcsdv))])
  return(avdf)
}

# make composite barplots
get_lbp <- function(avdf, pcsd, fname = "newplot.pdf", ncol.legend = 1,
                    plot.width = 20, plot.height = 4.5, 
                    lm = matrix(c(1, 1, 1, 2, 2, 2, 3, 3), nrow = 1),
                    axis.text.size = 20, axis.title.size = 25,
                    legend.text.size = 20, legend.title.size = 25,
                    xlab.str = paste0("Principal component (% total variance)", 
                                      paste0(rep(" ", 30), collapse = ""),collapse = ""),
                    varorder = NA,
                    varv = c("group", "platform", "predage", "predcell.Bcell",
                             "predcell.CD4T", "predcell.CD8T", "predcell.Gran",
                             "predcell.Mono", "predcell.NK", "predsex", "gse", 
                             "glint.epi.pc1", "glint.epi.pc2", "Residuals"),
                    varlabv = c("Blood subgroup", "Array platform", "Predicted age",
                                "Predicted B-cell fraction", "Predicted CD4+ T-cell fraction",
                                "Predicted CD8+ T-cell fraction", "Predicted Granylocyte fraction",
                                "Predicted Moncyte fraction", "Predicted NK cell fraction",
                                "Predicted sex", "Study/GSE ID", "EPISTRUCTURE PC1", 
                                "EPISTRUCTURE PC2", "Residuals")){
  dfp <- format_varlab(avdf = avdf, pcsdv = pcsd, varv = varv, varlabv = varlabv)
  dfp <- dfp[order(dfp$pc),]
  num.eigv <- length(which(dfp$pc=="1"))
  dfp$mssq.scale <- dfp$mssq.percv*rep(eigv[1:10], each = num.eigv) # get scaled mssq
  dfp$pc <- factor(dfp$pc, levels = seq(10))
  dfp$Variable <- dfp$var.label # get pca total var
  message("Ordering variables...")
  orderv <- unique(dfp$Variable)
  if(is.na(varorder)){
    orderv <- orderv[order(dfp[dfp$pc == 1,]$mssq.percv)] 
  } else{
    orderv <- varorder
  }
  dfp$Variable <- factor(dfp$Variable, levels = orderv)
  message("Getting plot objects...")
  lbp <- list(dfp = dfp)
  lbp[["bp.avdf"]] <- ggplot(dfp, aes(x = pca.label, y = mssq.scale, fill = Variable)) +
    geom_bar(stat = "identity", color = "black") + theme_bw() + ylab("Eigenvalue") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = axis.text.size),
          axis.text.y = element_text(size = axis.text.size),
          axis.title.x = element_blank(), axis.title.y = element_text(size = axis.title.size),
          legend.position = "none")
  lbp[["bp.avdf.perc"]] <- ggplot(dfp, aes(x = pca.label, y = mssq.scale, fill = Variable)) +
    geom_bar(stat = "identity", position = "fill", color = "black") + 
    theme_bw() + xlab("") + ylab("Eigenvalue (%)") + scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = axis.text.size),
          axis.text.y = element_text(size = axis.text.size),
          axis.title.x = element_blank(), axis.title.y = element_text(size = axis.title.size),
          legend.position = "none")
  # get legend
  legend.plot <- ggplot(dfp, aes(x = pca.label, y = mssq.scale, fill = Variable)) +
    geom_bar(stat = "identity", color = "black") + guides(fill=guide_legend(ncol=ncol.legend)) +
    theme(legend.text=element_text(size=legend.text.size),
          legend.title=element_text(size=legend.title.size))
  lbp[["bp.avdf.legend"]] <- get_legend(legend.plot)
  # save plot
  message("Saving new plot...")
  pdf(fname, plot.width, plot.height)
  grid.arrange(lbp[["bp.avdf"]], lbp[["bp.avdf.perc"]], 
               as_ggplot(lbp[["bp.avdf.legend"]]), layout_matrix = lm, 
               bottom = textGrob(xlab.str, gp=gpar(fontsize=axis.title.size, font = 8)))
  dev.off(); return(lbp)
}

#----------
# load data
#----------
# get pca anovas
lpca1.fpath <- "lpca-mpc10-sd-av_fh10k-2plat-autodnam.rda"
lpca <- get(load(lpca1.fpath))
md.anova <- lpca$md; dfp <- as.data.frame(t(lpca$mpc10))

# get eigenvalues percent vector
eigv <- lpca$pcsd^2; eigv.perc <- round(100*eigv/sum(eigv),2)

# format subgroup var
dfp$subgroup <- md.anova$group
dfp$Subgroup <- ifelse(dfp$subgroup == "other/NOS", "Other/NOS",
                       ifelse(dfp$subgroup == "whole_blood", "Whole blood",
                              ifelse(dfp$subgroup == "cord_blood", "Cord blood", "PBMC")))
# format platform var
dfp$Platform <- md.anova$platform
dfp$Platform <- ifelse(dfp$Platform == "epic", "EPIC", "HM450K")

#------------------------
# make composite barplots
#------------------------
lbp <- get_lbp(lpca$avdf, lpca$pcsd, 
               fname = "fig2_pca-bp-anova.pdf",
               plot.height = 5.5)

lbp.nogse <- get_lbp(lpca$avdf.nogse, lpca$pcsd, 
                     fname = "sfig_pca-anova-nogse.pdf",
                     plot.width = 10, plot.height = 5.5,
                     lm = matrix(c(1,1,2,2,3,3,3), nrow = 1),
                     xlab.str = paste0("Principal component (percent variance)", 
                                       paste0(rep(" ", 30), collapse = ""),collapse = ""))

lbp.noglint <- get_lbp(lpca$avdf.noglint, lpca$pcsd, 
                       fname = "sfig_pca-anova-noglint.pdf",
                       plot.width = 10, plot.height = 5.5,
                       lm = matrix(c(1,1,2,2,3,3,3), nrow = 1),
                       xlab.str = paste0("Principal component (percent variance)", 
                                         paste0(rep(" ", 30), collapse = ""),collapse = ""))

#-----------------------
# make filtered barplots
#-----------------------
# filter on top vars
varvf = c("glint.epi.pc1", "predcell.CD4T", "gse")
varvf <- c(varvf, "Residuals")
avdf <- lpca$avdf
avdf.inc <- avdf[avdf$var %in% varvf,]
avdf.oth <- avdf[!avdf$var %in% varvf,]
# get cumulative perc
oth.perc <- unlist(lapply(seq(10), function(ii){
  sum(avdf.oth[avdf.oth$pc==ii,]$mssq.percv)
}))
# bind final df
avdf.filt <- rbind(avdf.inc, 
              data.frame(mssq.percv = oth.perc, 
                         var = rep("other", 10),
                         pc = seq(10)))

# make barplots
pcsd <- lpca$pcsd
lbp <- get_lbp(avdf.filt, lpca$pcsd, 
               fname = "fig2_pca-bp-anova-filtvar.pdf",
               plot.height = 5.5,
               varv = c("gse", "predcell.CD4T", "glint.epi.pc1", "other", "Residuals"), 
               varlabv = c("Study ID", "CD4+ T-cells", "EPC1", "Other (10 variables)", "Residuals"),
               varorder = c("EPC1", "CD4+ T-cells", "Study ID", "Other (10 variables)", "Residuals"))