#!/usr/bin/env R

# Author: Sean Maden
#
# Set analysis of filtered probes across blood groups and platforms

library(UpSetR)
library(ggplot2)

# load data
lcg.fpath <- paste0("list-cgv-final-filt_anova-and-detp001_",
                    "2platforms_blood-all-groups.rda")
lcg <- get(load(lcg.fpath))
# format set names
names(lcg)[grepl("^peripheral.*", names(lcg))] <- "PBMC"

# define color palette
cbp2 <- rev(c("#000000", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

#-----------------
# stacked barplots
#-----------------
dfp <- do.call(rbind, lapply(seq(names(lcg)), function(i){
  groupi <- names(lcg)[i]; lgroupi <- lcg[[groupi]]
  num.shared <- length(lgroupi$overlap)
  which.epic.only <- !lgroupi$epic %in% lgroupi$overlap
  which.hm450k.only <- !lgroupi$hm450k %in% lgroupi$overlap
  num.epic.only <- length(lgroupi$epic[which.epic.only])
  num.hm450k.only <- length(lgroupi$hm450k[which.hm450k.only])
  dfpi <- data.frame(platform = c("HM450k-only", "Shared", "EPIC-only"), 
                     probes = as.numeric(c(num.hm450k.only, num.shared, 
                                           num.epic.only)), 
                     stringsAsFactors = F)
  dfpi$group <- groupi; return(dfpi)}))
dfp$platform <- factor(dfp$platform, levels = c("HM450k-only", 
                                                "Shared", "EPIC-only"))
colnames(dfp) <- c("Platform", "CpG probes", "Blood group")

bp.cgfilt <- ggplot(dfp, aes(x = `Blood group`, y = `CpG probes`, fill = Platform)) +
  geom_bar(stat = "identity") + theme_bw() + scale_fill_manual(values = cbp2) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 0.95))

bp.cgfilt.fpath <- "fig1b_bp-cgfilt_blood-groups_2platforms.pdf"
pdf(bp.cgfilt.fpath, 4, 3); print(bp.cgfilt); dev.off()

#----------------
# dodged barplots
#----------------
dfp <- do.call(rbind, lapply(seq(names(lcg)), function(i){
  groupi <- names(lcg)[i]; lgroupi <- lcg[[groupi]]
  dati <- unlist(lapply(lgroupi, function(x){length(x)}))
  dfpi <- data.frame(platform = names(dati), probes = as.numeric(dati), 
                     stringsAsFactors = F)
  dfpi$group <- groupi; return(dfpi)}))

bp.cgfilt <- ggplot(dfp, aes(x = group, y = probes, fill = platform)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = cbp2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.25, hjust = 0.95))

bp.cgfilt.fpath <- "bp-cgfilt_blood-groups_2platforms.pdf"
pdf(bp.cgfilt.fpath, 3, 2); print(bp.cgfilt); dev.off()

#----------------
# make upset plot
#----------------
pdf("fig1c_upset_cgfilt_blood-groups_2platforms.pdf", 4, 3.4)
upset(fromList(lcg), order.by = "freq", nsets = length(lcg),
      nintersects = 6, set_size.angles = 90, scale.sets = "log10",
      mainbar.y.label = "", sets.x.label = "Set size (log10)")
dev.off()

#----------------------
# get bval diffs as dfp
#----------------------
lstat.fname <- "lstat-cgfilt_bval-mean-var_2platforms_blood-groups.rda"
lstat <- get(load(lstat.fname))

dfp.cgdiff <- do.call(rbind, lapply(names(lstat), function(groupi){
  lgroupi <- lstat[[groupi]]
  df.hm450k.mean <- lgroupi$hm450k$mean; df.epic.mean <- lgroupi$epic$mean
  df.hm450k.var <- lgroupi$hm450k$var; df.epic.var <- lgroupi$epic$var
  names(df.hm450k.var) <- names(df.hm450k.mean)
  names(df.epic.var) <- names(df.epic.mean)
  cg.ol <- names(lgroupi$shared$mean)
  dff.hm450k.mean <- df.hm450k.mean[cg.ol];dff.epic.mean <- df.epic.mean[cg.ol]
  dff.hm450k.var <- df.hm450k.var[cg.ol];dff.epic.var <- df.epic.var[cg.ol]
  cond.mean <- identical(names(dff.hm450k.mean), names(dff.epic.mean))
  cond.var <- identical(names(dff.hm450k.var), names(dff.epic.var))
  if(!cond){stop("Couldn't match df's on cg.ol...")}
  dffv.mean <- dff.hm450k.mean - dff.epic.mean
  dffv.var <- dff.hm450k.var - dff.epic.var
  dfi <- data.frame(cgid = names(dffv.mean), mean.diff = as.numeric(dffv.mean),
                    var.diff = as.numeric(dffv.var), group = rep(groupi, length(dffv.mean)),
                    stringsAsFactors = F); return(dfi)}))

#------------------------------
# fig1d, bval diffs storm plots
#------------------------------
library(ggdist); library(gridExtra)

color <- "blue"; ptshape <- 16

dfp <- dfp.cgdiff
dfp[grepl("^peripheral.*", dfp$group),]$group <- "PBMC"

# single plot, mean diffs
sp.cgdiff <- ggplot(dfp, aes(x = group, y = mean.diff)) + theme_bw() + 
  ggdist::stat_halfeye(adjust = .3, width = 0.5, .width = 0, 
                       color = color, fill = color, justification = -.3, 
                       point_colour = color) + 
  geom_boxplot(width = .18, outlier.shape = NA, color = color) +
  gghalves::geom_half_point(side = "l", range_scale = 0, 
                            alpha = .1, color = color, shape = ptshape) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 0.95)) +
  ylab("Difference in Beta-value means\n(HM450K - EPIC)") +
  xlab("Blood group")

pdf.fname <- "fig1d_stormplot-cgdiff_2platforms_blood-groups.pdf"
pdf(pdf.fname, 3, 4); print(sp.cgdiff); dev.off()

# composite plot, mean and var diffs
sp.cgdiff.means <- ggplot(dfp, aes(x = group, y = mean.diff)) + theme_bw() + 
  ggdist::stat_halfeye(adjust = .3, width = 0.5, .width = 0, 
                       color = color, fill = color, justification = -.3, 
                       point_colour = color) + 
  geom_boxplot(width = .18, outlier.shape = NA, color = color) +
  gghalves::geom_half_point(side = "l", range_scale = 0, 
                            alpha = .1, color = color, shape = ptshape) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ylab("Mean difference")

sp.cgdiff.vars <- ggplot(dfp, aes(x = group, y = var.diff)) + theme_bw() + 
  ggdist::stat_halfeye(adjust = .3, width = 0.5, .width = 0, 
                       color = color, fill = color, justification = -.3, 
                       point_colour = color) + 
  geom_boxplot(width = .18, outlier.shape = NA, color = color) +
  gghalves::geom_half_point(side = "l", range_scale = 0, 
                            alpha = .1, color = color, shape = ptshape) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 0.95),
        axis.title.x = element_blank()) +
  ylab("Variance difference") + xlab("Blood group")

pdf.fname <- paste0("fig1d_stormplots-means-vars_cgdiffs-2platforms_blood-groups.pdf")
pdf(pdf.fname, 2, 4)
grid.arrange(sp.cgdiff.means, sp.cgdiff.vars, 
             layout_matrix = matrix(c(1,1,1,2,2,2,2), ncol = 1))
dev.off()

#-----------------------------------
# scatterplot, mean diff by var diff
#-----------------------------------

pdf("scatter-diffs.pdf", 15, 5)
ggplot(dfp, aes(x = mean.diff, y = var.diff, color = group)) +
  geom_point(alpha = 0.5) +
  facet_grid(mean.diff + group ~ var.diff)
dev.off()
