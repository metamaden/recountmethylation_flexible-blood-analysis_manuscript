#!/usr/bin/env R

# Author: Sean Maden
#
# Plots of the search index nearest neighbors querying Inoshita et al 2015
# whole blood samples. Nearest neighbors were obtained from querying an 
# autosomal DNAm search index of normal blood DNAm array samples run on 
# either the HM450K or EPIC platforms. Metadata label frequencies were then 
# plotted among nearest neighbors returned.
#

library(ggplot2)
library(ggdist)

#----------
# load data
#----------
# load search index results
sit.fpath <- "sitest_results_blood-groups-2platforms.csv"
sit <- data.table::fread(sit.fpath, sep = ",", header = T, data.table = F)
sit <- sit[,c(2:ncol(sit))]
colnames(sit) <- c("gsm", "k=1000", "platform")

# get metadata labels
save.dpath <- "." 
md.fname <- "mdf-epic-hm450k-merge_normal-blood_0-0-3.rda"
md.fpath <- file.path(save.dpath, md.fname)
md <- get(load(file.path(md.fpath)))

# get color palette
colv <- c("PBMC" = "#CD534CFF", "whole_blood" = "#868686FF", "other/NOS" = "#3B3B3BFF")

#---------------
# get si results
#---------------
# subset si results
gsmv <- md$gsm
sif <- sit[sit[,1] %in% gsmv,]
# get fract group by pos
md$group <- md$blood.subgroup
md[md$group %in% c("other/NOS", "cord_blood"),]$group <- "other/NOS"

#------------------------
# get raincloud plot data
#------------------------
# get plot data
rownames(md) <- md$gsm
mpos <- do.call(rbind, lapply(seq(nrow(sif)), function(ii){
  labv <- unlist(strsplit(sif[ii,2], ";"))
  gsmv <- gsub("\\..*", "", labv)
  md[gsmv,]$group}))
dfp <- do.call(rbind, lapply(unique(md$group), function(groupi){
  num.groupi <- apply(mpos, 1, function(ri){length(ri[ri==groupi])})
  data.frame(subgroup = rep(groupi, length(num.groupi)),
             num.subgroup = num.groupi, stringsAsFactors = F)}))
# format dfp vars
dfp$subgroup <- factor(dfp$subgroup, 
                       levels = c("whole_blood", "other/NOS", "PBMC"))


# get new plot object
text.size <- 15
title.size <- 18
ptsize <- 2
linesize <- 0.7
ggrain <- ggplot(dfp, aes(x = subgroup, y = num.subgroup, 
                          fill = "white", colour = subgroup)) + 
  stat_halfeye(aes(fill = subgroup), width = 5, 
               justification = -0.06, point_colour = NA, .width = 0) + 
  geom_boxplot(lwd = linesize, width = 0.45, alpha = 1, 
               lwd = 0.35, outlier.size = ptsize, fill = "white",
               outlier.alpha = 0.1, outlier.shape = 1) + 
  theme_bw() +  guides(fill = "none") + scale_color_manual(values = colv) +
  scale_fill_manual(values = colv) + ylab("Samples") +
  coord_flip() + theme(legend.position = "none") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = text.size),
        axis.text.x = element_text(size = text.size, angle = 45, 
                                   hjust = 1, vjust = 1),
        axis.title.x = element_text(size = title.size))

# save new pdf
pdf.fname <- "fig4_sex-dmp_nn1k-raincloud.pdf"
pdf(pdf.fname, 4, 3.5)
print(ggrain)
dev.off()