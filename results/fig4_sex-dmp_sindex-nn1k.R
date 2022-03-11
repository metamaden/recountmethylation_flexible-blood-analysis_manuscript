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

# get new plot object
ggrain <- ggplot(dfp, aes(x = subgroup, y = num.subgroup, 
                          fill = "white", colour = subgroup)) + 
  stat_halfeye(aes(fill = subgroup), width = 4, 
               justification = -0.08, point_colour = NA, .width = 0) + 
  geom_boxplot(width = 0.5, alpha = 1, lwd = 0.35, 
               outlier.size = 1.1, outlier.alpha = 0.1,
               outlier.shape = 1,
               fill = "white") + 
  theme_bw() +  guides(fill = "none") + scale_color_manual(values = colv) +
  scale_fill_manual(values = colv) + ylab("Number of samples (k = 1,000 neighbors)") +
  coord_flip() + theme(legend.position = "none") +
  theme(axis.title.y = element_blank())

# save new pdf
pdf.fname <- "fig3_sex-dmp_nn1k-raincloud.pdf"
pdf(pdf.fname, 4.2, 2.8)
print(ggrain)
dev.off()