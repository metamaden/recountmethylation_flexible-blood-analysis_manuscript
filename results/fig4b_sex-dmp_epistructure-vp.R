#!/usr/bin/env R

# Author: Sean Maden
#
# Plots of EPISTRUCTURE components. Show top population genetic ancestry
# component values between independent samples from PBMC, whole blood, and 
# Inoshita et al 2015.
# 

library(ggplot2)
library(gridExtra)

#----------
# load data
#----------
save.dpath <- "."
md.fname <- "mdf-epic-hm450k-merge_normal-blood_0-0-3.rda"
md.fpath <- file.path(save.dpath, md.fname)
md <- get(load(file.path(md.fpath)))

#------------------
# make plot objects
#------------------
gseid <- "GSE67393"
md.filt <- md$gse == gseid | md$blood.subgroup == "whole_blood" | md$blood.subgroup == "PBMC"
mde <- md[md.filt,]; mde[mde$gse == gseid,]$blood.subgroup <- "Inoshita_et_al_2015"
mde$blood.subgroup <- factor(mde$blood.subgroup, levels = c("Inoshita_et_al_2015", "whole_blood", "PBMC"))

# make composite violin plots
lgg <- list()
# define the color palette
colv <- c("PBMC" = "#CD534CFF", "whole_blood" = "#868686FF", "Inoshita_et_al_2015" = "#3B3B3BFF")
lgg[["vp.epi.pc1"]] <- ggplot(mde, aes(x = blood.subgroup, y = glint.epi.pc1, fill = blood.subgroup)) + 
  geom_violin(color = "blue", draw_quantiles = 0.5) + theme_bw() + scale_fill_manual(values = colv) + 
  theme(legend.position = "none", axis.title.y = element_blank()) + ylab("PC1 (eigenvalue)") +
  coord_flip() + scale_y_continuous(limits = c(-40, 60))
lgg[["vp.epi.pc2"]] <- ggplot(mde, aes(x = blood.subgroup, y = glint.epi.pc2, fill = blood.subgroup)) + 
  geom_violin(color = "blue", draw_quantiles = 0.5) + theme_bw() + scale_fill_manual(values = colv) + 
  theme(legend.position = "none", axis.text.y = element_blank(),
        axis.title.y = element_blank()) + 
  ylab("PC2 (eigenvalue)") +
  coord_flip() + scale_y_continuous(limits = c(-40, 60))
mde$Subgroup <- mde$blood.subgroup

#--------------
# save new plot
#--------------
pdf.fname <- "fig3_vp-composite_sex-dmp-samples.pdf"
# layout matrix
lm <- matrix(c(1,1,1,1,1,1,1,1,1,1,
               2,2,2,2,2,2),nrow = 1)
# save new pdf
pdf(pdf.fname, 4.2, 1)
grid.arrange(lgg[["vp.epi.pc1"]], lgg[["vp.epi.pc2"]], layout_matrix = lm)
dev.off()