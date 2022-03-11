#!/usr/bin/env R

# Author: Sean Maden
#
# Make barplot summaries for studies and samples, with fill by blood subgroup.

library(ggplot2);library(gridExtra);library(ggpubr)

#-----------
# load data
#-----------
md.merge.fname <- "mdf-epic-hm450k-merge_normal-blood_0-0-3.rda"
md <- get(load(md.merge.fname))

#---------------
# get plot data
#---------------
dfp <- do.call(rbind, lapply(unique(md$blood.subgroup), function(gi){
  mf <- md[md$blood.subgroup==gi,]
  df <- do.call(rbind, lapply(unique(mf$platform), function(pi){
    mff <- mf[mf$platform==pi,]
    data.frame(platform = pi, subgroup = gi, ngsm = nrow(mff),
               ngse = length(unique(mff$gse)), stringsAsFactors = F)}))
  df <- rbind(df, data.frame(platform = "combined", subgroup = gi,
                             ngsm = nrow(mf), ngse = length(unique(mf$gse)),
                             stringsAsFactors = F)); return(df)}))
pal <- c("Other/NOS" = "#3B3B3BFF", "Whole blood" = "#868686FF", 
         "PBMC" = "#CD534CFF", "Cord blood" = "#EFC000FF")
dfp$platform.lab <- ifelse(dfp$platform == "combined", "Combined",
                           ifelse(dfp$platform == "epic", "EPIC", "HM450K"))
dfp$platform.lab <- factor(dfp$platform.lab, levels = c("Combined", "HM450K", "EPIC"))
dfp$Subgroup <- ifelse(dfp$subgroup == "other/NOS", "Other/NOS",
                       ifelse(dfp$subgroup == "whole_blood", "Whole blood",
                              ifelse(dfp$subgroup == "cord_blood", "Cord blood", "PBMC")))
dfp$Subgroup <- factor(dfp$Subgroup, levels = c("Other/NOS", "Whole blood", "PBMC", "Cord blood"))

#-----------------
# get plot objects
#-----------------
bp.gsm.val <- ggplot(dfp, aes(x = platform.lab, y = ngsm, fill = Subgroup)) + 
  geom_bar(stat = "identity", color = "#0073C2FF") + theme_bw() + scale_fill_manual(values = pal) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank()) + ylab("Samples")
bp.gsm.perc <- ggplot(dfp, aes(x = platform.lab, y = ngsm, fill = Subgroup)) + 
  geom_bar(stat = "identity", position = "fill", color = "#0073C2FF") + 
  theme_bw() + scale_fill_manual(values = pal) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ylab("Samples (%)")
bp.gse.val <- ggplot(dfp, aes(x = platform.lab, y = ngse, fill = Subgroup)) + 
  geom_bar(stat = "identity", color = "#0073C2FF") + theme_bw() + scale_fill_manual(values = pal) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_blank()) + ylab("Studies\n")
bp.gse.perc <- ggplot(dfp, aes(x = platform.lab, y = ngse, fill = Subgroup)) + 
  geom_bar(stat = "identity", position = "fill", color = "#0073C2FF") + 
  theme_bw() + scale_fill_manual(values = pal) +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ylab("Studies (%)")
plot.legend <- ggplot(dfp, aes(x = platform.lab, y = ngsm, fill = Subgroup)) + 
  geom_bar(stat = "identity", color = "#0073C2FF") + scale_fill_manual(values = pal)
plot.legend <- get_legend(plot.legend)

#--------------------
# save composite plot
#--------------------
fig1.fname <- "fig1_bp-gsm-gse_blood-subgroups_2platforms.pdf"
lm <- matrix(c(rep(c(1,1,3,3,5), 12), 
               rep(c(2,2,4,4,5), 18)), 
             byrow = T, nrow = 30)
xlab.str <- paste0("Platform", paste0(rep(" ", 18), collapse = ""), collapse = "")

# save pdf
pdf(fig1.fname, 5.8, 2.5)
grid.arrange(bp.gsm.val, bp.gsm.perc, bp.gse.val, bp.gse.perc, 
             plot.legend, layout_matrix = lm, bottom = xlab.str)
dev.off()