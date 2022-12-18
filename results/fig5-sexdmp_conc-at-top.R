#!/usr/bin/env R

# Author: Sean Maden
#
# Make composite concordance at the top plot for sex differential DNAm.
# In short, show the sex DMPs from Inoshita et al 2015 which are validated
# among probes ordered on their significance/magnitude of differential 
# DNAm between sexes for PBMC and cord blood compilations.
#

library(ggsci); library(patchwork); library(png)
library(ggplot2)

#----------
# load data
#----------
# ttest validation results
tdf.wb.fname <- "ttest-df-results-mvalfit_wholeblood-2platforms_inoshita-2015-validate.rda"
tdf.pbmc.fname <- "ttest-df-results-mvalfit_pbmc-2platforms_inoshita-2015-validate.rda"
tdf.wb <- get(load(tdf.wb.fname))
tdf.pbmc <- get(load(tdf.pbmc.fname))
# study dmps
dmp.study.fname <- "table1_sex-dmp_inoshita-2015.csv"
dmp.study <- read.csv(dmp.study.fname)

#--------------
# get plot data
#--------------
# get concordance at top data
colv <- c("PBMC" = "#CD534CFF", "whole_blood" = "#868686FF")
ltdf <- list("whole_blood" = tdf.wb, "PBMC" = tdf.pbmc)
num.dmp <- c(seq(0, 950, 50), seq(1000, nrow(tdf), 1000))
dfp.cat <- do.call(rbind, lapply(seq(2), function(ii){
  tdf <- ltdf[[ii]]; groupi <- names(ltdf)[ii]
  tdf <- tdf[order(as.numeric(unlist(tdf$ttest.pnom))),]
  num.validated <- sapply(num.dmp, function(x){
    length(intersect(tdf[c(1:x),1], dmp.study[,2]))})
  dfp <- data.frame(num.dmp = num.dmp, 
                    num.validated = num.validated,
                    stringsAsFactors = F)
  dfp$color <- colv[groupi]; dfp$Subgroup <- groupi
  return(dfp)
}))

#------------------------------------------
# summarize mean diffs, pvals by probe sets
#------------------------------------------
# format vars
for(c in c(2,3,4)){
  tdf.wb[,c] <- as.numeric(tdf.wb[,c])
  tdf.pbmc[,c] <- as.numeric(tdf.pbmc[,c])}

# replace pval = 0
pmin.wb <- min(tdf.wb[tdf.wb$ttest.pnom>0,]$ttest.pnom)
pmin.pbmc <- min(tdf.pbmc[tdf.pbmc$ttest.pnom>0,]$ttest.pnom)
tdf.wb[tdf.wb$ttest.pnom==0,]$ttest.pnom <- pmin.wb
tdf.pbmc[tdf.pbmc$ttest.pnom==0,]$ttest.pnom <- pmin.pbmc
# order on pval
tdf.wb <- tdf.wb[order(tdf.wb$ttest.pnom),]
tdf.pbmc <- tdf.pbmc[order(tdf.pbmc$ttest.pnom),]
# get adjusted pvalues
tdf.wb$pbh <- p.adjust(tdf.wb$ttest.pnom, method = "BH")
tdf.pbmc$pbh <- p.adjust(tdf.pbmc$ttest.pnom, method = "BH")
# get -1 * log10 pbh
tdf.wb$yaxis <- -1*log10(tdf.wb$pbh + min(tdf.wb$pbh))
tdf.pbmc$yaxis <- -1*log10(tdf.pbmc$pbh)
# get diffs
tdf.wb$sex.diff <- tdf.wb$mean.male-tdf.wb$mean.female
tdf.pbmc$sex.diff <- tdf.pbmc$mean.male-tdf.pbmc$mean.female

# base r plots
# volcano plots
pdf("volcano_wb.pdf", 4, 5)
plot(tdf.wb$sex.diff, tdf.wb$yaxis, main = "Whole blood",
     col = rgb(0,0,0,0.4), xlab = "Mean diff. (M - F)",
     ylab = "-1*log10(P-adj.)")
points(tdf.wb[seq(1000),]$sex.diff, tdf.wb[seq(1000),]$yaxis, 
       pch = 16, col = rgb(1, 0, 0, 1))
abline(v = 0, col = "blue")
dev.off()
pdf("volcano_pbmc.pdf", 4, 5)
plot(tdf.pbmc$sex.diff, tdf.pbmc$yaxis, main = "PBMC",
     col = rgb(0,0,0,0.4), xlab = "Mean diff. (M - F)",
     ylab = "-1*log10(P-adj.)")
points(tdf.pbmc[seq(1000),]$sex.diff, tdf.pbmc[seq(1000),]$yaxis, 
       pch = 16, col = rgb(1, 0, 0, 1))
abline(v = 0, col = "blue")
dev.off()

# low-res volcano plots
jpeg("volcano_wb.jpg", width = 4, height = 5, units = "in", res = 300)
plot(tdf.wb$sex.diff, tdf.wb$yaxis, main = "Whole blood",
     col = rgb(0,0,0,0.4), xlab = "Mean diff. (M - F)",
     ylab = "-1*log10(P-adj.)")
points(tdf.wb[seq(1000),]$sex.diff, tdf.wb[seq(1000),]$yaxis, 
       pch = 16, col = rgb(1, 0, 0, 1))
abline(v = 0, col = "blue")
dev.off()
pdf("volcano_pbmc.jpg", width = 4, height = 5, units = "in", res = 300)
plot(tdf.pbmc$sex.diff, tdf.pbmc$yaxis, main = "PBMC",
     col = rgb(0,0,0,0.4), xlab = "Mean diff. (M - F)",
     ylab = "-1*log10(P-adj.)")
points(tdf.pbmc[seq(1000),]$sex.diff, tdf.pbmc[seq(1000),]$yaxis, 
       pch = 16, col = rgb(1, 0, 0, 1))
abline(v = 0, col = "blue")
dev.off()

# scatterplots -- means, male vs. female
pdf("scatter-diff_wb.pdf", 5, 5)
plot(tdf.wb$mean.male, tdf.wb$mean.female, main = "Whole blood")
points(tdf.wb[seq(1000),]$mean.male, 
       tdf.wb[seq(1000),]$mean.female, 
       col = "red", pch = 16)
abline(a = 0, b = 1, col = "blue")
dev.off()
pdf("scatter-diff_pbmc.pdf", 5, 5)
plot(tdf.pbmc$mean.male, tdf.pbmc$mean.female, main = "PBMC")
points(tdf.pbmc[seq(1000),]$mean.male, 
       tdf.pbmc[seq(1000),]$mean.female, 
       col = "red", pch = 16)
abline(a = 0, b = 1, col = "blue")
dev.off()

# abs diff cutoff -- wb
dfp <- tdf.wb
dfp$cat <- "all"
dfp.top <- tdf.wb[seq(1000),]
dfp.top$cat <- "top_1k"
dfp <- rbind(dfp, dfp.top)
pdf("violin-absdiff_wb.pdf", width = 2, height = 2)
ggplot(dfp, aes(x = cat, y = abs.diff)) + ggtitle("Whole blood") + theme_bw() + 
  geom_violin(draw_quantiles = 0.5) + ylab("Mean Abs. diff.") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, 
                                   vjust = 1))
dev.off()

# abs diff cutoff -- pbmc
dfp <- tdf.pbmc
dfp$cat <- "all"
dfp.top <- tdf.pbmc[seq(1000),]
dfp.top$cat <- "top_1k"
dfp <- rbind(dfp, dfp.top)
pdf("violin-absdiff_pbmc.pdf", width = 2, height = 2)
ggplot(dfp, aes(x = cat, y = abs.diff)) + ggtitle("PBMC") + theme_bw() + 
  geom_violin(draw_quantiles = 0.5) + ylab("Mean Abs. diff.") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, 
                                   vjust = 1))
dev.off()


#-----------------
# get plot objects
#-----------------
ndmp.tot.pbmc <- max(dfp.cat[dfp.cat$Subgroup == 'PBMC',]$num.validated)
ndmp.tot.wb <- max(dfp.cat[dfp.cat$Subgroup == 'whole_blood',]$num.validated)
ymax.zoom <- 100; xmax.zoom <- 1e3; ymin.zoom <- xmin.zoom <- 0

# main scatterplot
text.size.main <- 21
title.size <- 25
text.size.inset <- 20
label.size <- 8
ptsize <- 2
linesize <- 2
plot.main <- ggplot(dfp.cat, aes(x = num.dmp, y = num.validated, 
                                 group = Subgroup, color = Subgroup)) +
  annotate("rect", xmin = -15000, xmax = 20000, ymin = 0, ymax = 100, 
           alpha = 0.45, fill = "#a9a9a9") + scale_color_manual(values = colv) +
  ylab("Validated DMPs") + xlab("Total DMPs") + 
  geom_hline(aes(size = linesize), yintercept = ndmp.tot.pbmc, color = "#CD534CFF", 
             linetype = "longdash") +
  geom_hline(aes(size = linesize), yintercept = ndmp.tot.wb, color = "#868686FF", 
             linetype = "longdash") +
  annotate("text", color = "#CD534CFF", x = 1300, y = 288, 
           label= ndmp.tot.pbmc, size = label.size) +
  annotate("text", color = "#868686FF", x = 1300, y = 247, 
           label= ndmp.tot.wb, size = label.size) +
  geom_point(size = ptsize) + geom_line() + theme_bw() + ylim(0, 300) +
  annotation_raster(readPNG("magnifying_glass_bgtransparent.png"), 
                    xmin = 8000, xmax = 58000, ymin = 20, ymax = 80) +
  annotation_raster(readPNG("rightarrow_bgtransparent.png"), 
                    xmin = 58500, xmax = 135000, ymin = 20, ymax = 80) +
  theme(axis.text.x = element_text(size = text.size.main, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = text.size.main),
        axis.title.x = element_text(size = title.size),
        axis.title.y = element_text(size = title.size),
        legend.text = element_text(size = text.size.main),
        legend.title = element_text(size = title.size))
# zoomed inset
dfp.cat$Subgroup <- factor(dfp.cat$Subgroup, levels = c("PBMC", "whole_blood"))
plot.inset <- ggplot(dfp.cat, aes(x = num.dmp, y = num.validated, 
                                  group = Subgroup, color = Subgroup)) +
  geom_point(size = ptsize) + geom_line(aes(size = linesize)) + geom_line() + theme_bw() + 
  scale_color_manual(values = colv) +
  theme(axis.text.x = element_text(size = text.size.inset, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = text.size.inset),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ylim(ymin.zoom, ymax.zoom) + xlim(xmin.zoom, xmax.zoom)

# final composite plot object
plot.composite <- plot.main + 
  inset_element(plot.inset, left = 0.45, 
                bottom = 0.009, right = 0.98, 
                top = 0.68)

#--------------------
# save composite plot
#-------------------- 
# plot concodrance at the top data
pdf.fname <- "fig4_sex-dmp_conc-at-top-comp.pdf"
pdf(pdf.fname, width = 7.6, height = 4.3)
print(plot.composite)
dev.off()