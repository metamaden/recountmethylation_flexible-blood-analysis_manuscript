#!/usr/bin/env R

# Author: Sean Maden
#
# Make supplemental figures, volcano plots and violin plots of whole blood and PBMC 
# sex DMP analysis results.
#

libv <- c("ggplot2")
sapply(libv, library, character.only = TRUE)

#----------
# load data
#----------
# ttest validation results
tdf.wb.fname <- "ttest-df-results-mvalfit_wholeblood-2platforms_inoshita-2015-validate.rda"
tdf.pbmc.fname <- "ttest-df-results-mvalfit_pbmc-2platforms_inoshita-2015-validate.rda"
tdf.wb <- get(load(tdf.wb.fname))
tdf.pbmc <- get(load(tdf.pbmc.fname))

#-------------------
# get plot variables
#-------------------
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
# get abs.diff
tdf.wb$abs.diff <- abs(tdf.wb$sex.diff)
tdf.pbmc$abs.diff <- abs(tdf.pbmc$sex.diff)

#------------------
# new volcano plots
#------------------
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

jpeg("volcano_pbmc.jpg", width = 4, height = 5, units = "in", res = 300)
plot(tdf.pbmc$sex.diff, tdf.pbmc$yaxis, main = "PBMC",
     col = rgb(0,0,0,0.4), xlab = "Mean diff. (M - F)",
     ylab = "-1*log10(P-adj.)")
points(tdf.pbmc[seq(1000),]$sex.diff, tdf.pbmc[seq(1000),]$yaxis, 
       pch = 16, col = rgb(1, 0, 0, 1))
abline(v = 0, col = "blue")
dev.off()

#----------------------
# save new violin plots
#----------------------
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


ggplot(dfp, aes(x = cat, y = abs.diff)) + geom_violin(draw_quantiles = 0.5) +
  theme_bw() + ggtitle("PBMC") + ylab("Mean Abs. diff.") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))