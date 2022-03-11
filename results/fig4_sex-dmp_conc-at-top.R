#!/usr/bin/env R

# Author: Sean Maden
#
# Make composite concordance at the top plot for sex differential DNAm.
# In short, show the sex DMPs from Inoshita et al 2015 which are validated
# among probes ordered on their significance/magnitude of differential 
# DNAm between sexes for PBMC and cord blood compilations.
#

library(ggsci); library(patchwork); library(png)

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

#-----------------
# get plot objects
#-----------------
ndmp.tot.pbmc <- max(dfp.cat[dfp.cat$Subgroup == 'PBMC',]$num.validated)
ndmp.tot.wb <- max(dfp.cat[dfp.cat$Subgroup == 'whole_blood',]$num.validated)
ymax.zoom <- 100; xmax.zoom <- 1e3; ymin.zoom <- xmin.zoom <- 0

# main scatterplot
plot.main <- ggplot(dfp.cat, aes(x = num.dmp, y = num.validated, 
                                 group = Subgroup, color = Subgroup)) +
  annotate("rect", xmin = -8000, xmax = 15000, ymin = 0, ymax = 100, 
           alpha = 0.45, fill = "grey") + scale_color_manual(values = colv) +
  ylab("Validated DMPs") + xlab("Total DMPs") + theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = ndmp.tot.pbmc, color = "#CD534CFF", linetype = "longdash") +
  geom_hline(yintercept = ndmp.tot.wb, color = "#868686FF", linetype = "longdash") +
  annotate("text", color = "#CD534CFF", x = 1000, y = 285, label= ndmp.tot.pbmc) +
  annotate("text", color = "#868686FF", x = 1000, y = 250, label= ndmp.tot.wb) +
  geom_point() + geom_line() + theme_bw() + ylim(0, 300) +
  annotation_raster(readPNG("magnifying_glass_bgtransparent.png"), 
                    xmin = 8000, xmax = 48000, ymin = 20, ymax = 80) +
  annotation_raster(readPNG("rightarrow_bgtransparent.png"), 
                    xmin = 53000, xmax = 110000, ymin = 20, ymax = 80)
# zoomed inset
dfp.cat$Subgroup <- factor(dfp.cat$Subgroup, levels = c("PBMC", "whole_blood"))
plot.inset <- ggplot(dfp.cat, aes(x = num.dmp, y = num.validated, 
                                  group = Subgroup, color = Subgroup)) +
  geom_point() + geom_line() + geom_line() + theme_bw() + 
  scale_color_manual(values = colv) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ylim(ymin.zoom, ymax.zoom) + xlim(xmin.zoom, xmax.zoom)

# final composite plot object
plot.composite <- plot.main + 
  inset_element(plot.inset, left = 0.3, bottom = 0.02, 
                right = 0.98, top = 0.65)

#--------------------
# save composite plot
#-------------------- 
# plot concodrance at the top data
pdf.fname <- "fig3_sex-dmp_conc-at-top-comp.pdf"
pdf(pdf.fname, 6, 3)
print(plot.composite)
dev.off()