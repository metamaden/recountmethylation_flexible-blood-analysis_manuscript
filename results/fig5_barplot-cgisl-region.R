#!/usr/bin/env R

# Author: Sean Maden
#
# Barplots of CpG Island-relation regions across DMP sets.

library(ggplot2)
library(minfiData)
library(minfiDataEPIC)

#-----------
# load data
#-----------
st.fname <- "st_gest-age-dmp_all-dmp-info.rda"
dmpdf <- get(load(file.path(st.fname)))

# load the probe annotations for background
data(RGsetEPIC)
anno.bg <- getAnnotation(RGsetEPIC)
anno.bg <- anno.bg[anno.bg$Methyl450_Loci==T,]

#--------------
# get plot data
#--------------
cnv <- colnames(dmpdf)[2:4]
isl.typev <- c("Island", "Shore", "Shelf", "OpenSea")
dmp.labv <- c("Background", "cord_blood", 
              "Merid et al 2020", "Haftorn et al 2021")
dfp <- do.call(rbind, lapply(isl.typev, function(isl.type){
  # background fraction
  num.region.bg <- nrow(anno.bg[grepl(isl.type, anno.bg$Relation_to_Island),])
  fract.region.bg <- num.region.bg/nrow(anno.bg)
  # dmp fractions
  fract.dmpv <- unlist(lapply(cnv, function(ci){
    dmpdff <- dmpdf[dmpdf[,ci]==T,]
    is.region <- grepl(isl.type, dmpdff$relation.to.island)
    num.region <- length(which(is.region))
    num.region/nrow(dmpdff)}))
  data.frame(region = rep(isl.type, 4), dmp_set = dmp.labv,
                     fract_region = c(fract.region.bg, fract.dmpv),
                     stringsAsFactors = FALSE)}))
dfp$`DMP set` <- dfp$dmp_set

#-----------------
# make plot object
#-----------------
dfp$hline <- rep(dfp[dfp$dmp_set=="Background",]$fract_region, each = 4)
# get plot object
ggbp <- ggplot(dfp, aes(x = dmp_set, y = fract_region, fill = `DMP set`)) + 
  geom_bar(stat = "identity") + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(), legend.position = "none") + 
  ylab("DMP fraction") +
  geom_hline(df = dfp, aes(yintercept = hline), 
             size = 0.8, color = "red", alpha = 0.4)
bpfinal <- ggbp + facet_wrap(~region)

#--------------
# save new plot
#--------------
pdf.fname <- "barplot_islregions-anno_ga-dmps-3sets.pdf"
pdf(pdf.fname, 3, 2.8); print(bpfinal); dev.off()