#!/usr/bin/env R

# Author: Sean Maden
#
# Gestational age DMP genes upset plot.

library(UpSetR)

#----------
# load data
#----------
st.fname <- "st_gest-age-dmp_all-dmp-info.rda"
dmpdf <- get(load(file.path(st.fname)))

#----------------------
# get plot data as list
#----------------------
dmpdff <- dmpdf[!dmpdf$ucsc.gene.names=="",]
namev <- c("Haftorn_et_al_2021", "Merid_et_al_2020", "cord_blood"); cnv <- colnames(dmpdff)[2:4]
lgene <- lapply(cnv, function(ci){
  unique(unlist(strsplit(dmpdff[dmpdff[,ci]==T,]$ucsc.gene.names, ";")))})
names(lgene) <- namev

#--------------------
# save the upset plot
#--------------------
# get plot vars
pdf.fname <- "fig5_upset-ga-dmp-genes.pdf"
ylab.str <- paste0(paste0(rep("\n", 10),collapse =""), 
                   "Overlapping DMP genes")
xlab.str <- "Total DMP genes"

# save new pdf
pdf(pdf.fname, 3.6, 2.6)
upset(fromList(lgene), order.by = "freq", 
      mainbar.y.label = ylab.str, sets.x.label = xlab.str)
dev.off()