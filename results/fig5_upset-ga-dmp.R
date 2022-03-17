#!/usr/bin/env R

# Author: Sean Maden
#
# Gestational age DMP upset plot.

library(UpSetR)

#----------
# load data
#----------
st.fname <- "st_gest-age-dmp_all-dmp-info.rda"
dmpdf <- get(load(file.path(st.fname)))

#--------------
# get plot data
#--------------
cnv <- colnames(dmpdf)[2:4]
ldmp <- lapply(cnv,function(ci){dmpdf[dmpdf[,ci]==TRUE,1]})
names(ldmp) <- c("Haftorn_et_al_2021", "Merid_et_al_2021", "cord_blood")

#--------------------
# save new upset plot
#--------------------
# make dmp set upset plot -- figure 4a
pdf.fname <- "upset_gest-age-dmps_3-studies.pdf"
ylab.str <- paste0(paste0(rep("\n", 8),collapse =""), 
                   "Overlapping DMPs")
xlab.str <- "Total DMPs"
  
# plot results
pdf(pdf.fname, 5.3, 2.7)
upset(fromList(ldmp), order.by = "freq",
      mainbar.y.label = ylab.str, sets.x.label = xlab.str, 
      text.scale = 1.1)
dev.off()

