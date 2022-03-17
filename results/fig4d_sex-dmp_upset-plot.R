#!/usr/bin/env R

# Author: Sean Maden
#
# Upset plot for sex DMP sets. Get set overlaps for 4 independent sex DMP
# sets, from: (1) Inoshita et al 2015; (2) whole blood compilation; 
# (3) PBMC compilation; and (4) Grant et al 2021.
# 

library(UpSetR)

#----------
# load data
#----------
# grant.fpath <- file.path("grant_2021_tables", "media-1.csv")
# grant.dmp <- read.csv(grant.fpath)
grant.dmp <- read.csv("grant-2021_wb-sex-dmp.csv")
# load dmp info
new.dmp <- get(load(file.path("st_sex-dmp_all-dmp-info.rda")))

#-----------------
# make probes list
#-----------------
lup <- list()
lup[["Inoshita_et_al_2015"]] <- new.dmp[new.dmp$is.inoshita.2015.dmp==T,]$cgid
lup[["whole_blood"]] <- new.dmp[new.dmp$is.whole.blood.dmp==T,]$cgid
lup[["PBMC"]] <- new.dmp[new.dmp$is.pbmc.dmp==T,]$cgid
lup[["Grant_et_al_2021"]] <- grant.dmp$X

#--------------
# save new plot
#--------------
# get plot vars
pdf.fname <- "fig4d_upset-sexdmp-4sets.pdf"
xlab.str <- "Total DMPs"
ylab.str <- paste0(paste0(rep('\n', 10),collapse = ""),
                   "Overlapping DMPs", collapse = "")

# save new pdf
pdf(pdf.fname, 5.8, 3.3)
print(upset(fromList(lup), order.by = "freq", 
            text.scale = 1.4,set_size.angles = 45,
            set_size.numbers_size = 1,
            mainbar.y.label = ylab.str,
            sets.x.label = xlab.str))
dev.off()