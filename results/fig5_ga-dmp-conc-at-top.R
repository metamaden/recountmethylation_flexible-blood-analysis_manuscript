#!/usr/bin/env R

# Author: Sean Maden
#
# Make the concordance at the top plot for getational age DMPs.

#----------
# load data
#----------
st.fname <- "st_gest-age-dmp_all-dmp-info.rda"
dmpdf <- get(load(file.path(st.fname)))

#--------------
# get plot data
#--------------
cnv <- colnames(dmpdf)[2:3]
dmp.cb <- dmpdf[dmpdf$is.cord.blood.dmp==T,]
dmp.cb <- dmp.cb[order(dmp.cb$cord.blood.ttest.pnominal),]
dfp <- do.call(rbind, lapply(seq(nrow(dmp.cb)), function(si){
  dmpf <- dmp.cb[seq(si),]
  do.call(rbind, lapply(cnv, function(ci){
    data.frame(dmp.index = si, 
               num.dmp = nrow(dmpf[dmpf[,ci]==T,]), 
               dmp.set = ci)
  }))
}))

# get plot object
dfp$Study <- ifelse(dfp$dmp.set == "is.haftorn.2021.dmp", "Haftorn et al 2021", "Merid et al 2020")
col.merid <- "#26AD34"; col.haftorn <- "#B570DF"
lab.haftorn <- max(dfp[dfp$dmp.set == "is.haftorn.2021.dmp",]$num.dmp)
lab.merid <- max(dfp[dfp$dmp.set == "is.merid.2020.dmp",]$num.dmp)
colv <- c("Merid et al 2020" = col.merid, "Haftorn et al 2021" = col.haftorn)

gg.cat <- ggplot(dfp, aes(x = dmp.index, y = num.dmp, color = Study)) +
  geom_point(alpha = 0.1, cex = 0.5) + geom_line() + theme_bw() + scale_color_manual(values = colv) +
  geom_hline(yintercept = lab.merid, color = col.merid, linetype = "longdash") +
  geom_hline(yintercept = lab.haftorn, color = col.haftorn, linetype = "longdash") +
  annotate("text", color = col.merid, x = 50, y = lab.merid - 80, label= lab.merid) +
  annotate("text", color = col.haftorn, x = 20, y = lab.haftorn + 80, label= lab.haftorn) +
  xlab("Total DMPs (cord blood)") + ylab("Validated DMPs")

#--------------
# save new plot
#--------------
pdf.fname <- "fig5_ga-dmp-conc-at-top.pdf"
pdf(pdf.fname, 4.2, 1.5); print(gg.cat); dev.off()