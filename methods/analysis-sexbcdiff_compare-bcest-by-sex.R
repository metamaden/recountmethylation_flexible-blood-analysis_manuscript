#!/usr/bin/env R

# Author: Sean Maden
#
# Compare estimated blood cell proportions between males and females.
#
#

libv <- c("HDF5Array", "minfi", "ggplot2")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
# load grset
gr.fpath <- file.path("gr-gseadj_h5se_hm450k-epic-merge_0-0-3")
gr <- loadHDF5SummarizedExperiment(gr.fpath)

# compare blood cell abundances
pd <- pData(gr)
cnv <- colnames(pd)
bcv <- cnv[grepl("predcell", cnv)]

#-----------------
# helper functions
#-----------------
get_dfstat <- function(dfp, ctv, rounding = FALSE, csv.fname = NULL){
  # make df of summary statistics
  dfstat <- do.call(rbind, lapply(ctv2, function(cti){
    dff <- dfp[dfp$celltype==cti,]
    grp1 <- dff[dff$sex=="M",]$cellprop; grp2 <- dff[dff$sex=="F",]$cellprop
    wti <- wilcox.test(x = grp1, y = grp2, alternative = "two.sided")
    tti <- t.test(x = grp1, y = grp2, alternative = "two.sided")
    data.frame(celltype = cti, 
               blood.subgroup = ifelse("blood.subgroup" %in% colnames(dfp), 
                                       unique(dfp$blood.subgroup), "NA"),
               mean.male = mean(grp1), median.male = median(grp1),
               sd.male = sd(grp1), var.male = var(grp1),
               mean.female = mean(grp2), median.female = median(grp2),
               sd.female = sd(grp2), var.female = var(grp2),
               wcox.stat = wti$statistic, wcox.punadj = wti$p.value,
               ttest.stat = tti$statistic, ttest.punadj = tti$p.value)
  }))
  # get adj pvalues
  dfstat$wcox.padj.bh <- p.adjust(dfstat$wcox.punadj, method = "BH")
  dfstat$ttest.padj.bh <- p.adjust(dfstat$ttest.punadj, method = "BH")
  if(rounding){
    # format large num
    for(c in c(3:10, 12)){
      dfstat[,c] <- round(dfstat[,c], 3)
    } 
    # format small num
    for(c in c(12,14:16)){
      dfstat[,c] <- format(dfstat[,c], digits = 3) 
    }
  }
  if(!is(csv.fname, "NULL")){
    write.csv(dfstat, file = csv.fname, row.names = F)}
  return(dfstat)
}

#---------------------------
# summarize, stats and tests
#---------------------------
# get dfstat for each blood sample type
pd.all <- pd; pd.all$blood.subgroup <- "all"
pd <- rbind(pd, pd.all)
grpv <- unique(pd$blood.subgroup)

# iterate on groups
dfstat <- do.call(rbind, lapply(grpv, function(grpi){
  pdf <- pd[pd$blood.subgroup==grpi,]
  # get plot df
  dfp <- do.call(rbind, lapply(bcv, function(bi){
    dfi <- as.data.frame(pdf[,c("predsex", cnv[grepl(bi, cnv)])])
    dfi$type <- gsub("predcell\\.", "", bi)
    colnames(dfi) <- c("sex", "cellprop", "celltype")
    dfi[,2] <- as.numeric(dfi[,2]); dfi
  }))
  dfp$blood.subgroup <- grpi
  get_dfstat(dfp, unique(dfp$celltype))
}))
# save dfstat
save(dfstat, file = paste0(dfstat.fname, ".rda"))

#-------------------------------------------
# summarize, stats and tests -- use rounding
#-------------------------------------------
# iterate on groups -- rounded
dfstat.round <- do.call(rbind, lapply(grpv, function(grpi){
  pdf <- pd[pd$blood.subgroup==grpi,]
  # get plot df
  dfp <- do.call(rbind, lapply(bcv, function(bi){
    dfi <- as.data.frame(pdf[,c("predsex", cnv[grepl(bi, cnv)])])
    dfi$type <- gsub("predcell\\.", "", bi)
    colnames(dfi) <- c("sex", "cellprop", "celltype")
    dfi[,2] <- as.numeric(dfi[,2]); dfi
  }))
  dfp$blood.subgroup <- grpi
  get_dfstat(dfp, unique(dfp$celltype), rounding = T)
}))

# write csv -- rounded
dfstat.round.fname <- "supptable-rounded_sexdiff-blood-cell-est"
write.csv(dfstat.round, file = paste0(dfstat.round.fname, ".csv"), row.names = F)

#-------------
# violin plots
#-------------
# get plot df
dfp <- do.call(rbind, lapply(grpv, function(grpi){
  pdf <- pd[pd$blood.subgroup==grpi,]
  dfpi <- do.call(rbind, lapply(bcv, function(bi){
    dfi <- as.data.frame(pdf[,c("predsex", cnv[grepl(bi, cnv)])])
    dfi$type <- gsub("predcell\\.", "", bi)
    colnames(dfi) <- c("sex", "cellprop", "celltype")
    dfi[,2] <- as.numeric(dfi[,2]); dfi
  }))
  dfpi$blood.subgroup <- grpi
  dfpi
}))

# make plot object
ggvp <- ggplot(dfp, aes(x = sex, y = cellprop, color = sex)) +
  geom_violin(draw_quantiles = 0.5) + theme_bw() + 
  ylab("Proportion") + xlab("Sex")
ggvp <- ggvp + facet_wrap(~blood.subgroup + celltype, nrow = 5)

# save new plot
vp.fname <- "ggvp-groups_sexdiff-blood-cell.pdf"
pdf(vp.fname, width = 7, height = 7); ggvp; dev.off()

#-----------------------
# scatter plots of means
#-----------------------
# make plot object
ggpt <- ggplot(dfstat, aes(x = mean.male, y = mean.female, color = celltype)) +
  geom_point() + geom_abline(intercept = 0, slope = 1, col = "red") + 
  theme_bw() + xlab("Male (mean prop.)") + ylab("Female (mean prop.)")
ggpt <- ggpt + facet_wrap(~blood.subgroup)

# save new plot
pt.fname <- "ggpt-groups_mean-prop_sexdiff-blood-cell.pdf"
pdf(pt.fname, width = 5, height = 3); ggpt; dev.off()

#----------------------
# granulocytes by study
#----------------------
# get plot data
pdf <- as.data.frame(pd[pd$blood.subgroup=="PBMC",])
gse.filt <- as.data.frame(table(pdf$gse))
gse.filt <- gse.filt[gse.filt[,2] > 10,]
pdf <- pdf[pdf$gse %in% gse.filt[,1],]
pdf$gran.num <- as.numeric(pdf$predcell.Gran)
# format variables
lvl <- unlist(lapply(gse.filt[,1], function(gsei){median(pdf[pdf$gse==gsei,]$gran.num)}))
names(lvl) <- gse.filt[,1]
pdf$gse.fact <- factor(pdf$gse, levels = names(lvl[order(lvl)]))
# make summary df
pdf.summary <- aggregate(gran.num ~ gse.fact, data = pdf, median)

# get new plot object
ggjitter <- ggplot(pdf, aes(x = gse.fact, y = gran.num)) + 
  geom_boxplot(color = "black", outlier.alpha = 0) +
  geom_jitter(alpha = 0.3, color = "blue") + theme_bw() +
  geom_crossbar(data=pdf.summary, aes(ymin = gran.num, ymax = gran.num),
                size=1, col="red", width = 1) +
  theme(axis.text.x = element_blank()) +
  ylab("Granulocytes\n(proportion)")
ggjitter <- ggjitter + facet_grid(. ~ gse.fact, scales = "free", space = "free")

# save new plot
pdf.fname <- "ggdist_pbmc-gran-bygse2.pdf"
pdf(pdf.fname, width = 5, height = 2.5)
print(ggjitter); dev.off()



