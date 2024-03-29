#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize Granulocyte proportions in PBMCs by Study/GSE ID.
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
# get pheno data
pd <- pData(gr)

#---------------------------------
# filter samples, format variables
#---------------------------------
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

#----------------------
# granulocytes by study
#----------------------
# make summary df
pdf.summary <- aggregate(gran.num ~ gse.fact, data = pdf, median)

# get new plot object
ggjitter <- ggplot(pdf, aes(x = gse.fact, y = gran.num)) + 
  geom_boxplot(color = "black", outlier.alpha = 0) +
  geom_jitter(alpha = 0.3, color = "blue") + theme_bw() +
  geom_crossbar(data=pdf.summary, aes(ymin = gran.num, ymax = gran.num),
                size=1, col="red", width = 1) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ylab("Granulocytes\n(proportion)")
ggjitter + facet_grid(. ~ gse.fact, scales = "free", space = "free")

ggjitter + facet_grid(~ gse.fact, scales = "free")

# save new plot
pdf.fname <- "ggdist_pbmc-gran-bygse2.pdf"
pdf(pdf.fname, width = 5, height = 2.5)
ggjitter + facet_grid(~ gse.fact, scales = "free")
dev.off()

#------------------------
# cumulative distribution
#------------------------
# plot(ecdf(pdf$gran.num))

ggecdf <- ggplot(pdf, aes(gran.num)) + 
  stat_ecdf(geom = "step", size = 1.2) + theme_bw() +
  geom_vline(xintercept = 0.1, color = "red", size = 1.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Percent of samples") + xlab("Granulocute (proportion)")
  
ggecdf + facet_grid(~gse.fact, scales = "free")

#---------------------------------
# compare qc metrics by gran group
#---------------------------------
# correlation tests
cnv.qc <- colnames(pdf)[c(18:34, 52, 53)]
df.corr <- as.data.frame(do.call(rbind, lapply(cnv.qc, function(cni){
  cti <- cor.test(pdf[,cni], pdf[,"gran.num"], method = "spearman")
  c(cti$estimate, cti$p.value)
})))
colnames(df.corr) <- c("corr.spear.rho", "corr.spear.punadj")
df.corr$p.adj.bh <- p.adjust(df.corr[,2], method = "BH")
df.corr$qc.metric <- cnv.qc
df.corr <- df.corr[,c(4,1:3)]
df.corr <- df.corr[order(df.corr$p.adj.bh),]

# get rounded table
df.corr[,2] <- round(as.numeric(df.corr[,2]), digits = 3)
for(c in c(3:4)){df.corr[,c] <- format(as.numeric(df.corr[,c]), digits = 3)}
# save new table
csv.fname <- "stable_corr_pbmc-gran-qcmetric.csv"
write.csv(df.corr, file = csv.fname, row.names = F)
