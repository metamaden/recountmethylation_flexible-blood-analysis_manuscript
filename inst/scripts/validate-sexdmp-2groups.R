#!/usr/bin/env R

# Author: Sean Maden
#
# Get results for the Inoshita et al 2015 sex DMP validation in
# whole blood and PBMCs, compare, etc.

library(minfi); library(sva); library(gridExtra); 
library(HDF5Array); library(gaston); library(ggrepel)
library(methyPre); library(UpSetR)
library(ggdist)

#----------
# load data
#----------
# get all metadata for 2 platforms
md.all <- get(load("si1_all-md-2platforms.rda"))
md.blood <- get(load("si2_blood-md-2platforms.rda"))
# get study samples
study.acc <- "GSE67393"; mdf <- md.all[md.all$gse == study.acc,]
mdf <- md.blood[md.blood$gse == study.acc,]

# inoshita et al 2015 dmps
dmp.study.fname <- "table1_sex-dmp_inoshita-2015.csv"
dmp.study <- read.csv(dmp.study.fname)

# ttest validation results
tdf.wb.fname <- "ttest-df-results-mvalfit_wholeblood-2platforms_inoshita-2015-validate.rda"
tdf.pbmc.fname <- "ttest-df-results-mvalfit_pbmc-2platforms_inoshita-2015-validate.rda"
tdf.wb <- get(load(tdf.wb.fname))
tdf.pbmc <- get(load(tdf.pbmc.fname))
# format tdf
for(c in c(2:4)){
  tdf.wb[,c] <- as.numeric(tdf.wb[,c])
  tdf.pbmc[,c] <- as.numeric(tdf.pbmc[,c])}

# load search index results
sit.fpath <- "sitest_results_blood-groups-2platforms.csv"
sit <- data.table::fread(sit.fpath, sep = ",", header = T, data.table = F)
sit <- sit[,c(2:ncol(sit))]
colnames(sit) <- c("gsm", "k=1000", "platform")

#-------------------
# supp table -- dmps
#-------------------
# get df for unique cpgs
tdf.pbmc <- tdf.pbmc[order(tdf.pbmc$ttest.pnom),]
tdf.wb <- tdf.wb[order(tdf.wb$ttest.pnom),]
cgv <- unique(c(tdf.pbmc[c(1:1000),1], tdf.wb[c(1:1000),1], dmp.study$NAME))
cgdf <- data.frame(cgid = cgv, stringsAsFactors = F)
# get cpg categories
cgdf$is.whole.blood.dmp <- cgdf$is.pbmc.dmp <- cgdf$is.inoshita.2015.dmp <- FALSE
cgdf[cgdf[,1] %in% tdf.pbmc[c(1:1000),1],]$is.pbmc.dmp <- TRUE
cgdf[cgdf[,1] %in% tdf.wb[c(1:1000),1],]$is.whole.blood.dmp <- TRUE
cgdf[cgdf[,1] %in% dmp.study$NAME,]$is.inoshita.2015.dmp <- TRUE
# get the dnam directions (male - female), and ttest pvalues
cgdf$dnam.direction.whole.blood.dmp <- cgdf$dnam.direction.pbmc.dmp <-
  cgdf$dnam.direction.inoshita.2015.dmp <- cgdf$ttest.pvalue.whole.blood <-
  cgdf$ttest.pvalue.pbmc <- 0
for(ri in seq(nrow(cgdf))){
  message(cgid)
  cgid <- cgdf[ri,1]
  if(cgid %in% tdf.wb[c(1:1000), 1]){
    tdf.wb.filt <- tdf.wb[tdf.wb$cgid == cgid,]
    dnam.dir <- tdf.wb.filt$mean.male - tdf.wb.filt$mean.female
    cgdf[ri,]$dnam.direction.whole.blood.dmp <- round(dnam.dir, 3)
    cgdf[ri,]$ttest.pvalue.whole.blood <- format(tdf.wb.filt[2], 
                                                 scientific = T, digits = 3)}
  if(cgid %in% tdf.pbmc[c(1:1000),1]){
    tdf.pbmc.filt <- tdf.pbmc[tdf.pbmc$cgid == cgid,]
    dnam.dir <- tdf.pbmc.filt$mean.male - tdf.pbmc.filt$mean.female
    cgdf[ri,]$dnam.direction.pbmc.dmp <- round(dnam.dir, 3)
    cgdf[ri,]$ttest.pvalue.pbmc <- format(tdf.pbmc.filt[2], 
                                          scientific = T, digits = 3)}
  if(cgid %in% dmp.study$NAME){
    dmpi <- dmp.study[dmp.study$NAME == cgid,]
    dnam.dir <- as.numeric(dmpi[12])
    cgdf[ri,]$dnam.direction.inoshita.2015.dmp <- round(dnam.dir, 2)}
}

# append probe annotations
library(minfiData); library(minfiDataEPIC)
data("RGsetEx"); data("RGsetEPIC")
# bind anno for unique probes
anno1 <- getAnnotation(RGsetEx); anno1 <- anno1[,c(19, 24, 25, 26)]
anno2 <- getAnnotation(RGsetEPIC); anno2 <- anno2[,c(19, 22, 23, 24)]
anno2 <- anno2[!rownames(anno2) %in% rownames(anno1),]
anno <- rbind(anno1, anno2)
anno <- anno[rownames(anno) %in% cgdf$cgid,]
anno <- anno[order(match(rownames(anno), cgdf$cgid)),]
identical(rownames(anno), cgdf$cgid) # TRUE
cgdf$relation.to.island <- anno$Relation_to_Island
cgdf$ucsc.gene.names <- anno$UCSC_RefGene_Name
cgdf$ucsc.gene.groups <- anno$UCSC_RefGene_Group
cgdf$ucsc.gene.acc <- anno$UCSC_RefGene_Accession

# save
st.fname <- "st_sex-dmp_all-dmp-info"
write.csv(as.matrix(cgdf), file = paste0(st.fname, ".csv"))
write.table(as.matrix(cgdf), file = paste0(st.fname, ".tsv"))
save(cgdf, file = paste0(st.fname, ".rda"))

#---------------
# get si results
#---------------
# subset si results
gsmv <- mdf$gsm
sif <- sit[sit[,1] %in% gsmv,]
dim(sif) # [1] 117   3
# get fract group by pos
md <- md.blood; md$group <- md$blood_subgroup
md[md$group %in% c("all", "NA", "blood_spot", "cord_blood"),]$group <- "other/NOS"
md[grepl("^peripheral.*", md$group),]$group <- "PBMC"
table(md$group)
# other/NOS   PBMC whole_blood 
# 6236         642        5980 

# make new raincloud plots -- fig 3a
pdf.fname <- "ggrain-sexdmp-sitest_2groups-2platforms.pdf"
# make dfp data for plotting
rownames(md) <- md$gsm
mpos <- do.call(rbind, lapply(seq(nrow(sif)), function(ii){
  labv <- unlist(strsplit(sif[ii,2], ";"))
  gsmv <- gsub("\\..*", "", labv)
  md[gsmv,]$group}))
dfp <- do.call(rbind, lapply(unique(md$group), function(groupi){
  num.groupi <- apply(mpos, 1, function(ri){length(ri[ri==groupi])})
  data.frame(subgroup = rep(groupi, length(num.groupi)),
             num.subgroup = num.groupi, stringsAsFactors = F)}))
# raincloud plots
# color palette
colv <- c("PBMC" = "#CD534CFF", "whole_blood" = "#868686FF", 
          "other/NOS" = "#3B3B3BFF")
# order subgroups
dfp$subgroup <- factor(dfp$subgroup, levels = c("whole_blood", "other/NOS","PBMC"))
# make plot
pdf(pdf.fname, 5, 2)
ggplot(dfp, aes(x = subgroup, y = num.subgroup, 
                colour = subgroup, fill = subgroup)) + 
  stat_halfeye(width = 2, justification = -0.16, point_colour = NA, .width = 0) + 
  stat_dots(side = "left", justification = 1.4, binwidth = 2) +
  geom_boxplot(width = 0.4, outlier.colour = NA, alpha = 0.5) + 
  theme_bw() +  guides(fill = "none") + scale_color_manual(values = colv) +
  scale_fill_manual(values = colv) + ylab("Number of samples (k = 1,000 neighbors)") +
  xlab("Subgroup") + coord_flip() + theme(legend.position = "none")
dev.off()

#------------
# upset plots
#------------
library(UpSetR)

# sort tdfs
tdf.wb <- tdf.wb[order(tdf.wb$ttest.pnom),]
tdf.pbmc <- tdf.pbmc[order(tdf.pbmc$ttest.pnom),]

# make the upset plot
pdf.fname <- "upset-sexdmp_validate-2groups.pdf"
lup <- list() # get upset data as list
lup[["Whole blood DMPs"]] <- tdf.wb[c(1:1000),1]
lup[["PBMC DMPs"]] <- tdf.pbmc[c(1:1000),1]
lup[["Inoshita et. al. 2015 DMPs"]] <- dmp.study$NAME
# save the new plot
pdf(pdf.fname, 4.3, 2.6)
upset(fromList(lup), order.by = "freq",
      mainbar.y.label = paste0(paste0(rep("\n", 15),collapse =""), 
                               "Intersection size (probes)"),
      sets.x.label = "Total set size (probes)")
dev.off()

#----------------------------------
# concordance at the top, composite
#----------------------------------
library(ggsci); library(patchwork); library(png)

# plot concodrance at the top data
pdf.fname <- "ggcat-comp_sexdmp-val_2groups-2platforms.pdf"
# get concordance at top data
ltdf <- list("whole_blood" = tdf.wb, "PBMC" = tdf.pbmc)
num.dmp <- c(seq(0, 950, 50), seq(1000, nrow(tdf), 1000))
dfp.cat <- do.call(rbind, lapply(seq(2), function(ii){
  tdf <- ltdf[[ii]]; groupi <- names(ltdf)[ii]
  tdf <- tdf[order(tdf$ttest.pnom),]
  num.validated <- sapply(num.dmp, function(x){
    length(intersect(tdf[c(1:x),1], dmp.study[,2]))})
  dfp <- data.frame(num.dmp = num.dmp, 
                    num.validated = num.validated,
                    stringsAsFactors = F)
  dfp$color <- colv[ii]; dfp$group <- groupi
  return(dfp)}))
dfp.cat$Subgroup <- dfp.cat$group
# composite plot
colv <- c("PBMC" = "#CD534CFF", "whole_blood" = "#868686FF")
ndmp.tot.pbmc <- max(dfp.cat[dfp.cat$Subgroup == 'PBMC',]$num.validated)
ndmp.tot.wb <- max(dfp.cat[dfp.cat$Subgroup == 'whole_blood',]$num.validated)
ymax.zoom <- 100; xmax.zoom <- 1e3; ymin.zoom <- xmin.zoom <- 0

plot.main <- ggplot(dfp.cat, aes(x = num.dmp, y = num.validated, 
                                 group = Subgroup, color = Subgroup)) +
  annotate("rect", xmin = -8000, xmax = 15000, 
           ymin = 0, ymax = 100, 
           alpha = 0.45, fill = "grey") +
  scale_color_manual(values = colv) +
  ylab("Validated DMPs") + xlab("Total DMPs") +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = ndmp.tot.pbmc, color = "#CD534CFF",
             linetype = "longdash") +
  geom_hline(yintercept = ndmp.tot.wb, color = "#868686FF",
             linetype = "longdash") +
  annotate("text", color = "#CD534CFF", x = 1000, 
           y = 285, label= ndmp.tot.pbmc) +
  annotate("text", color = "#868686FF", x = 1000, 
           y = 250, label= ndmp.tot.wb) +
  geom_point() + geom_line() + theme_bw() + ylim(0, 300) +
  annotation_raster(readPNG("magnifying_glass_bgtransparent.png"), 
                    xmin = 8000, xmax = 48000, ymin = 20, ymax = 80) +
  annotation_raster(readPNG("rightarrow_bgtransparent.png"), 
                    xmin = 53000, xmax = 110000, ymin = 20, ymax = 80)
  

plot.inset <- ggplot(dfp.cat, aes(x = num.dmp, y = num.validated, 
                                  group = Subgroup, color = Subgroup)) +
  geom_point() + geom_line() + geom_line() + theme_bw() + 
  scale_color_manual(values = colv) +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ylim(ymin.zoom, ymax.zoom) + xlim(xmin.zoom, xmax.zoom)
  
plot.composite <- plot.main + 
  inset_element(plot.inset, left = 0.3, bottom = 0.02, 
                right = 0.98, top = 0.65)

pdf(pdf.fname, 6, 3); print(plot.composite); dev.off()
  
#-------------------
# compare directions
#-------------------
# inoshita dmps
table(dmp.study$Mean.β..value.of.male > dmp.study$Mean.β.value.of.female)
# FALSE  TRUE 
# 237    55
237/292

# whole blood
tdff.wb <- tdf.wb[c(1:1000),]
table(tdff.wb$mean.male > tdff.wb$mean.female)
# FALSE  TRUE 
# 641   359

# pbmcs
tdff.pbmc <- tdf.pbmc[c(1:1000),]
table(tdff.pbmc$mean.male > tdff.pbmc$mean.female)
# FALSE  TRUE 
# 347   653

# store results as list
ldir <- list()

# check directions
# pbmc
dmp.study$dif <- dmp.study$Mean.β..value.of.male - dmp.study$Mean.β.value.of.female
tdff.wb <- tdf.wb[tdf.wb$cgid %in% dmp.study$NAME,]
dim(tdff.wb) # [1] 263   4
tdff.pbmc <- tdf.pbmc[tdf.pbmc$cgid %in% dmp.study$NAME,]
dim(tdff.pbmc) # [1] 270   4
tdff.pbmc$dif <- tdff.pbmc$mean.male - tdff.pbmc$mean.female
tdff.wb$dif <- tdff.wb$mean.male - tdff.wb$mean.female
table(tdff.wb$dif > 0 & dmp.study$dif > 0)
# FALSE  TRUE 
# 279    13

# corr tests -- whole blood, all DMPs
dmp.int <- intersect(dmp.study$NAME, tdf.wb$cgid)
dmpf <- dmp.study[dmp.study$NAME %in% dmp.int,]
tdff.wb <- tdf.wb[tdf.wb$cgid %in% dmpf$NAME,]
tdff.wb <- tdff.wb[order(match(tdff.wb$cgid, dmpf$NAME)),]
identical(tdff.wb$cgid, dmpf$NAME)
tdff.wb$dif <- as.numeric(tdff.wb$mean.male)-as.numeric(tdff.wb$mean.female)
dmpf$dif <- as.numeric(dmpf$Mean.β..value.of.male) - as.numeric(dmpf$Mean.β.value.of.female)
ldir[["whole_blood"]] <- data.frame(cgid = tdff.wb$cgid, diff.meta = tdff.wb$dif,
                                    diff.inoshita = dmpf$dif, stringsAsFactors = F)
# check direction
table(tdff.wb$dif > 0, dmpf$dif > 0)
#       FALSE TRUE
# FALSE   212    2
# TRUE      1   48
nrow(tdff.wb) # 263
num.agree <- (212+48) # 260
num.agree/nrow(tdff.wb) # 0.9885932
cor.test(tdff.wb$dif, dmpf$dif, method = "spearman")$estimate # 0.9029886

# corr tests -- whole blood, validated DMPs
tdff.wb <- tdf.wb[c(1:1000),]
dmp.int <- intersect(dmp.study$NAME, tdff.wb$cgid)
dmpf <- dmp.study[dmp.study$NAME %in% dmp.int,]
tdff.wb <- tdff.wb[tdff.wb$cgid %in% dmpf$NAME,]
tdff.wb <- tdff.wb[order(match(tdff.wb$cgid, dmpf$NAME)),]
identical(tdff.wb$cgid, dmpf$NAME)
tdff.wb$dif <- as.numeric(tdff.wb$mean.male)-as.numeric(tdff.wb$mean.female)
dmpf$dif <- as.numeric(dmpf$Mean.β..value.of.male) - as.numeric(dmpf$Mean.β.value.of.female)
cor.test(tdff.wb$dif, dmpf$dif, method = "spearman")$estimate # 0.9189878

# corr tests -- PBMC, all DMPs
dmp.int <- intersect(dmp.study$NAME, tdf.pbmc$cgid)
dmpf <- dmp.study[dmp.study$NAME %in% dmp.int,]
tdff.pbmc <- tdf.pbmc[tdf.pbmc$cgid %in% dmpf$NAME,]
tdff.pbmc <- tdff.pbmc[order(match(tdff.pbmc$cgid, dmpf$NAME)),]
identical(tdff.pbmc$cgid, dmpf$NAME)
tdff.pbmc$dif <- as.numeric(tdff.pbmc$mean.male)-as.numeric(tdff.pbmc$mean.female)
dmpf$dif <- as.numeric(dmpf$Mean.β..value.of.male) - as.numeric(dmpf$Mean.β.value.of.female)
ldir[["PBMC"]] <- data.frame(cgid = tdff.pbmc$cgid, diff.meta = tdff.pbmc$dif,
                                    diff.inoshita = dmpf$dif, stringsAsFactors = F)
# check direction
table(tdff.pbmc$dif > 0, dmpf$dif > 0)
#       FALSE TRUE
# FALSE   210    2
# TRUE     10   48
nrow(tdff.pbmc) # 270
num.agree <- (210+48) # 258
num.agree/nrow(tdff.pbmc) # 0.9555556
cor.test(tdff.pbmc$dif, dmpf$dif, method = "spearman")$estimate # 0.8241015

# corr tests -- PBMC, validated DMPs
tdff.pbmc <- tdf.pbmc[c(1:1000),]
dmp.int <- intersect(dmp.study$NAME, tdff.pbmc$cgid)
dmpf <- dmp.study[dmp.study$NAME %in% dmp.int,]
tdff.pbmc <- tdff.pbmc[tdff.pbmc$cgid %in% dmpf$NAME,]
tdff.pbmc <- tdff.pbmc[order(match(tdff.pbmc$cgid, dmpf$NAME)),]
identical(tdff.pbmc$cgid, dmpf$NAME)
tdff.pbmc$dif <- as.numeric(tdff.pbmc$mean.male)-as.numeric(tdff.pbmc$mean.female)
dmpf$dif <- as.numeric(dmpf$Mean.β..value.of.male) - as.numeric(dmpf$Mean.β.value.of.female)
cor.test(tdff.pbmc$dif, dmpf$dif, method = "spearman")$estimate # 0.9386937

#--------------------
# dmp directions plot
#--------------------
pdf.fname <- "scatter_dnam-directions_sex-dmp-2tx.pdf"
# format dfp plot data
for(txi in names(ldir)){ldir[[txi]]$Subgroup <- txi}
dfp <- rbind(ldir[[1]], ldir[[2]])
# get the direction status var
status.cond <- dfp$diff.meta > 0 & dfp$diff.inoshita > 0
status.cond <- status.cond | dfp$diff.meta < 0 & dfp$diff.inoshita < 0
dfp$`Direction\nagreement` <- ifelse(status.cond, "Agree", "Disagree")
# get color palette
colv <- c("PBMC" = "#CD534CFF", "whole_blood" = "#868686FF")
# make the new plot
ggpt_main <- ggplot(dfp, aes(x = diff.inoshita, y = diff.meta, 
                group = Subgroup, color = Subgroup, 
                shape = `Direction\nagreement`)) +
  annotate("rect", xmin = 0, xmax = 0.2, 
           ymin = 0, ymax = 0.15, alpha = 0.2, 
           fill = "goldenrod") +
  annotate("rect", xmin = -0.25, xmax = 0, 
           ymin = -0.15, ymax = 0, alpha = 0.2, 
           fill = "goldenrod") +
  annotate("rect", xmin = -0.25, xmax = 0, 
           ymin = 0, ymax = 0.15, alpha = 0.2, 
           fill = "blue") +
  annotate("rect", xmin = 0, xmax = 0.2, 
           ymin = -0.15, ymax = 0, alpha = 0.2, 
           fill = "blue") +
  geom_point(alpha = 0.5, size = 2) + theme_bw() + 
  scale_color_manual(values = colv) +
  xlab("Mean DNAm difference\n(male - female, Inoshita et al 2015)") +
  ylab("Mean DNAm difference\n(male - female, new analysis)")

# save the new plot
pdf(pdf.fname, 4, 2.8); print(ggpt_main); dev.off()


#-------------------
# analyze probe sets
#-------------------
library(minfiData); library(minfiDataEPIC)
data("RGsetEx"); data("RGsetEPIC")
anno.hm450k <- getAnnotation(RGsetEx)
anno.epic <- getAnnotation(RGsetEPIC)
cgidv.merge <- intersect(rownames(anno.hm450k),rownames(anno.epic))
anno.merge <- anno.hm450k[cgidv.merge,]
# get all probes anno
cnv <- c("Relation_to_Island", "UCSC_RefGene_Group", "UCSC_RefGene_Name")
annof <- anno.hm450k[,cnv]
annof <- rbind(annof, anno.epic[!anno.epic$Name %in% rownames(annof),cnv])
dim(annof)
# get region types
annof$is.promoter <- annof$is.body <- FALSE
which.promoter <- which(grepl(".*TSS.*|.*1stExon.*|.*5'UTR.*", annof[,2]))
which.body <- which(grepl(".*Body.*", annof[,2]))
annof[which.promoter,]$is.promoter <- TRUE
annof[which.body,]$is.body <- TRUE
# get unique gene names
unique.genev <- unlist(lapply(annof$UCSC_RefGene_Name, function(x){
  paste0(unique(unlist(strsplit(x, ";"))), collapse = ";")}))
annof$unique.gene <- unique.genev

# get gene info
geneanno <- function(annof, annobg){
  message("Getting genes among annof probes...")
  tgene <- as.data.frame(table(unlist(strsplit(annof$unique.gene, ";"))))
  tgene <- as.data.frame(tgene[rev(order(tgene[,2])),], stringsAsFactors = F)
  colnames(tgene) <- c("gene.id", "num.dmp")
  tgene$dmp.fract <- tgene[,2]/nrow(annof)
  message("Getting gene background info...")
  tfrat.genev <- do.call(rbind, lapply(seq(nrow(tgene)), function(ii){
    genei <- tgene[ii,1]; ndmpi <- tgene[ii, 2];genei.str <- paste0("(^|;)", genei, "($|;)")
    ncg.annof.genei <- nrow(annobg[grepl(genei.str, annobg$UCSC_RefGene_Name),])
    return(c(ncg.annof.genei, ndmpi/ncg.annof.genei, ncg.annof.genei/nrow(annobg)))
  })); colnames(tfrat.genev) <- c("ncg_gene", "fract.dmp.over.ncg_gene", "fract.bg")
  tgene <- cbind(tgene, tfrat.genev); return(tgene)
}

# get the probe sets
dmp <- dmp.study$NAME
dmp.wb <- tdf.wb$cgid[1:1000]
dmp.pbmc <- tdf.pbmc$cgid[1:1000]
dmp.allol <- intersect(dmp, intersect(dmp.wb, dmp.pbmc))

# get gene anno for datasets
lgene <- list()
lgene[["dmp"]] <- geneanno(annof[dmp,], anno.hm450k)
lgene[["wb"]] <- geneanno(annof[dmp.wb,], anno.merge)
lgene[["pbmc"]] <- geneanno(annof[dmp.pbmc,], anno.merge)
lgene[["all.ol"]] <- geneanno(annof[dmp.allol,], anno.merge)
# summarize overlaps
length(intersect(lgene$dmp[,1], lgene$wb[,1])) # 63
length(intersect(lgene$dmp[,1], lgene$pbmc[,1])) # 58
dim(lgene$all.ol) # [1] 33  6

# anno for inoshita 2015 dmps
annof.dmp <- annof[dmp.study$NAME,]
dim(annof.dmp) # [1] 292   5
tgene.dmp <- as.data.frame(table(unlist(strsplit(annof.dmp$unique.gene, ";"))))
tgene.dmp <- as.data.frame(tgene.dmp[rev(order(tgene.dmp[,2])),], stringsAsFactors = F)
head(tgene.dmp)
# Var1 Freq
# 163    SPEG    5
# 157 SLC4A11    3
# 134  PPFIA3    3
# 131   PEX10    3
# 124   NUPL1    3
# 108     LTK    3
colnames(tgene.dmp) <- c("gene.id", "num.dmp")
tgene.dmp$dmp.fract <- tgene.dmp[,2]/nrow(annof.dmp)
tfrat.genev <- do.call(rbind, lapply(seq(nrow(tgene)), function(ii){
  genei <- tgene.dmp[ii,1]; ndmpi <- tgene.dmp[ii, 2]
  genei.str <- paste0("(^|;)", genei, "($|;)")
  ncg.annof.genei <- nrow(anno.hm450k[grepl(genei.str, anno.hm450k$UCSC_RefGene_Name),])
  return(c(ncg.annof.genei, 
           ndmpi/ncg.annof.genei, 
           ncg.annof.genei/nrow(anno.hm450k)))
})); colnames(tfrat.genev) <- c("ncg_gene", "fract.dmp.over.ncg_gene", "fract.bg")
tgene.dmp <- cbind(tgene.dmp, tfrat.genev)
head(tgene.dmp)

# anno for whole blood
tdf.wb <- tdf.wb[order(tdf.wb$ttest.pnom),]
annof.wb <- annof[tdf.wb$cgid,]
dim(annof.wb) # [1] 292   5
tgene.wb <- as.data.frame(table(unlist(strsplit(annof.wb$unique.gene, ";"))))
tgene.wb <- as.data.frame(tgene[rev(order(tgene[,2])),], stringsAsFactors = F)
head(tgene)

