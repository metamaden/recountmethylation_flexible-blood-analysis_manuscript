#!/usr/bin/env R

# Author: Sean Maden
#
# Validation of meta-analysis of differential methylation in 
# cord blood (Merid et al 2020).
#
# Note: no GEO records provided, so assume none of the hm450k 
# cord blood group samples were contained.
#
# Methods:
# # 1. filter samples
# filter samples with maternal age >= 42 weeks/294 days
# define term children as having gest age > 37 weeks/259 days & < 42 weeks
# exclude multiple births
# filter out sex-chr probes
#
# 2. do diff DNAm analysis in cord blood
# multiple regression adjustment vars -- 
# * sex, 
# * maternal age (yrs), 
# * maternal social class, 
# * maternal smoking status (3-state), 
# * parity (0 or >1 prior children),
# * birth weight (grams),
# * batch or SVs
# * blood cell type est from minfi
#
# 3. interpret results
# Filter on I^2 heterogeneity
# Bonferroni adj p-values
# Select regions with at least 3 adjacent DMPs
# note: p-adj < 0.01, at least 2 probes/region, probes within 1kbp apart
# comp-p and DMRcate results

#------------------------------------
# potential confounding vars by study
#------------------------------------
# GSE66459
# includes some pre-term (<37 weeks) neonates
# all neonates: 
# * born after normotensive, 
# * singleton pregnancies, 
# * presented head first, 
# * vaginal delivery
# no difference in clinical terms across neonates:
# * parity, 
# * sex,
# * gestational age,

# GSE69636
# low-income population
# mothers from Mexico City
# study tracks lead exposure

# GSE132181
# longitudinal samples, cord at birth and PBMC at age 7yr, for 193 children
# matched genotyping array data for 196 children
# deliveries < 34 weeks excluded
# at least 1 parent had history of asthma, allergic rhinitis, or eczema
# includes paired cord blood and PBMCs from 193 children
# maternal cortinine levels used to categorize mothers as smokers/non-smokers

# GSE141065
# genotypes determined by array 
# blood samples all buffy coat samples
# 14 technical replicates included
# obesity study, all mothers had BMI >= 30 kg/m^2
# all singleton pregnancies
# mothers with underlying health conditions excluded
# mothers randomized to an intervention, either:
# * low GI diet
# * reduced saturated fat intake
# * increased physical activity)
# non-intervention mothers received standard antenatal care
#
#---------------------------
# Data availability by study
#---------------------------
# hm450k
# GSE66459 # soft contains gestational age; table 2, saved as csv, match on gsm_title
# GSE69636 # soft file contains gestational age
# GSE85042 # gest age not in soft data; all samples at term, not pre-term or post-term
# GSE104376 # gest age not in soft data; all children newborn # gest age distributions give as: 'Gestational age (days (SD)): (low_anxiety) 277.3 (7.2); (high_anxiety) 277.1 (10.2)'
# GSE104778 # gest age not in soft data; neonates only; gest age given as 41.76 ± 9.48 weeks
# GSE151042 # gest age not in soft data; neonates only; gest age reported as Gestational age at delivery (mean±SD): (male) 39.17 (1.58); (female) 39.22 (1.49)
# GSE152380 # gest age not in soft data; neonates only; gest age at birth reported as: age: (preterm) 33.27 (6.49); (full term) 34.77 (5.03)
#
# epic
# GSE132181 # gestational age in soft data; race in soft data; sex in soft data
# GSE141065 # gestational age not in soft data; sex in soft data; obesity study

library(HDF5Array); library(methyPre)
library(ggrepel); library(ggplot2)

#----------
# load data
#----------
# get the samples metadata
si2 <- get(load("si2_blood-md-2platforms.rda"))
which.cord.hm450k <- si2$blood_subgroup == "cord_blood"
mdv <- si2[which.cord.hm450k,]

# load the study cg info from analyses
cginfo.study.fname <- "table1_validation2-cordblood-meta_merid-et-al-2020.csv"
cginfo.study <- read.csv(cginfo.study.fname)

# load epic clock study info
lf <- list.files(); lff <- lf[grepl("haftorn", lf)]
lht <- lapply(lff, function(fn){read.csv(fn)})

#-----------------------------------
# assign gestational ages, term info
#-----------------------------------
# read in gestational ages from GSE66459, GSE69636, GSE132181
gsev <- c("GSE66459", "GSE69636", "GSE132181")
gsm.id.line.str <- "'!Sample_geo_accession'"
gestage.str <- ".*(gestational_age|gestational age|gestationalage_birth).*"
which.soft.exp <- grepl(paste0(gsev, collapse = "|"), lf) & 
                          grepl(".*\\.soft$", lf)
lf <- list.files(); lff <- lf[which.soft.exp]; lga <- list()
for(fn in lff){
  softi.lines <- readLines(fn)
  which.gsm.start <- which(grepl("\\^SAMPLE.*", softi.lines))
  which.gsm.end <- which(grepl("\\!Sample_data_row_count", softi.lines))
  mga <- do.call(rbind, lapply(seq(length(which.gsm.start)), function(ii){
    indexv <- which.gsm.start[ii]:which.gsm.end[ii]
    gsm.lines <- softi.lines[indexv]
    gsm.title <- gsm.lines[grepl("\\!Sample_title.*", gsm.lines)]
    which.ga <- which(grepl(gestage.str, gsm.lines))
    gsm.gestage <- ifelse(length(which.ga)==1, gsm.lines[which.ga], "NA")
    matrix(c(gsm.lines[1], gsm.title, gsm.gestage), nrow = 1)}))
  for(c in seq(ncol(mga))){mga[,c] <- gsub(".* \\= ", "", mga[,c])}
  colnames(mga) <- c("gsm", "gsm_title", "gestational_age");lga[[fn]] <- mga}

lga.fname <- "l-gestage-soft_2platforms.rda"
save(lga, file = lga.fname)

#-----------------------------
# get pheno data with gest age
#-----------------------------
# load data
lga.fname <- "l-gestage-soft_2platforms.rda"; lga <- get(load(lga.fname))
mdf <- mdv[mdv$gse %in% gsub("\\_.*", "", names(lga)),]

# harmonize gest age
dfga <- as.data.frame(do.call(rbind, lga), stringsAsFactors = F)
dfga$ga <- as.numeric(gsub(".* ", "", dfga[,3]))
which.ga <- dfga[,4] > 100 & !is.na(dfga[,4])
dfga[which.ga,]$ga <- round(dfga[which.ga,]$ga/7, 0)
dfga <- dfga[!is.na(dfga[,4]),]; dim(dfga) # [1] 505   4
# append to mdf
mdf <- mdf[mdf$gsm %in% dfga[,1],]
dfga <- dfga[dfga[,1] %in% mdf$gsm,]
dfga <- dfga[order(match(dfga[,1], mdf$gsm)),]
cond <- identical(dfga[,1], mdf$gsm)
if(!cond){stop("Couldn't match GSM IDs in mdf and dfga.")}
mdf$ga <- dfga$ga
table(mdf$platform, mdf$predsex)
#         F  M
# epic   64 73
# hm450k 40 41
summary(mdf[mdf$platform == "hm450k",]$ga)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 26.00   37.00   39.00   37.83   40.00   42.00
summary(mdf[mdf$platform == "epic",]$ga)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 34.00   38.00   39.00   39.08   40.00   42.00
# identify term vs not term
term.thresh <- 37 # weeks thresh, term vs. not-term pregnancies
table(mdf$gse, mdf$ga <= term.thresh)
# FALSE TRUE
# GSE132181   121   16
# GSE66459     11   11
# GSE69636     49   10

paste0(paste0("'",mdf$gsm, "'"), collapse = ",")
paste0(paste0("'",mdf$ga, "'"), collapse = ",")

# load metadata
save.dpath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_results")
pheno.fname <- 'si2_blood-md-2platforms.rda'
pheno <- get(load(file.path(save.dpath, pheno.fname)))
# get ga data
gsmv <- c('GSM1622768','GSM1622769','GSM1622770','GSM1622771','GSM1622772','GSM1622773','GSM1622774','GSM1622775','GSM1622776','GSM1622777','GSM1622778','GSM1622779','GSM1622780','GSM1622781','GSM1622782','GSM1622783','GSM1622784','GSM1622785','GSM1622786','GSM1622787','GSM1622788','GSM1622789','GSM1704996','GSM1704997','GSM1704999','GSM1705004','GSM1705005','GSM1705006','GSM1705007','GSM1705008','GSM1705009','GSM1705012','GSM1705013','GSM1705014','GSM1705015','GSM1705017','GSM1705018','GSM1705020','GSM1705021','GSM1705022','GSM1705023','GSM1705024','GSM1705026','GSM1705027','GSM1705028','GSM1705029','GSM1705030','GSM1705031','GSM1705036','GSM1705037','GSM1705038','GSM1705039','GSM1705042','GSM1705044','GSM1705045','GSM1705046','GSM1705047','GSM1705048','GSM1705050','GSM1705052','GSM1705053','GSM1705054','GSM1705055','GSM1705056','GSM1705057','GSM1705058','GSM1705059','GSM1705062','GSM1705063','GSM1705064','GSM1705065','GSM1705066','GSM1705067','GSM1705072','GSM1705073','GSM1705074','GSM1705075','GSM1705080','GSM1705081','GSM1705082','GSM1705083','GSM3852241','GSM3852245','GSM3852247','GSM3852249','GSM3852251','GSM3852255','GSM3852257','GSM3852259','GSM3852263','GSM3852265','GSM3852271','GSM3852273','GSM3852275','GSM3852277','GSM3852279','GSM3852281','GSM3852283','GSM3852285','GSM3852287','GSM3852289','GSM3852291','GSM3852293','GSM3852297','GSM3852299','GSM3852301','GSM3852303','GSM3852307','GSM3852311','GSM3852315','GSM3852319','GSM3852321','GSM3852323','GSM3852325','GSM3852329','GSM3852341','GSM3852343','GSM3852349','GSM3852353','GSM3852355','GSM3852357','GSM3852361','GSM3852363','GSM3852365','GSM3852369','GSM3852371','GSM3852373','GSM3852375','GSM3852377','GSM3852381','GSM3852383','GSM3852389','GSM3852391','GSM3852395','GSM3852397','GSM3852401','GSM3852403','GSM3852405','GSM3852409','GSM3852411','GSM3852413','GSM3852415','GSM3852417','GSM3852419','GSM3852421','GSM3852423','GSM3852425','GSM3852427','GSM3852429','GSM3852431','GSM3852433','GSM3852435','GSM3852437','GSM3852439','GSM3852443','GSM3852445','GSM3852447','GSM3852449','GSM3852461','GSM3852463','GSM3852467','GSM3852469','GSM3852471','GSM3852473','GSM3852481','GSM3852483','GSM3852485','GSM3852487','GSM3852489','GSM3852491','GSM3852493','GSM3852495','GSM3852499','GSM3852501','GSM3852503','GSM3852505','GSM3852507','GSM3852509','GSM3852511','GSM3852513','GSM3852515','GSM3852519','GSM3852525','GSM3852527','GSM3852533','GSM3852535','GSM3852537','GSM3852539','GSM3852541','GSM3852543','GSM3852545','GSM3852547','GSM3852549','GSM3852551','GSM3852555','GSM3852559','GSM3852563','GSM3852565','GSM3852569','GSM3852571','GSM3852573','GSM3852577','GSM3852579','GSM3852583','GSM3852585','GSM3852587','GSM3852589','GSM3852591','GSM3852593','GSM3852597','GSM3852603','GSM3852607','GSM3852609','GSM3852619','GSM3852625','GSM3852627','GSM3852629','GSM3852631')
gav <- c('26','27','28','29','30','30','30','31','36','36','37','38','38','39','39','39','39','40','40','40','41','42','36','40','40','39','39','37','39','39','40','39','39','39','40','39','40','40','38','39','40','39','40','37','37','40','40','39','37','40','41','40','36','40','36','40','39','39','38','39','39','41','38','37','38','38','37','39','38','38','38','38','38','38','39','41','40','40','39','39','37','39','39','40','41','40','40','40','41','42','39','39','38','41','40','35','40','40','40','41','41','41','40','41','39','39','39','40','39','39','39','39','35','39','41','40','39','38','41','41','38','40','39','41','39','39','37','39','38','40','39','38','40','36','38','40','36','37','40','37','40','39','40','39','40','38','38','39','38','41','38','38','40','35','38','38','40','38','34','40','40','39','40','39','41','40','38','39.3','39.2','38.1','39','39','40.3','40.3','41','40','40','40','39','39','39.6','37','40','37','38','38','39','40','39','39','39','39','40','39','39','38','41','39','39','37','38','37','39','39','39','40','38','40','41','39','40','37','39','37','37','40','39','40')
dfga <- data.frame(gsm = gsmv, ga = gav, stringsAsFactors = F)
pheno <- pheno[pheno$gsm %in% dfga$gsm,]
dfga <- dfga[order(match(dfga$gsm, pheno$gsm)),]
cond <- identical(pheno$gsm, dfga$gsm)
if(!cond){stop("Couldn't match GSM IDs in dfga and pheno")}
pheno$ga <- dfga$ga

# save pheno
pheno.fname <- "df-pheno-ga_cord-blood-val.rda"
pheno.fpath <- file.path(save.dpath, pheno.fname)
save(pheno, file = pheno.fpath)
# save(pheno, file = pheno.fname)

#--------------
# get dnam data
#--------------
# load gr sets containing noob-norm dnam
gr.hm450k.fpath <- file.path("eternity", "recount-methylation", "recount-methylation-hm450k", "rmi",
                             "remethdb_hm450k_h5se_gr_merged_1619736651-1607018051_0-0-3")
gr.epic.fpath <- file.path("eternity", "recount-methylation", "recount-methylation-epic", "rmi",
                           "remethdb_h5se_gr_epic-hm850k_merged_1621537799-1607018051_0-0-3")
gr.hm450k <- loadHDF5SummarizedExperiment(gr.hm450k.fpath)
gr.epic <- loadHDF5SummarizedExperiment(gr.epic.fpath)
gr.hm450k <- gr.hm450k[,gr.hm450k$gsm %in% pheno$gsm]
gr.epic <- gr.epic[,gr.epic$gsm %in% pheno$gsm]
dim(gr.hm450k) # [1] 485512     81
dim(gr.epic) # [1] 866836    137
# combine gr sets
gr.combined <- combineArrays(gr.epic, gr.hm450k, outType = "IlluminaHumanMethylation450k")
# save the combined gr set
gr.combined.savepath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_results",
                                  "remethdb_h5se_gr_cord-blood-hm450k-epic-merged_0-0-3")
saveHDF5SummarizedExperiment(x = gr.combined, dir = gr.combined.savepath)

# filtered probes
anno <- getAnnotation(gr.combined)
cg.keep <- rownames(anno[!anno$chr %in% c("chrY", "chrX"),])
length(cg.keep) # 442474
# filter on detp values
ptable.path <- file.path(save.dpath, "detp-sstats_hm450k-blood-groups.rda")
ptable <- get(load(ptable.path)); varname <- "perc_above_05.cord_blood"
cg.keep <- cg.keep[cg.keep %in% ptable[ptable[, varname] == 0, 1]]
length(cg.keep) # 347276
# filter cross-reactive probes from Chen et al 2013
data(chen_crxcg); cg.keep <- cg.keep[!cg.keep %in% chen.crxcg]
length(cg.keep) # 324407
# remove filtered probes
grf <- gr.combined[rownames(gr.combined) %in% cg.keep,]
dim(grf) # 324407    218

#--------------------------------------
# do multiple reg, corr tests, lm tests
#--------------------------------------
# get SVs, preserving gest. age
pheno <- pheno[order(match(rownames(pheno), colnames(grf))),]
identical(colnames(grf), rownames(pheno)); cnv <- colnames(pheno)
mval <- na.omit(as.matrix(getM(grf))); pheno$ga <- as.numeric(pheno$ga)
mod <- model.matrix(~as.factor(ga), data = pheno) # design model, preserving ga
mod0 <- model.matrix(~1, data = pheno) # null model
sva.results <- sva(mval, mod, mod0) # significant SVs
# Number of significant surrogate variables is:  14
# Iteration (out of 5 ):1  2  3  4  5
sva.fpath <- file.path(save.dpath, "sva-results-ga_cord-blood-218-2platforms.rda")
save(sva.results, file = sva.fpath)
message("Finished!")

# get the formatted pheno data
svam <- sva.results[[1]]; colnames(svam) <- paste0("sv", seq(ncol(svam)))
pheno.all <- cbind(pheno, svam); pheno.all <- cbind(pheno.all, t(mval))
cnv <- colnames(pheno.all)
which.numeric <- which(grepl("^predcell.*", cnv) | grepl("^sv.*", cnv))
for(cname in cnv[which.numeric]){
  pheno.all[,cname] <- as.numeric(pheno.all[,cname]); message(cname)}
pheno.all$platform <- as.factor(pheno.all$platform)
pheno.all$gse <- as.factor(pheno.all$gse)

# get the model matrix for multiple regressions
model.vars.str <- paste0("gse + platform + predsex + ",
                     paste0(cnv[grepl("^predcell.*", cnv)], collapse = " + "), 
                     " + ", paste0(cnv[grepl("^sv.*", cnv)], collapse = " + "))
mod.str <- paste0("model.matrix(~ ",model.vars.str,", data = pheno.all)") # design model
mod <- eval(parse(text = mod.str))
dim(mod) # [1] 218  25

# do spearman corr's
dfct <- do.call(rbind, lapply(rownames(mval), function(cgid){
  mval.cgid <- pheno.all[,cgid] # get the dnam
  mval.cgid.fitted <- lm.fit(mod, mval.cgid)$fitted # get adj dnam
  cti <- cor.test(mval.cgid.fitted, pheno.all$ga, method = "spearman")
  mcori <- matrix(c(cgid, cti$estimate, cti$p.value), nrow = 1)
  return(mcori)})); colnames(dfct) <- c("cgid", "rho", "p.unadj")
dfct <- as.data.frame(dfct, stringsAsFactors = F)
for(c in c(2,3)){dfct[,c] <- as.numeric(dfct[,c])}
dfct$p.bonferroni <- p.adjust(dfct[,3], method = "bonferroni")
dfct.save.fname <- "df-corr-spear_mval-adj_cord-blood-val.rda"
dfct.save.dpath <- file.path(save.dpath, dfct.save.fname)
save(dfct, file = dfct.save.dpath)
message("done")

# Get lm adjusted pval, etc.
dflm <- do.call(rbind, lapply(rownames(mval), function(cgid){
  lmi.str <- paste0("lm(",cgid,"~ ga + ",model.vars.str, ", data = pheno.all)")
  lmi <- eval(parse(text = lmi.str)); coeffv <- summary(lmi)$coefficients
  return(coeffv[2,])})); rownames(dflm) <- rownames(mval)
dflm.fname <- "df-lm-ga_mval-2platform_cord-blood-val.rda"
dflm.fpath <- file.path(save.dpath, dflm.fname)
save(dflm, file = dflm.fpath)
message("done")

#----------------------------------
# supp table -- gestational age dmp 
#----------------------------------
# get cord blood dmps
dflm <- get(load("df-lm-ga_mval-2platform_cord-blood-val.rda"))
dflm <- as.data.frame(dflm); dflm$padj <- p.adjust(dflm[,4], method = "bonferroni")
nrow(dflm[dflm$padj <= 0.05,]) # [1] 1097
dmpdf.cb <- dflm[dflm$padj <= 0.05,]

# get the merid dmps
dmp1.fname <- "dmps1-nocomplications_cord-validation_merid-et-al-2020.csv"
dmp2.fname <- "dmps2-allbirths_cord-validation_merid-et-al-2020.csv"
dmp1.study <- read.csv(dmp1.fname); dmp2.study <- read.csv(dmp2.fname)
dmpv.merid <- unique(c(dmp1.study$CpGID, dmp2.study$CpGID))

# get the haftorn dmps
lf <- list.files(); lff <- lf[grepl("haftorn", lf)]
lht <- lapply(lff, function(fn){read.csv(fn)})
dmpv.haf <- unique(c(lht[[1]][,1], lht[[2]][,1], lht[[3]][,1]))

# get cgdf
cgv.unique <- unique(c(rownames(dmpdf.cb), dmpv.merid, dmpv.haf))
cgdf <- data.frame(cgid = cgv.unique)
cgdf$is.cord.blood.dmp <- cgdf$is.merid.2020.dmp <- 
  cgdf$is.haftorn.2021.dmp <- FALSE
cgdf$cord.blood.ttest.estimate <- cgdf$cord.blood.ttest.std.err <-
  cgdf$cord.blood.ttest.tvalue <- cgdf$cord.blood.ttest.pnominal <-
  cgdf$cord.blood.ttests.padjusted <- 'NA'
for(ri in seq(nrow(cgdf))){
  cgi <- cgdf[ri,1]; message(cgi[1])
  if(cgi %in% rownames(dmpdf.cb)){
    dmpi <- dmpdf.cb[rownames(dmpdf.cb) == cgi,]
    cgdf[ri,]$cord.blood.ttest.estimate <- round(dmpi$Estimate, 3)
    cgdf[ri,]$cord.blood.ttest.std.err <- round(dmpi$`Std. Error`, 3)
    cgdf[ri,]$cord.blood.ttest.tvalue <- round(dmpi$`t value`, 3)
    cgdf[ri,]$cord.blood.ttest.pnominal <- format(dmpi$`Pr(>|t|)`,
                                                  scientific = 3, digits = 3)
    cgdf[ri,]$cord.blood.ttests.padjusted <- format(dmpi$padj,
                                                    scientific = 3, digits = 3)
    cgdf[ri,]$is.cord.blood.dmp <- TRUE}
  is.cgi <- cgdf$cgid==cgi
  if(cgi %in% dmpv.merid){cgdf[is.cgi,]$is.merid.2020.dmp <- TRUE}
  if(cgi %in% dmpv.haf){cgdf[is.cgi,]$is.haftorn.2021.dmp <- TRUE}
}

# append anno
# append probe annotations
library(minfiData); library(minfiDataEPIC)
data("RGsetEx"); data("RGsetEPIC")
# bind anno for unique probes
anno1 <- getAnnotation(RGsetEx); anno1 <- anno1[,c(19, 24, 25, 26)]
anno2 <- getAnnotation(RGsetEPIC); anno2 <- anno2[,c(19, 22, 23, 24)]
anno2 <- anno2[!rownames(anno2) %in% rownames(anno1),]
anno <- rbind(anno1, anno2)
anno <- anno[rownames(anno) %in% cgdf$cgid,]
cgdf <- cgdf[cgdf$cgid %in% rownames(anno),]
anno <- anno[order(match(rownames(anno), cgdf$cgid)),]
identical(rownames(anno), cgdf$cgid) # TRUE
cgdf$relation.to.island <- anno$Relation_to_Island
cgdf$ucsc.gene.names <- anno$UCSC_RefGene_Name
cgdf$ucsc.gene.groups <- anno$UCSC_RefGene_Group
cgdf$ucsc.gene.acc <- anno$UCSC_RefGene_Accession

# save
st.fname <- "st_gest-age-dmp_all-dmp-info"
write.csv(as.matrix(cgdf), file = paste0(st.fname, ".csv"))
write.table(as.matrix(cgdf), file = paste0(st.fname, ".tsv"))
save(cgdf, file = paste0(st.fname, ".rda"))

#-----------------
# compare dmp sets
#-----------------
library(UpSetR)

# load the probe annotations
library(minfiData); library(minfiDataEPIC)
data(RGsetEx); data(RGsetEPIC); 
anno.hm450k <- getAnnotation(RGsetEx)
anno.epic <- getAnnotation(RGsetEPIC)
annof.inter <- anno.hm450k[rownames(anno.hm450k) %in% rownames(anno.epic),]

# make dmp set upset plot -- figure 4a
pdf.fname <- "upset_gest-age-dmps_3-studies.pdf"
# filter significant probes
dflm <- get(load("df-lm-ga_mval-2platform_cord-blood-val.rda"))
dflm <- as.data.frame(dflm); dflm$padj <- p.adjust(dflm[,4], method = "bonferroni")
nrow(dflm[dflm$padj <= 0.05,]) # [1] 1097
dmp.new <- dflm[dflm$padj <= 0.05,]
# halftorn dmps
ht.cgidv <- unique(unlist(lapply(lht, function(lhi){lhi[,1]})))
ht.cgidv <- intersect(ht.cgidv, rownames(annof.inter))
length(ht.cgidv)
length(intersect(ht.cgidv, rownames(dmp.new))) # [1] 58
# merid dmps
# get the validation matrix
dmp1.fname <- "dmps1-nocomplications_cord-validation_merid-et-al-2020.csv"
dmp2.fname <- "dmps2-allbirths_cord-validation_merid-et-al-2020.csv"
dmp1.study <- read.csv(dmp1.fname); dmp2.study <- read.csv(dmp2.fname)
dmpv.combined.study <- unique(c(dmp1.study[,1], dmp2.study[,1]))
# get set data as list
ldmp <- list()
ldmp[["cord_blood DMPs"]] <- rownames(dmp.new)
ldmp[["Clock DMPs (Haftorn et al 2021)"]] <- ht.cgidv
ldmp[["Combined DMPs (Merid et al 2020)"]] <- dmpv.combined.study
# plot results
pdf(pdf.fname, 5.5, 2.7)
upset(fromList(ldmp), order.by = "freq",
      mainbar.y.label = paste0(paste0(rep("\n", 20),collapse =""), 
                               "Intersection size (probes)"),
      sets.x.label = "Total set size (probes)")
dev.off()

# make the probe type barplots -- figure 4b
pdf.fname <- "barplot_islregions-anno_ga-dmps-3sets.pdf"
# load the probe annotations
annof.hm450k <- anno.hm450k[!rownames(anno.hm450k) %in% rownames(anno.epic),c(19,24,26)]
annof.epic <- anno.epic[,c(19,22,24)]
anno.all <- rbind(annof.hm450k, annof.epic)
# get plot data for barplots
isl.typev <- c("Island", "Shore", "Shelf", "OpenSea")
dfp <- do.call(rbind, lapply(isl.typev, function(isl.type){
  # isl.type <- "Island"
  num.isl.bg <- nrow(annof.inter[grepl(isl.type, annof.inter$Relation_to_Island),])
  fract.isl.bg <- num.isl.bg/nrow(annof.inter) # 0.305754
  # new dmps
  annoi.all <- anno.all[rownames(anno.all) %in% ldmp[[1]],]
  num.isl <- nrow(annoi.all[grepl(isl.type, annoi.all$Relation_to_Island),])
  fract.isl.dmp <- num.isl/length(ldmp[[1]])
  # merid et al dmps
  annoi.all <- anno.all[rownames(anno.all) %in% ldmp[[3]],]
  num.isl.merid <- nrow(annoi.all[grepl(isl.type, annoi.all$Relation_to_Island),])
  fract.isl.dmp.merid <- num.isl.merid/length(ldmp[[3]])
  # haftorn et al dmps
  ht.dmp <- intersect(ldmp[[2]], rownames(annof.inter))
  annoi.all <- anno.all[rownames(anno.all) %in% ht.dmp,]
  num.isl.haftorn <- nrow(annoi.all[grepl(isl.type, annoi.all$Relation_to_Island),])
  fract.isl.dmp.haftorn <- num.isl.haftorn/length(ht.dmp)
  # return results
  dfpi <- data.frame(region = rep(isl.type, 4),
                     dmp_set = c("Background", "cord_blood", 
                                 "Merid et al 2020", "Haftorn et al 2021"),
                     fract_region = c(fract.isl.bg, fract.isl.dmp, 
                                      fract.isl.dmp.merid, 
                                      fract.isl.dmp.haftorn),
                     stringsAsFactors = FALSE)
  return(dfpi)}))
# make new region barplots
dfp$`Probe set` <- dfp$dmp_set
ggbp <- ggplot(dfp, aes(x = dmp_set, y = fract_region, fill = `Probe set`)) + 
  geom_bar(stat = "identity") + theme_bw() + theme(axis.text.x = element_blank()) + 
  xlab("Probe set") + ylab("Fraction of probes")
bpfinal <- ggbp + facet_wrap(~region)
pdf(pdf.fname, 5, 3); print(bpfinal); dev.off()

# make dmp gene set upset plot -- figure 4c
# get the plot name
pdf.fname <- "upset_ga-dmp-genes_3-sources.pdf"
# get plot data as list
lgene <- list()
lgene[["cord_blood DMP genes"]] <- unique(unlist(strsplit(anno.all[rownames(anno.all) %in% rownames(dmp.new),]$UCSC_RefGene_Name, ";")))
lgene[["Clock DMP genes (Haftorn et al 2021)"]] <- unique(unlist(strsplit(anno.all[rownames(anno.all) %in% ht.cgidv,]$UCSC_RefGene_Name, ";"))) 
lgene[["Combined DMP genes (Merid et al 2020)"]] <- unique(unlist(strsplit(anno.all[rownames(anno.all) %in% dmpv.combined.study,]$UCSC_RefGene_Name, ";")))   
# save the upset plot
pdf(pdf.fname, 5, 2.7)
upset(fromList(lgene), order.by = "freq",
      mainbar.y.label = paste0(paste0(rep("\n", 25),collapse =""), 
                               "Intersection size (genes)"),
      sets.x.label = "Total set size (genes)")
dev.off()

# check for region dmp enrichment
dmpv <- intersect(ldmp$`cord_blood DMPs`, dmp1.study$CpGID)
dmp1.3adj <- dmp1.study[dmp1.study[,21] == "yes",1]
dmpv.3adj <- intersect(dmpv, dmp1.3adj)
length(dmpv.3adj) # 246
length(dmp1.3adj) # 1276
nrow(dmp1.study) # 8899
length(dmp1.3adj)/nrow(dmp1.study) # 0.1433869
length(dmpv.3adj)/length(ldmp$`cord_blood DMPs`) # 0.2242479
bt <- binom.test(x=length(dmpv.3adj), n=length(ldmp$`cord_blood DMPs`),
           p=length(dmp1.3adj)/nrow(dmp1.study))
format(bt$p.value, scientific=T)

length(dmpv)/nrow(dmp1.study)
length(dmpv)/nrow()

#--------------
# assess probes
#--------------
# filter significant probes
dflm <- get(load("df-lm-ga_mval-2platform_cord-blood-val.rda"))
dflm <- as.data.frame(dflm); dflm$padj <- p.adjust(dflm[,4], method = "bonferroni")
nrow(dflm[dflm$padj <= 0.05,]) # [1] 1097
nrow(dflm[dflm$padj <= 0.01,]) # [1] 818

# get the validation matrix
dmp1.fname <- "dmps1-nocomplications_cord-validation_merid-et-al-2020.csv"
dmp2.fname <- "dmps2-allbirths_cord-validation_merid-et-al-2020.csv"
dmp1.study <- read.csv(dmp1.fname); dmp2.study <- read.csv(dmp2.fname)
length(unique(dmp1.study[,1])) # [1] 8899
length(unique(dmp2.study[,1])) # [1] 17095
length(unique(c(dmp1.study[,1], dmp2.study[,1]))) # [1] 17686
# get region-wise info for set dmp1
853/8899 # 0.09585347
table(dmp1.study$Three.or.more.consecutive.CpG.sites)
# no  yes 
# 7623 1276
1276/8899 # 0.1433869
table(dmp1.study[dmp1.study$CpGID %in% dmp.dn,]$Three.or.more.consecutive.CpG.sites)
# no yes 
# 607 246
246/853 # 0.2883939
bt <- binom.test(c(246, 607), p = 1276/8899)
bt$p.value # 1.399976e-27

dmpv.study <- unique(c(dmp1.study[,1], dmp2.study[,1]))
dmp.dn <- rownames(dflm[dflm$padj <= 0.05,])
length(dmp.dn) # [1] 1097
length(intersect(dmp.dn, dmpv.study)) # [1] 1005
length(intersect(dmp.dn, dmp1.study[,1])) # [1] 853
length(intersect(dmp.dn, dmp2.study[,1])) # [1] 994
dmp.pos1.pos2 <- dmp.dn[dmp.dn %in% dmpv.study]
dmp.pos1.neg2 <- dmpv.study[!dmpv.study %in% dmp.dn]
dmp.neg1.pos2 <- dmp.dn[!dmp.dn %in% dmpv.study]
dmp.neg1.neg2 <- rownames(dflm[!rownames(dflm) %in% dmp.dn & 
                        !rownames(dflm) %in% dmpv.study,])
valm <- matrix(c(length(dmp.pos1.pos2), length(dmp.pos1.neg2), 
         length(dmp.neg1.pos2), length(dmp.neg1.neg2)), nrow = 2)
colnames(valm) <- c("pos1", "neg1"); rownames(valm) <- c("pos2", "neg2")
valm
#       pos1   neg1
# pos2  1005     92
# neg2 16681 310234

#--------------
# results plots
#--------------
# plot top significant probes
dflm$cgid <- rownames(dflm); dflmf <- dflm[dflm$cgid %in% dmp.dn,]
pnom <- max(dflmf[,4]) # get the max nominal pvalue
# define the cg group colors, points alpha
col.neg1.neg2 <- "#A7A7A7"; col.neg1.pos2 <- "#DFD970"; 
col.pos1.neg2 <- "#FF4B4B"; col.pos1.pos2 <- "#339FFF"
alpha <- 0.4 # point transparency value
# assign color labels
dflm$cgid.category <- ifelse(dflm$cgid %in% dmp.pos1.pos2, "pos_pos", 
                                 ifelse(dflm$cgid %in% dmp.neg1.neg2, "neg_neg",
                                        ifelse(dflm$cgid %in% dmp.pos1.neg2, 
                                               "pos_neg", "neg_pos")))
# get volcano plot axes
dflm$neg.l10.pval <- -1*log10(dflm[,4])
# get cgid label plot
num.probes <- 8
dflabel <- dflm[order(dflm$padj, 1-abs(dflm$Estimate)),]
dflabel <- dflabel[c(1:num.probes),]

# get the tile plot
dfp <- data.frame(color = c(col.pos1.pos2, col.pos1.neg2, col.neg1.pos2, col.neg1.neg2), 
                  var1 = c("+", "+", "-", "-"), var2 = c("+", "-", "+", "-"),
                  num_cg = c(length(dmp.pos1.pos2), length(dmp.pos1.neg2), 
                             length(dmp.neg1.pos2), length(dmp.neg1.neg2)),
                  stringsAsFactors = F)
colv <- c(col.pos1.pos2, col.neg1.neg2, col.neg1.pos2, col.pos1.neg2)
ggtile <- ggplot(dfp, aes(x = var1, y = var2, fill = color)) + 
  geom_tile() + scale_fill_manual(values = colv) + geom_text(aes(label=num_cg)) + 
  theme_bw() + theme(plot.margin=unit(rep(1, 4),"cm"), legend.position = "none") +
  xlab("Original DMPs\n(Merid et al 2020)") + ylab("Cord blood DMPs")

# volcano plot
gg.vp <- ggplot() + ylim(0, 40) +
  geom_point(data = dflm[dflm$cgid.category == "neg_neg",], 
             aes(x = Estimate, y = neg.l10.pval), 
             color = col.neg1.neg2, alpha = alpha) +
  geom_point(data = dflm[dflm$cgid.category == "pos_neg",], 
             aes(x = Estimate, y = neg.l10.pval), 
             color = col.pos1.neg2, alpha = alpha) +
  geom_point(data = dflm[dflm$cgid.category == "pos_pos",], 
             aes(x = Estimate, y = neg.l10.pval), 
             color = col.pos1.pos2, alpha = alpha) +
  geom_point(data = dflm[dflm$cgid.category == "neg_pos",], 
             aes(x = Estimate, y = neg.l10.pval), 
             color = col.neg1.pos2, alpha = alpha) +
  geom_hline(yintercept = -1*log10(pnom), color = "black") + 
  geom_vline(xintercept = 0, color = "black") + theme_bw() +
  xlab("Gest. age coeff.\n(cord blood multiple reg.)") + 
  ylab("-log10(p-value)") + 
  geom_text_repel(data = dflabel, aes(label = cgid, x = Estimate, y = neg.l10.pval),
                  segment.color = 'grey50', size = 4, force = 100, segment.size = 0.2, 
                  box.padding = 1, direction = "both", max.overlaps = 100)

# scatter plot of coefficients for overlappig probes
cgidv <- intersect(dmp.dn, dmpv.study)
cnv.ol <- intersect(colnames(dmp1.study), colnames(dmp2.study))
dmp12.study <- rbind(dmp1.study[,cnv.ol], dmp2.study[,cnv.ol])
dmp.studyf <- dmp12.study[dmp12.study[,1] %in% cgidv,]
dmp.studyf <- dmp.studyf[!duplicated(dmp.studyf$CpGID),]
dflmf <- dflm[dflm$cgid %in% cgidv,]
dflmf <- dflmf[order(match(dflmf$cgid, dmp.studyf[,1])),]
identical(dflmf$cgid, dmp.studyf[,1])
coeff.study <- dmp.studyf[,8]; coeff.cord <- dflmf$Estimate
dfp.pt <- data.frame(cgid = dflmf$cgid, coeff.study = coeff.study, 
                     coeff.cord = coeff.cord, stringsAsFactors = F)
dflabel.pt <- dfp.pt[dfp.pt$cgid %in% dflabel$cgid,]
gg.pt <- ggplot(dfp.pt, aes(x = coeff.study, y = coeff.cord)) +
  geom_point(alpha = 0.5, color = col.pos1.pos2) + 
  geom_hline(yintercept = 0, col = "black") +
  geom_vline(xintercept = 0) + theme_bw() +
  xlim(-0.004, 0.003) + ylim(-0.3, 0.2) +
  ylab("Gest. age coeff.\n(cord blood multiple reg.)") +
  xlab("Gest. age coeff.\n(Merid et al 2020 multiple reg.)") +
  geom_text_repel(data = subset(dfp.pt, cgid %in% dflabel$cgid),
                  aes(x = coeff.study, y = coeff.cord, label = cgid),
                  size = 4, box.padding = 1.5, force = 100, segment.size = 0.4, 
                  segment.color = "grey50", direction = "both", max.overlaps = 50)

# make the qqplot
# note: borrows code from gaston::qqplot.pvalues()
p <- dflm$`Pr(>|t|)`; names(p) <- dflm$cgid
p <- p[!is.na(p)]; w <- (p == 0); p <- p[!w]; args <- list()
args$xlab <- expression(paste("expected ", -log[10](p)))
args$ylab <- expression(paste("observed ", -log[10](p)))
args$main <- "QQ plot of p-values"; n <- length(p); args$type <- "n"
expected <- -log10((n:1)/(n + 1)); observed <- sort(-log10(p))
w <- gaston:::manhattan.thinning(expected, observed, 10000, 10000)
# append labeled probes in w
w <- c(w, which(names(p) %in% dflabel$cgid & !names(p) %in% names(p[w])))
args$x <- expected[w]; args$y <- observed[w]; args$type <- "p"
dflmi <- dflm[order(match(dflm$cgid, names(observed))),]
identical(dflmi$cgid, names(observed)); dflmi.thinned <- dflmi[w,]
# get plot data
dfp.qq <- data.frame(cgid = names(args$y), expected = args$x, 
                     observed = args$y, stringsAsFactors = F)
dfp.qq$color <- ifelse(dfp.qq$cgid %in% dmp.neg1.neg2, col.neg1.neg2,
                       ifelse(dfp.qq$cgid %in% dmp.pos1.neg2, col.pos1.neg2,
                              ifelse(dfp.qq$cgid %in% dmp.neg1.pos2, col.neg1.pos2,
                                     ifelse(dfp.qq$cgid %in% dmp.pos1.pos2, col.pos1.pos2, "NA"))))
df.segment <- data.frame(x1 = 0, x2 = -log10(1/n), y1 = 0, y2 = -log10(1/n),
                         stringsAsFactors = F)
# make qqplot
gg.qq <- ggplot() + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
                                 data = df.segment, color = "black") +
  geom_segment(x = 0, y = 0, xend = -log10(1/n), yend = -log10(1/n), colour = "black") +
  geom_point(data = dfp.qq[dfp.qq$color == col.neg1.neg2,], aes(x = expected, y = observed), 
             color = col.neg1.neg2, alpha = 0.4) +
  geom_point(data = dfp.qq[dfp.qq$color == col.pos1.neg2,], aes(x = expected, y = observed), 
             color = col.pos1.neg2, alpha = 0.4) +
  geom_point(data = dfp.qq[dfp.qq$color == col.pos1.pos2,], aes(x = expected, y = observed), 
             color = col.pos1.pos2, alpha = 0.4) + 
  geom_point(data = dfp.qq[dfp.qq$color == col.neg1.pos2,], aes(x = expected, y = observed), 
             color = col.neg1.pos2, alpha = 0.4) + theme_bw() + 
  xlab("Expected\n(-log10[p-value])") + ylab("Observed \n(-log10[p-value])") +
  geom_text_repel(data = subset(dfp.qq, cgid %in% dflabel$cgid),
                  aes(x = expected, y = observed, label = cgid),
                  #x = subset(dfp.qq, cgid %in% dflabel$cgid)$expected,
                  #y = subset(dfp.qq, cgid %in% dflabel$cgid)$observed,
                  #label = subset(dfp.qq, cgid %in% dflabel$cgid)$cgid,
                  size = 4, box.padding = 1.8, point.padding = 0, 
                  force = 200, segment.size  = 0.4, segment.color = "grey50", 
                  direction = "both", max.overlaps = 25)

# get concordance at top plot
dflm <- dflm[order(dflm$padj),]
num.dmp <- c(seq(0, 950, 50), seq(1000, nrow(dflm), 1000))
num.validated <- sapply(num.dmp, function(x){
  length(intersect(dflm[c(1:x),]$cgid, dmp12.study[,1]))})
dfp.cat <- data.frame(num.dmp = num.dmp, num.validated = num.validated,
                      stringsAsFactors = F)
dfp.cat$color <- ifelse(dfp.cat$num.validated <= length(intersect(dmp.dn, dmpv.study)), 
                        col.pos1.pos2, "black")
# full plot
max.yval <- max(dfp.cat$num.validated); max.xval <- max(dfp.cat$num.dmp)
sig.xline <- length(dmp.dn); sig.yline <- max(dfp.cat[dfp.cat$num.dmp <= sig.xline, 2])
gg.cat1 <- ggplot(dfp.cat, aes(x = num.dmp, y = num.validated, colour = color)) +
  geom_point() + geom_line() + theme_bw() + 
  scale_color_manual(name = "color", values = c(col.pos1.pos2, "black")) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
  ylab("Validated DMPs") + xlab("Total DMPs") +
  annotate("text", label = as.character(max.yval), 
           x = 2.5e4, y = max.yval - 450) +
  geom_hline(yintercept = max.yval) +
  annotate("text", x = max.xval-3.2e4, y=-10, 
           label = format(max.xval, scientific = T, digits = 3)) +
  geom_vline(xintercept = max.xval) +
  geom_vline(xintercept = sig.xline, linetype = 2) +
  geom_hline(yintercept = sig.yline, linetype = 2) +
  annotate("text", x = sig.xline+38000, y = sig.yline+400,
           label = paste0("(",format(sig.xline, scientific = T, digits = 3), 
                          ", ", sig.yline,")"))
# get zoom plot and add grob to annotation_custom()
gg.cat2 <- ggplot(dfp.cat[dfp.cat$num.dmp <= length(dmp.dn),], 
                  aes(x = num.dmp, y = num.validated)) +
  geom_point(color = col.pos1.pos2) + geom_line(color = col.pos1.pos2) + theme_bw() +
  ylim(0, 1000) + xlab("") + ylab("") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.x = element_text(angle = 90))
# full plot
gg.cat <- gg.cat1 + 
  annotation_custom(ggplotGrob(gg.cat2), 
                    xmin=1.68e5, xmax=3.2e5, 
                    ymin=920, ymax=7.6e3)

#------------
# store plots
#------------
# store individual plots
pdf("fig3a-cat.pdf", 4.5, 4); gg.cat; dev.off()
pdf("fig3b-grid.pdf", 3, 3); ggtile; dev.off()
pdf("fig3c-vp.pdf", 4, 5); gg.vp; dev.off()
pdf("fig3d-qq.pdf", 4, 4); gg.qq; dev.off()
pdf("fig3e-pt.pdf", 4, 4); gg.pt; dev.off()

# store the composite plot
val1.comp.fname <- "ggcomp-tile-vp-mvalfilt_pbmc-hm450k_val1-inoshita.pdf"
pdf(val1.comp.fname, 7, 6); grid.arrange(ggtile, gg.vp, gg.qq, gg.pt, nrow = 2)
dev.off()

#----------------------
# assess probe features
#----------------------
# get array probe annotations
library(minfiData)
library(minfiDataEPIC)
data(RGsetEx); data(RGsetEPIC)
dim(RGsetEx); dim(RGsetEPIC)
anno <- getAnnotation(RGsetEx)
anno <- anno[anno$Name %in% getAnnotation(RGsetEPIC)$Name,]
dim(anno)

# assess the validated DMPs
annof <- anno[anno$Name %in% dmp.pos1.pos2,]
nrow(annof) # 1005
table(annof$Relation_to_Island)
# Island N_Shelf N_Shore OpenSea S_Shelf S_Shore 
# 151      73     171     415      50     145
nrow(annof[annof$UCSC_RefGene_Group == "",]) # 223
415/1005 # percent opensea = 41%
151/1005 # perc island = 15%
(73+50)/1005 # perc shelf = 12%
(171+145)/1005 # perc shore = 31%


table(grepl(".*Body.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 566   439
table(grepl(".*TSS.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 769   236 
table(grepl(".*5'UTR.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 861   144
table(grepl(".*3'UTR.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 958    47
table(grepl(".*1stExon.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 956    49

table(grepl(".*1stExon.*|.*5'UTR.*|.*TSS.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 639   366

table(grepl(".*Body.*", annof$UCSC_RefGene_Group) & 
        annof$Relation_to_Island == "OpenSea")
# FALSE  TRUE 
# 826   179

table(annof[grepl(".*Body.*", annof$UCSC_RefGene_Group),]$Relation_to_Island)
# Island N_Shelf N_Shore OpenSea S_Shelf S_Shore 
# 70      37      66     179      31      56
(37+66+31+56) # 190
190/1005

# most common genes
genev <- unlist(sapply(annof$UCSC_RefGene_Name, function(x){
  unique(unlist(strsplit(x, ";")))}))
dt.gene <- as.data.frame(table(genev))
dt.gene <- dt.gene[rev(order(dt.gene[,2])),]
nrow(dt.gene) # 604
sum(dt.gene[,2]) # 823
table(annof$UCSC_RefGene_Name=="")
# FALSE  TRUE 
# 782   223
782/1005 # perc gene-mapping = 78%
head(dt.gene)
# genev Freq
# 445 RAP1GAP2    7
# 340    MCF2L    7
# 592  ZC3H12D    6
# 424    PSMB8    6
# 455  RNF144A    5
# 421    PRRT1    5

genei <- "RAP1GAP2"
which.gene <- which(grepl(paste0(".*", genei,".*"), annof$UCSC_RefGene_Name))
annof.gene <- annof[which.gene,]
unique(unlist(strsplit(annof.gene$UCSC_RefGene_Group, ";"))) # [1] "Body"

genei <- "MCF2L"
which.gene <- which(grepl(paste0(".*", genei,".*"), annof$UCSC_RefGene_Name))
annof.gene <- annof[which.gene,]
unique(unlist(strsplit(annof.gene$UCSC_RefGene_Group, ";"))) # "Body"   "TSS200"

genei <- "ZC3H12D"
which.gene <- which(grepl(paste0(".*", genei,".*"), annof$UCSC_RefGene_Name))
annof.gene <- annof[which.gene,]
unique(unlist(strsplit(annof.gene$UCSC_RefGene_Group, ";"))) # "TSS200"  "5'UTR"   "1stExon"

genei <- "PSMB8"
which.gene <- which(grepl(paste0(".*", genei,".*"), annof$UCSC_RefGene_Name))
annof.gene <- annof[which.gene,]
unique(unlist(strsplit(annof.gene$UCSC_RefGene_Group, ";"))) # "Body"    "TSS1500" "3'UTR"

# assess the 8 top validated DMPs
annof <- annof[annof$Name %in% dflabel$cgid,]
annof[,c(19, 26, 24)]
# DataFrame with 8 rows and 3 columns
# Relation_to_Island     UCSC_RefGene_Group UCSC_RefGene_Name
# <character>            <character>       <character>
#   cg04347477             Island            5'UTR;5'UTR       NCOR2;NCOR2
# cg02001279             Island                   Body            ARID3A
# cg18598117             Island                   Body            ARID3A
# cg20334115            N_Shelf                  3'UTR             PYCR2
# cg20033652            OpenSea 5'UTR;5'UTR;TSS1500;..   NIN;NIN;NIN;NIN
# cg08817867            OpenSea                                         
# cg09915396            OpenSea              Body;Body RAP1GAP2;RAP1GAP2
# cg06870470            N_Shelf                   Body             DOCK6
annof[annof$Name == "cg20033652",26] # "5'UTR;5'UTR;TSS1500;5'UTR"

#-------------------------------------------------
# do lookup analysis in infant/child blood samples
#-------------------------------------------------
# get gestational ages for pbmcs, 7yo longitudinal
md.ga <- get(load(file.path(save.dpath, "df-pheno-ga_cord-blood-val.rda")))

md.ga <- pheno
gsei <- "GSE132181"; mdf.ga <- md.ga[md.ga$gse == gsei,]
mdf.ga$subject_id <- gsub("\\_.*", "", mdf.ga$gsm_title)

md <- get(load(file.path(save.dpath, "si2_blood-md-2platforms.rda")))
mdf <- md[md$gse == gsei,]
which.pb <- which(!mdf$gsm %in% md.ga$gsm & 
                    !grepl(".*cord.*", mdf$tissue))
mdf.pb <- mdf[which.pb,]
mdf.pb$subject_id <- gsub("\\_.*", "", mdf.pb$gsm_title)
mdf.pb$ga <- "NA"
for(sid in unique(mdf.pb$subject_id)){
  if(sid %in% mdf.ga$subject_id){
    which.mdf.ga <- which(mdf.ga$subject_id == sid)
    which.mdf.pb <- which(mdf.pb$subject_id == sid)
    mdf.pb[which.mdf.pb,]$ga <- mdf.ga[which.mdf.ga,]$ga}
}; mdf.pb <- mdf.pb[!is.na(mdf.pb$ga),]

# get the dmps for lookup
dflm.fname <- "df-lm-ga_mval-2platform_cord-blood-val.rda"
dflm <- get(load(file.path(save.dpath, dflm.fname)))
dflm <- as.data.frame(dflm)
dflm$padj <- p.adjust(dflm[,4], method = "bonferroni")
dmp.dn <- rownames(dflm[dflm$padj <= 0.05,])

# Get lm adjusted pval, etc.
model.vars.str <- paste0("gse + platform + predsex + predage + ",
                         paste0(cnv[grepl("^predcell.*", cnv)], collapse = " + "), 
                         " + ", paste0(cnv[grepl("^sv.*", cnv)], collapse = " + "))
mod.str <- paste0("model.matrix(~ ",model.vars.str,", data = pheno.all)") # design model
mod <- eval(parse(text = mod.str))
dim(mod) # [1] 218  25

dflm <- do.call(rbind, lapply(rownames(mval), function(cgid){
  lmi.str <- paste0("lm(",cgid,"~ ga + ",model.vars.str, ", data = pheno.all)")
  lmi <- eval(parse(text = lmi.str)); coeffv <- summary(lmi)$coefficients
  return(coeffv[2,])}))
rownames(dflm) <- rownames(mval)
dflm.fname <- "df-lm-ga_mval-2platform_cord-blood-val.rda"
dflm.fpath <- file.path(save.dpath, dflm.fname)
save(dflm, file = dflm.fpath)
message("done")