#!/usr/bin/env R

# Author: Sean Maden
#
# Validation of Inoshita et al 2015, sex differences in blood using PBMCs
# compiled from public HM450K and EPIC data.
#
#

# get dependencies for se objects and plots
library(minfi); library(sva); library(gridExtra); 
library(HDF5Array); library(gaston); library(ggrepel)
# get cross-reactive probes from methyPre
devtools::install_github("metamaden/methyPre"); library(methyPre)

#----------
# load data
#----------
# get all metadata for 2 platforms
md.all <- get(load("si1_all-md-2platforms.rda"))
md.blood <- get(load("si2_blood-md-2platforms.rda"))

# get study samples
study.acc <- "GSE67393"; mdf <- md.all[md.all$gse == study.acc,]

#-------------------------------------------------
# validation of male vs female cell type fractions
#-------------------------------------------------
mdv <- md.blood[grepl("^peripheral.*", md.blood$blood_subgroup),]
dim(mdv) # [1] 642  56

# get cell vars, convert to numeric
cellv <- colnames(mdv)[12:17]
for(cvi in cellv){mdv[,cvi] <- as.numeric(mdv[,cvi])}

# get list of outcomes
ltv.mdv <- lapply(cellv, function(x){
  cellv.male <- mdv[mdv$predsex == "M",x]; 
  cellv.female <- mdv[mdv$predsex == "F",x]; 
  t.test(cellv.male, cellv.female)}); names(ltv.mdv) <- cellv

for(cvi in names(ltv.mdv)){
  ltvi <- ltv.mdv[[cvi]]; mean.str <- ltvi$estimate
  names(mean.str) <- c("mean in male", "mean in female")
  message(cvi, " results: ")
  message("pval = ", ltvi$p.value)
  message(names(mean.str)[1], " = ", round(mean.str[1], 3))
  message(names(mean.str)[2], " = ", round(mean.str[2], 3))
}
# predcell.CD8T results: 
#   pval = 0.00122322371014183
# mean in male = 0.22
# mean in female = 0.197
# predcell.CD4T results: 
#   pval = 0.418961160868054
# mean in male = 0.333
# mean in female = 0.341
# predcell.NK results: 
#   pval = 0.000373548795124694
# mean in male = 0.14
# mean in female = 0.106
# predcell.Bcell results: 
#   pval = 0.0166006229564799
# mean in male = 0.146
# mean in female = 0.135
# predcell.Mono results: 
#   pval = 0.0110032055655048
# mean in male = 0.132
# mean in female = 0.151
# predcell.Gran results: 
#   pval = 2.85710169808682e-08
# mean in male = 0.038
# mean in female = 0.085

#-----------------------------
# check hm450k pbmcs sex, ages
#-----------------------------
md <- get(load("si2_blood-md-2platforms.rda"))
md <- md[grepl("^peripheral.*", md$blood_subgroup) & md$platform == "hm450k",]
summary(as.numeric(md$predage))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.606  13.457  28.165  32.168  48.541  89.203 
summary(as.numeric(md[md$predsex == "M",]$predage))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 5.758  10.510  35.733  33.447  53.045  89.203 
summary(as.numeric(md[md$predsex == "F",]$predage))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.606  18.924  27.567  31.605  44.289  80.117
tt <- t.test(as.numeric(md[md$predsex == "M",]$predage),
       as.numeric(md[md$predsex == "F",]$predage))
tt$p.value # 0.4504652

#-------------------------------------
# process hm450k pbmc's for validation
#-------------------------------------
# get the sample ids
md.fpath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_results",
                      "si2_blood-md-2platforms.rda")
md <- get(load(md.fpath)); md <- md[grepl("^peripheral.*", md$blood_subgroup),]
dim(md); table(md$platform)
# [1] 642  56
# epic hm450k
# 256    386

# load gr sets containing noob-norm dnam
gr.fpath.hm450k <- file.path("eternity", "recount-methylation", "recount-methylation-hm450k", "rmi",
                             "remethdb_hm450k_h5se_gr_merged_1619736651-1607018051_0-0-3")
gr.hm450k <- loadHDF5SummarizedExperiment(gr.fpath.hm450k)
gr.hm450k <- gr.hm450k[,gr.hm450k$gsm %in% md$gsm]; dim(gr.hm450k) # [1] 485512    386
gr.fpath.epic <- file.path("eternity", "recount-methylation", "recount-methylation-epic", "rmi",
                           "remethdb_h5se_gr_epic-hm850k_merged_1621537799-1607018051_0-0-3")
gr.epic <- loadHDF5SummarizedExperiment(gr.fpath.epic)
gr.epic <- gr.epic[,gr.epic$gsm %in% md$gsm]; dim(gr.epic) # [1] 866836    256

gr.combined <- combineArrays(gr.epic, gr.hm450k, outType = "IlluminaHumanMethylation450k")
gr.combined.savepath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_results",
                                  "remethdb_h5se_gr_pbmc-hm450k-epic-merged_0-0-3")
saveHDF5SummarizedExperiment(x = gr.combined, dir = gr.combined.savepath)

#-----------------------------
# process PBMCs for validation
#-----------------------------
# load the noob norm/study bias corrected grset
# gr.path <- "gr-adj_hm450k_blood-group-peripheral_blood_mononuclear_cells.rda"
gr.path <- "remethdb_h5se_gr_pbmc-hm450k-epic-merged_0-0-3"
gr <- loadHDF5SummarizedExperiment(gr.path)
dim(gr)

# append pheno data
pheno.fname <- 'si2_blood-md-2platforms.rda'
pheno <- get(load(pheno.fname))
pheno <- pheno[colnames(gr),]
pheno <- pheno[order(match(pheno$gsm, gr$gsm)),]
identical(rownames(pheno), colnames(gr)) # TRUE
colData(gr) <- DataFrame(pheno)

# get filtered probes
man.dpath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_results")
# filter for autosomal probes
anno <- getAnnotation(gr); 
cg.keep <- anno[!anno$chr %in% c("chrY", "chrX"),]$Name
length(cg.keep) # [1] 442474
# filter on detp values
ptable.path <- file.path(man.dpath, "detp-sstats_hm450k-blood-groups.rda")
ptable <- get(load(ptable.path))
varname <- "perc_above_05.peripheral_blood_mononuclear_cells"
cg.keep <- cg.keep[cg.keep %in% ptable[ptable[, varname] == 0, 1]]
length(cg.keep) # 428387
# filter cross-reactive probes
# note: use probes from Chen et al 2013
data(chen_crxcg); cg.keep <- cg.keep[!cg.keep %in% chen.crxcg]
length(cg.keep) # [1] 402470
# get the filtered data
grf <- gr[rownames(gr) %in% cg.keep,]; dim(grf) # [1] 402470    642
table(grf$platform, grf$predsex)
#       F   M
# epic   142 114
# hm450k 268 118
summary(as.numeric(grf$predage))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3.606   9.630  25.838  27.978  41.267  89.203
summary(as.numeric(grf[,grf$platform == "hm450k"]$predage))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3.606  13.457  28.165  32.168  48.541  89.203
summary(as.numeric(grf[,grf$platform == "epic"]$predage))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 4.463   7.629  15.230  21.660  31.922  63.551

# get the SVs
pheno <- colData(grf); cnv <- colnames(pheno)
mval <- na.omit(as.matrix(getM(grf)))
# design model to preserve predsex variable
mod <- model.matrix(~as.factor(predsex), data = pheno) 
# null model
mod0 <- model.matrix(~1, data = pheno) 
# get significant sv's using m-values
sva.results <- sva(mval, mod, mod0)
# Number of significant surrogate variables is:  24
# Iteration (out of 5 ):1  2  3  4  5  >

# sva.results.fname <- "sva-results_pbmc-hm450k_inoshita-2015-validation.rda"
save.dname <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_results")
sva.results.fname <- "sva-results_pbmc-642-2platform_inoshita-2015-validation.rda"
sva.results.fpath <- file.path(save.dname, sva.results.fname)
save(sva.results, file = sva.results.fpath)

#------------------------------------
# do sex t-tests on lm-corrected dnam
#------------------------------------
# get formatted pheno data
pheno <- colData(grf); 
pheno <- pheno[,!grepl("^sv.*", colnames(pheno))]
pheno$predsex <- as.factor(pheno$predsex)
pheno$predage <- as.numeric(pheno$predage)
pheno <- as.data.frame(pheno, stringsAsFactors = F)
# append sva results
svam <- sva.results[[1]]; colnames(svam) <- paste0("sv", seq(ncol(svam)))
pheno <- cbind(pheno, svam)
# append mval matrix to pheno as pheno.all
pheno.all <- cbind(pheno, t(mval)); cnv <- colnames(pheno.all)
pheno.all$predage <- as.numeric(pheno.all$predage)
which.numeric <- which(grepl("^predcell.*", cnv) | grepl("^sv.*", cnv))
for(cname in cnv[which.numeric]){
  pheno.all[,cname] <- as.numeric(pheno.all[,cname]); message(cname)}
pheno.all$platform <- as.factor(pheno.all$platform)
pheno.all$gse <- as.factor(pheno.all$gse)
# parse the design model for lm correction
mod.str <- paste0("model.matrix(~ gse + platform + predage + ",
                  paste0(cnv[grepl("^predcell.*", cnv)], collapse = " + "), 
                  " + ", paste0(cnv[grepl("^sv.*", cnv)], collapse = " + "),
                  ", data = pheno.all)")
mod <- eval(parse(text = mod.str))
cgidv <- rownames(mval); which.male <- pheno.all$predsex == "M"
which.female <- pheno.all$predsex == "F"
# get ttest results
tdf <- do.call(rbind, lapply(cgidv, function(cgid){
  fitted.mval <- lm.fit(mod, pheno.all[,cgid])$fitted.values
  fitted.bval <- minfi::ilogit2(fitted.mval) # convert fitted mvals to bvals
  ttesti <- t.test(fitted.bval[which.male], fitted.bval[which.female])
  tdfi <- as.data.frame(matrix(c(cgid, ttesti$p.value, ttesti$estimate[1],
                                 ttesti$estimate[2]), nrow = 1), 
                        stringsAsFactors = F)
  return(tdfi)}))
colnames(tdf) <- c("cgid", "ttest.pnom", "mean.male", "mean.female")
# save ttest results
tdf.fname <- "ttest-df-results-mvalfit_pbmc-2platforms_inoshita-2015-validate.rda"
tdf.fpath <- file.path(save.dname, tdf.fname); save(tdf, file = tdf.fpath)
message("done")

#-----------------
# analyze results
#-----------------
dmp.study.fname <- "table1_sex-dmp_inoshita-2015.csv"
dmp.study <- read.csv(dmp.study.fname)
tdf.fname <- "ttest-df-results-mvalfit_pbmc-2platforms_inoshita-2015-validate.rda"
tdf.meta <- get(load(tdf.fname))
for(c in c(2:4)){tdf.meta[,c] <- as.numeric(tdf.meta[,c])}
# get bonferroni-adj pvalues
tdf.meta$p.bf <- p.adjust(tdf.meta$ttest.pnom, method = "bonferroni")

# show DMP intersect with max.dmp top significant DMPs
max.dmp <- 300; tdf.meta <- tdf.meta[order(tdf.meta$p.bf),]
length(intersect(dmp.study$Target.ID, tdf.meta[c(1:300),1])) # [1] 48
max.dmp <- 1000; tdff.meta <- tdf.meta[c(1:max.dmp),]
summary(tdff.meta[,3] - tdff.meta[,4])
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.157345 -0.038896 -0.015107 -0.005023  0.030754  0.138592
summary(abs(tdff.meta[,3] - tdff.meta[,4]))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002866 0.024043 0.035663 0.038040 0.048390 0.157345 
length(intersect(dmp.study$Target.ID, tdff.meta$cgid)) # [1] 80
max.dmp <- 10000; tdff.meta <- tdf.meta[c(1:max.dmp),]
summary(tdff.meta[,3] - tdff.meta[,4])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.16016 -0.00428  0.01511  0.01199  0.02874  0.13859
summary(abs(tdff.meta[,3] - tdff.meta[,4]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00104 0.01227 0.02134 0.02510 0.03464 0.16016
length(intersect(dmp.study$Target.ID, tdff.meta$cgid)) # [1] 119

# make 2x2 results summary matrix
dmp.pos1 <- dmp.study$Target.ID; dmp.pos2 <- tdff.meta$cgid
dmp.neg1 <- tdf.meta$cgid[!tdf.meta$cgid %in% dmp.pos1]
dmp.neg2 <- tdf.meta$cgid[!tdf.meta$cgid %in% dmp.pos2]
dmp.pos1.pos2 <- intersect(dmp.pos1, dmp.pos2)
dmp.pos1.neg2 <- intersect(dmp.pos1, dmp.neg2)
dmp.neg1.pos2 <- intersect(dmp.neg1, dmp.pos2)
dmp.neg1.neg2 <- intersect(dmp.neg1, dmp.neg2)
mdmp.comp <- matrix(c(length(dmp.pos1.pos2), length(dmp.pos1.neg2),
                      length(dmp.neg1.pos2), length(dmp.neg1.neg2)), nrow = 2)
colnames(mdmp.comp) <- c("study.pos", "study.neg")
rownames(mdmp.comp) <- c("meta.pos", "meta.neg")
mdmp.comp
# study.pos study.neg
# meta.pos        79       250
# meta.neg       213    430707

#-------------------
# make results plots
#-------------------
library(ggrepel)
# plot top 10k significant probes
tdff.meta <- tdf.meta[c(1:10000),]
# get the max nominal pvalue
pnom <- max(tdff.meta$ttest.pnom)
# define the cg group colors, points alpha
col.neg1.neg2 <- "#A7A7A7"; col.neg1.pos2 <- "#DFD970"; 
col.pos1.neg2 <- "#FF4B4B"; col.pos1.pos2 <- "#339FFF"
alpha <- 0.4 # point transparency value
# assign color labels
tdf.meta$cgid.category <- ifelse(tdf.meta$cgid %in% dmp.pos1.pos2, "pos_pos", 
                                 ifelse(tdf.meta$cgid %in% dmp.neg1.neg2, "neg_neg",
                                        ifelse(tdf.meta$cgid %in% dmp.pos1.neg2, 
                                               "pos_neg", "neg_pos")))
# get volcano plot axes
tdf.meta$neg.l10.pval <- -1*log10(tdf.meta$ttest.pnom)
tdf.meta$sex.diff.minusf <- tdf.meta$mean.male - tdf.meta$mean.female
# get cgid label plot
cond.label <- tdf.meta$cgid.category == "pos_pos" & 
  abs(tdf.meta$sex.diff.minusf) > 0.08 & tdf.meta$neg.l10.pval > 50
dflabel <- tdf.meta[cond.label,]

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
  xlab("Original DMPs\n(Inoshita et al 2015)") + ylab("PBMC DMPs")

# volcano plot, labelling probes as in mdmp.comp
gg.vp <- ggplot() + ylim(0, 140) +
  geom_point(data = tdf.meta[tdf.meta$cgid.category == "neg_neg",], 
             aes(x = sex.diff.minusf, y = neg.l10.pval), 
             color = col.neg1.neg2, alpha = alpha) +
  geom_point(data = tdf.meta[tdf.meta$cgid.category == "neg_pos",], 
             aes(x = sex.diff.minusf, y = neg.l10.pval), 
             color = col.neg1.pos2, alpha = alpha) +
  geom_point(data = tdf.meta[tdf.meta$cgid.category == "pos_neg",], 
             aes(x = sex.diff.minusf, y = neg.l10.pval), 
             color = col.pos1.neg2, alpha = alpha) +
  geom_point(data = tdf.meta[tdf.meta$cgid.category == "pos_pos",], 
             aes(x = sex.diff.minusf, y = neg.l10.pval), 
             color = col.pos1.pos2, alpha = alpha) +
  geom_hline(yintercept = -1*log10(pnom), color = "black") + 
  geom_vline(xintercept = 0, color = "black") + theme_bw() +
  xlab("Mean Beta-value diff.\n(Male - Female, PBMCs)") + ylab("-log10(p-value)") + 
  geom_text_repel(data = dflabel, 
                  aes(label = cgid, x = sex.diff.minusf, y = neg.l10.pval),
                  segment.color = 'grey50', size = 4, force = 100, segment.size = 0.2, 
                  box.padding = 1, point_padding = 0, segment.color = "grey50", 
                  direction = "both")

# scatter plot of sex diffs for overlappig probes
cgidv <- intersect(tdff.meta$cgid, dmp.study$Target.ID)
dmp.studyf <- dmp.study[dmp.study$Target.ID %in% cgidv,]
tdff.metaf <- tdff.meta[tdff.meta$cgid %in% cgidv,]
tdff.metaf <- tdff.metaf[order(match(tdff.metaf$cgid, dmp.studyf$Target.ID)),]
identical(tdff.metaf$cgid, dmp.studyf$Target.ID)
mean.diff.ino <- dmp.studyf$Gender.average.difference.of.β.value.
mean.diff.meta <- tdff.metaf$sex.diff.minusf
dfp.pt <- data.frame(cgid = tdff.metaf$cgid, mdiff.ino = mean.diff.ino, 
                  mdiff.meta = mean.diff.meta, stringsAsFactors = F)
dflabel.pt <- dfp.pt[dfp.pt$cgid %in% dflabel$cgid,]
gg.pt <- ggplot(dfp.pt, aes(x = mdiff.ino, y = mdiff.meta)) +
  geom_point(alpha = 0.5, color = col.pos1.pos2) + 
  geom_abline(intercept = 0, slope = 1, col = "black") +
  geom_hline(yintercept = 0, col = "black") +
  geom_vline(xintercept = 0) + theme_bw() +
  xlim(-0.2, 0.2) + ylim(-0.2, 0.2) +
  ylab("Mean Beta-value diff.\n(Male - Female, PBMCs)") +
  xlab("Mean Beta-value diff.\n(Male - Female, Inoshita et al 2015)") +
  geom_text_repel(data = subset(dfp.pt, cgid %in% dflabel$cgid),
                label = subset(dfp.pt, cgid %in% dflabel$cgid)$cgid,
                size = 4, box.padding = 1.5, point_padding = 0, 
                force = 100, segment.size = 0.4, segment.color = "grey50", 
                direction = "both")

# make the qqplot
# note: borrows code from gaston::qqplot.pvalues()
p <- tdf.meta$ttest.pnom; names(p) <- tdf.meta[,1]
p <- p[!is.na(p)]; w <- (p == 0); p <- p[!w]; args <- list()
args$xlab <- expression(paste("expected ", -log[10](p)))
args$ylab <- expression(paste("observed ", -log[10](p)))
args$main <- "QQ plot of p-values"; n <- length(p); args$type <- "n"
expected <- -log10((n:1)/(n + 1)); observed <- sort(-log10(p))
w <- gaston:::manhattan.thinning(expected, observed, 10000, 10000)
args$x <- expected[w]; args$y <- observed[w]; args$type <- "p"
tdfi <- tdf.meta[order(match(tdf.meta[,1], names(observed))),]
identical(tdfi[,1], names(observed)); tdf.thinned <- tdfi[w,]
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
  geom_point(data = dfp.qq[dfp.qq$color == col.neg1.pos2,], aes(x = expected, y = observed), 
             color = col.neg1.pos2, alpha = 0.4) +
  geom_point(data = dfp.qq[dfp.qq$color == col.pos1.pos2,], aes(x = expected, y = observed), 
             color = col.pos1.pos2, alpha = 0.4) + theme_bw() + 
  xlab("Expected\n(-log10[p-value])") + ylab("Observed \n(-log10[p-value])") +
  geom_text_repel(data = subset(dfp.qq, cgid %in% dflabel$cgid),
                  x = subset(dfp.qq, cgid %in% dflabel$cgid)$expected,
                  y = subset(dfp.qq, cgid %in% dflabel$cgid)$observed,
                  label = subset(dfp, cgid %in% dflabel$cgid)$cgid,
                  size = 4, box.padding = 1.8, point.padding = 0, 
                  force = 100, segment.size  = 0.4, segment.color = "grey50", 
                  direction = "both")

# get concordance at top plot
tdf.meta <- tdf.meta[order(tdf.meta$p.bf),]
num.dmp <- c(seq(0, 950, 50), seq(1000, nrow(tdf.meta), 1000))
num.validated <- sapply(num.dmp, function(x){
  length(intersect(tdf.meta[c(1:x),1], dmp.study[,2]))})
dfp.cat <- data.frame(num.dmp = num.dmp, 
                      num.validated = num.validated,
                      stringsAsFactors = F)
dfp.cat$color <- ifelse(dfp.cat$num.validated <= 119, 
                       col.pos1.pos2, "black")
# full plot
max.yval <- max(dfp.cat$num.validated)
max.xval <- max(dfp.cat$num.dmp)
sig.xline <- 1e4
sig.yline <- dfp.cat[dfp.cat$num.dmp == sig.xline, 2]
gg.cat1 <- ggplot(dfp.cat, aes(x = num.dmp, y = num.validated, colour = color)) +
  geom_point() + geom_line() + theme_bw() + 
  scale_color_manual(name = "color", values = c(col.pos1.pos2, "black")) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) + 
  ylab("Validated DMPs") + xlab("Total DMPs") +
  annotate("text", label = as.character(max.yval), 
           x = 3500, y = max.yval - 10) +
  geom_hline(yintercept = max.yval) +
  annotate("text", x = max.xval - 55000, y = -10, label = format(max.xval, 
                          scientific = T, digits = 3)) +
  geom_vline(xintercept = max.xval) +
  geom_vline(xintercept = sig.xline, linetype = 2) +
  geom_hline(yintercept = sig.yline, linetype = 2) +
  annotate("text", x = sig.xline + 55000, y = sig.yline - 10,
           label = paste0("(",format(sig.xline, scientific = T, digits = 3), 
                          ", ", sig.yline,")"))
# get zoom plot and add grob to annotation_custom()
gg.cat2 <- ggplot(dfp.cat[dfp.cat$num.dmp <= 1e4,], 
                  aes(x = num.dmp, y = num.validated)) +
  geom_point(color = col.pos1.pos2) + geom_line(color = col.pos1.pos2) + theme_bw() +
  ylim(0, 130) + xlim(0, 1e4) + xlab("") + ylab("") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.x = element_text(angle = 90))
grobfun <- function(plot){annotation_custom(annotation_custom(ggplotGrob(plot)))}
gg.grob <- dplyr::mutate(grobs = purrr::pmap(gg.cat2, grobfun))
# full plot
gg.cat <- gg.cat1 + 
  annotation_custom(ggplotGrob(gg.cat2), 
                    xmin=1.3e5, xmax=3.8e5, 
                    ymin=5, ymax=195)

#------------
# store plots
#------------
# store individual plots
pdf("fig2a-cat.pdf", 4.5, 4); gg.cat; dev.off()
pdf("fig2b-grid.pdf", 3, 3); ggtile; dev.off()
pdf("fig2c-vp.pdf", 4, 5); gg.vp; dev.off()
pdf("fig2d-qq.pdf", 4, 4); gg.qq; dev.off()
pdf("fig2e-pt.pdf", 4, 4); gg.pt; dev.off()

# store the composite plot
val1.comp.fname <- "ggcomp-tile-vp-mvalfilt_pbmc-hm450k_val1-inoshita.pdf"
pdf(val1.comp.fname, 7, 6); grid.arrange(ggtile, gg.vp, gg.qq, gg.pt, nrow = 2)
dev.off()

#-----------------------------
# summary stats for manuscript
#-----------------------------
# fraction prior dmps that were validated
num.pos1.pos2 <- nrow(tdf.meta[tdf.meta$cgid.category == "pos_pos",]) # 119
num.pos1.pos2/nrow(dmp.study) # 0.4075342

# fraction of pos1:pos2 dmps with same direction
df1 <- tdf.meta[tdf.meta$cgid.category == "pos_pos",]
df2 <- dmp.study[dmp.study$Target.ID %in% df1$cgid,]
df1 <- df1[order(match(df1$cgid, df2$Target.ID)),]
identical(df1$cgid,df2$Target.ID) # TRUE
dir1 <- df1$sex.diff.minusf
dir2 <- df2$Gender.average.difference.of.β.value.
cond1 <- dir1 > 0 & dir2 > 0
cond2 <- dir1 < 0 & dir2 < 0
length(dir1[which(cond1)]) # 30
length(dir1[which(cond2)]) # 89
length(dir1[which(cond1 | cond2)]) # 119
cor.test(dir1, dir2, method = "spearman") 
# rho = 0.947267, pval < 2.2e-16

# fract validated dmps where mean female > male
length(dir1[which(cond2)])/nrow(df1) # 0.7478992

# fract overall dmps where mean female < male
df3 <- tdf.meta[tdf.meta$cgid.category %in% c("pos_pos", "neg_pos"),]
length(which(df3$sex.diff.minusf < 0))/nrow(df3) # 0.3005

# features at validated probes
library(minfiData)
data(MsetEx); anno <- getAnnotation(MsetEx)
cnv <- c("chr", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group")
annof <- anno[anno$Name %in% df1$cgid, cnv]
table(annof$Relation_to_Island)
# Island N_Shore OpenSea S_Shelf S_Shore 
# 65      11      28       4      11
table(grepl(".*Body.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 86    33
table(grepl(".*TSS.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 86    33
table(annof$UCSC_RefGene_Group == "")
# FALSE  TRUE 
# 76    43
table(grepl(".*5'UTR.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 104    15 
table(grepl(".*3'UTR.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 115     4
table(grepl(".*1stExon.*", annof$UCSC_RefGene_Group))
# FALSE  TRUE 
# 107    12

# features at top validated probes
annof2 <- annof[dflabel$cgid,]
unique(unlist(strsplit(annof2$UCSC_RefGene_Name, ";")))
# [1] "LOC644649" "RFTN1"     "SLC6A4"
# DataFrame with 8 rows and 4 columns
# chr Relation_to_Island UCSC_RefGene_Name UCSC_RefGene_Group
# <character>        <character>       <character>        <character>
#   cg04946709       chr16             Island         LOC644649               Body
# cg17238319        chr3            OpenSea             RFTN1               Body
# cg12691488        chr1             Island                                     
# cg20299935       chr17            OpenSea                                     
# cg21148594       chr14            OpenSea                                     
# cg19544707       chr14             Island                                     
# cg11092486        chr6            S_Shore                                     
# cg05951817       chr17            N_Shore            SLC6A4              5'UTR


