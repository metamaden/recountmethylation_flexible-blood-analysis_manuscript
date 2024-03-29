#!/usr/bin/env R

# Author: Sean Maden
#
# Validation of Inoshita et al 2015, sex differences in blood using whole blood
# compiled from public HM450K and EPIC data.
#
#

# get dependencies for se objects and plots
libv <- c("minfi", "sva", "gridExtra", "HDF5Array", "methyPre", "ggrepel", 
          "gaston", "data.table")
sapply(libv, library, character.only = T)

#----------
# load data
#----------
md.fname <- "mdf-epic-hm450k-merge_normal-blood_0-0-3.rda"
md <- get(load(file.path(md.fname)))
md <- md[md$blood.subgroup == "PBMC",]
dim(md) # 627  66

#-----------------------------
# check hm450k pbmcs sex, ages
#-----------------------------
summary(as.numeric(md$predage))
summary(as.numeric(md[md$predsex == "M",]$predage))
summary(as.numeric(md[md$predsex == "F",]$predage))

#----------------------------
# process pbmc for validation
#----------------------------
# get granges
gr.fname <- "gr-noob_h5se_hm450k-epic-merge_0-0-3"
gr <- loadHDF5SummarizedExperiment(gr.fname)
gr <- gr[,colnames(gr) %in% rownames(md)]
dim(gr) # 442474   6866

# filter cross-reactive probes
# note: use probes from Chen et al 2013
library(methyPre)
data(chen_crxcg)
cgfilt <- which(!rownames(gr) %in% chen.crxcg)
gr <- gr[cgfilt,]
dim(gr) # 416045   6866

# filter on granulocytes > 0.25
gran.limit <- 0.25
md$gran.num <- as.numeric(md$predcell.Gran)
qv <- quantile(md$gran.num, seq(0, 1, 1e-2)) # roughly 93rd quantile
which.gran.remove <- md$gran.num >= gran.limit
md <- md[!which.gran.remove,]
gr <- gr[,rownames(md)]
dim(gr) # 416045    580
table(gr$predsex)
# F   M 
# 357 223

#-----------------------
# write the mvals matrix
#-----------------------
mval.fname <- "mval-sexdmp_pbmc_2platforms.csv"
labels.fname <- "samplabels_pbmc_2platforms.csv"
sep.sym <- ","
mval <- t(getM(gr)) # get mvals
samp.labels.fpath <- labels.fname # file.path(save.dpath, labels.fname) # write labels
fwrite(matrix(rownames(mval), ncol = 1), col.names = F, sep = ",", file = samp.labels.fpath)
# write dnam mvals
mval.fpath <- mval.fname # file.path(save.dpath, mval.fname)
cnv <- matrix(colnames(mval), nrow = 1)
data.table::fwrite(cnv, file = mval.fpath, append = F, 
                   col.names = F, row.names = F)
num.samp <- 2000
indexv <- seq(1, nrow(mval), num.samp)
t1 <- Sys.time()
for(ri in indexv){
  start.index <- ri; end.index <- start.index + num.samp - 1
  end.index <- ifelse(end.index > nrow(mval), nrow(mval), end.index)
  # write.csv(dati, file = , append = T, row.names = F, col.names = F)
  data.table::fwrite(as.matrix(mval[start.index:end.index,]), 
                     file = mval.fpath, append = T, 
                     col.names = F, row.names = F)
  message("Finished writing row ", ri, ", time: ", Sys.time()-t1)
}; message("done")

# convert mval table to fh10k
# for fh conversion code, see script: make_fhtable.py

#------------
# get the SVs
#------------
# note: cannot process all of mvals in memory, so 
# this approach uses the 1k hashed features

sva.results.fname <- "sva-results-sexdmp_pbmc-fh10k_2-platform.rda"
# load the fh1k matrix
fh10k.fname <- "mval-fh10k_pbmc.csv"  # "mval-cgfilt-fh10k_pbmc_2platforms.csv"
fh10k.fpath <- fh10k.fname # file.path(save.dpath, fh10k.fname)
fh10k <- data.table::fread(fh10k.fpath, header = F, sep = ",", data.table = F)
rownames(fh10k) <- colnames(gr)
# get pheno data, design matrix
pheno <- colData(gr); cnv <- colnames(pheno)
identical(rownames(pheno), rownames(fh10k)) # TRUE
mod <- model.matrix(~as.factor(predsex), data = pheno)
mod0 <- model.matrix(~1, data = pheno) # null model
# get significant sv's using m-values
sva.results <- sva(t(fh10k), mod, mod0)
# Number of significant surrogate variables is:  17 
# Iteration (out of 5 ):1  2  3  4  5
# save results
save(sva.results, file = sva.results.fname)

#------------------------------------
# do sex t-tests on lm-corrected dnam
#------------------------------------
tdf.fname <- "ttest-df-results-mvalfit_pbmc-2platforms_inoshita-2015-validate.rda"
# get formatted pheno data
pheno <- colData(gr); 
pheno <- pheno[,!grepl("^sv.*", colnames(pheno))]
pheno$predsex <- as.factor(pheno$predsex)
pheno$predage <- as.numeric(pheno$predage)
pheno <- as.data.frame(pheno, stringsAsFactors = F)
# append sva results to pheno
svam <- sva.results[[1]]; colnames(svam) <- paste0("sv", seq(ncol(svam)))
pheno <- cbind(pheno, svam)
# append mval matrix to pheno as pheno.all
num.cg <- 5000; indexv <- seq(1, nrow(gr), num.cg)
# iterate on mvals
# note: impute medians for NAs
t1 = Sys.time(); ltdf <- list()
for(start.index in indexv){
  message("Starting at index ", start.index,
          ", time: ", Sys.time() - t1)
  end.index <- start.index + num.cg - 1
  end.index <- ifelse(end.index > nrow(gr), nrow(gr), end.index)
  # get filtered probe data
  mvali <- getM(gr[start.index:end.index,])
  cgidv <- rownames(mvali)
  # bind all pheno, dnam in pheno.all
  pheno.all <- cbind(pheno, t(mvali))
  cnv <- colnames(pheno.all)
  pheno.all$predage <- as.numeric(pheno.all$predage)
  which.numeric <- which(grepl("^predcell.*", cnv) | grepl("^sv.*", cnv))
  for(cname in cnv[which.numeric]){
    pheno.all[,cname] <- as.numeric(pheno.all[,cname])}
  pheno.all$platform <- as.factor(pheno.all$platform)
  pheno.all$gse <- as.factor(pheno.all$gse)
  # parse the design model for lm correction
  mod.str <- paste0("model.matrix(~ gse + platform + predage + ",
                    paste0(cnv[grepl("^predcell.*", cnv)], collapse = " + "), 
                    " + ", paste0(cnv[grepl("^sv.*", cnv)], collapse = " + "),
                    ", data = pheno.all)")
  mod <- eval(parse(text = mod.str))
  which.male <- pheno.all$predsex == "M"
  which.female <- pheno.all$predsex == "F"
  # get ttest results
  tdf <- do.call(rbind, lapply(cgidv, function(cgid){
    message(cgid)
    cgid.mval <- pheno.all[,cgid]
    cgid.mval[is.na(cgid.mval)] <- median(cgid.mval, na.rm = T)
    fitted.mval <- lm.fit(mod, cgid.mval)$fitted.values
    fitted.bval <- minfi::ilogit2(fitted.mval) # convert fitted mvals to bvals
    ttesti <- t.test(fitted.bval[which.male], fitted.bval[which.female])
    tdfi <- as.data.frame(matrix(c(cgid, ttesti$p.value, ttesti$estimate[1],
                                   ttesti$estimate[2]), nrow = 1), 
                          stringsAsFactors = F)
    return(tdfi)}))
  colnames(tdf) <- c("cgid", "ttest.pnom", "mean.male", "mean.female")
  ltdf[[as.character(start.index)]] <- tdf
  message("Finished index ", start.index, 
          ", time: ", Sys.time() - t1)
}
# save ttest results
tdf <- do.call(rbind, ltdf)
save(tdf, file = tdf.fname)

#----------------
# analyze results
#----------------
dmp.study.fname <- "table1_sex-dmp_inoshita-2015.csv"
dmp.study <- read.csv(dmp.study.fname)

tdf.fname <- "ttest-df-results-mvalfit_pbmc-2platforms_inoshita-2015-validate.rda"
tdf.meta <- get(load(tdf.fname))
for(c in c(2:4)){tdf.meta[,c] <- as.numeric(tdf.meta[,c])}

# get bonferroni-adj pvalues
min.pval <- min(tdf.meta[!tdf.meta$ttest.pnom == 0,]$ttest.pnom) # impute missing pval
which.zero <- which(tdf.meta$ttest.pnom == 0)
if(length(which.zero) > 0){tdf.meta[which.zero,]$ttest.pnom <- min.pval}
tdf.meta$p.bf <- p.adjust(tdf.meta$ttest.pnom, method = "bonferroni")

# show DMP intersect with max.dmp top significant DMPs
max.dmp <- 300; tdf.meta <- tdf.meta[order(tdf.meta$p.bf),]
length(intersect(dmp.study$Target.ID, tdf.meta[c(1:300),1]))
max.dmp <- 1000; tdff.meta <- tdf.meta[c(1:max.dmp),]
summary(tdff.meta[,3] - tdff.meta[,4])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.18097 -0.04147 -0.02474 -0.01203  0.02349  0.13657
summary(abs(tdff.meta[,3] - tdff.meta[,4])) 
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002568 0.024610 0.036712 0.039177 0.050270 0.180967
length(intersect(dmp.study$Target.ID, tdff.meta$cgid))
# [1] 70
max.dmp <- 10000; tdff.meta <- tdf.meta[c(1:max.dmp),]
summary(tdff.meta[,3] - tdff.meta[,4])
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.180967 -0.005567  0.012050  0.008786  0.023422  0.150665 
summary(abs(tdff.meta[,3] - tdff.meta[,4]))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.001361 0.009777 0.018602 0.023071 0.031143 0.180967 
length(intersect(dmp.study$Target.ID, tdff.meta$cgid))
# [1] 119

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
#          study.pos study.neg
# meta.pos       119      9881
# meta.neg       151    405894

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
  abs(tdf.meta$sex.diff.minusf) > 0.045 & tdf.meta$neg.l10.pval > 30
dflabel <- tdf.meta[cond.label,]
dflabel

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
gg.vp <- ggplot() + ylim(0, 320) +
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
  xlab("Mean Beta-value diff.\n(Male - Female, PBMC)") + ylab("-log10(p-value)") + 
  geom_text_repel(data = dflabel, 
                  aes(label = cgid, x = sex.diff.minusf, y = neg.l10.pval),
                  segment.color = 'grey50', size = 3, force = 100, segment.size = 0.2, 
                  box.padding = 1, point_padding = 0, segment.color = "grey50", 
                  direction = "both")

# scatter plot of sex diffs for overlappig probes
cgidv <- intersect(tdff.meta$cgid, dmp.study$Target.ID)
dmp.studyf <- dmp.study[dmp.study$Target.ID %in% cgidv,]
tdff.metaf <- tdff.meta[tdff.meta$cgid %in% cgidv,]
tdff.metaf <- tdff.metaf[order(match(tdff.metaf$cgid, dmp.studyf$Target.ID)),]
identical(tdff.metaf$cgid, dmp.studyf$Target.ID)
mean.diff.ino <- dmp.studyf$Gender.average.difference.of.β.value.
tdff.metaf$sex.diff.minusf <- tdff.metaf$mean.male-tdff.metaf$mean.female
mean.diff.meta <- tdff.metaf$sex.diff.minusf
dfp.pt <- data.frame(cgid = tdff.metaf$cgid, mdiff.ino = mean.diff.ino, 
                     mdiff.meta = mean.diff.meta, stringsAsFactors = F)
dflabel.pt <- dfp.pt[dfp.pt$cgid %in% dflabel$cgid,]
gg.pt <- ggplot(dfp.pt, aes(x = mdiff.ino, y = mdiff.meta)) +
  geom_point(alpha = 0.5, color = col.pos1.pos2) + 
  geom_abline(intercept = 0, slope = 1, col = "black") +
  geom_hline(yintercept = 0, col = "black") +
  geom_vline(xintercept = 0) + theme_bw() +
  xlim(-0.25, 0.25) + ylim(-0.1, 0.1) +
  ylab("Mean Beta-value diff.\n(Male - Female, PBMC)") +
  xlab("Mean Beta-value diff.\n(Male - Female, Inoshita et al 2015)") +
  geom_text_repel(data = subset(dfp.pt, cgid %in% dflabel$cgid),
                  label = subset(dfp.pt, cgid %in% dflabel$cgid)$cgid,
                  size = 3, box.padding = 1.5, point_padding = 0, 
                  force = 100, segment.size = 0.4, segment.color = "grey50", 
                  direction = "both", max.overlaps = 20)

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
  xlim(0, 6) + ylim(0, 350) +
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
  geom_text_repel(data = dfp.qq[dfp.qq$cgid %in% dflabel$cgid,],
                  x = dfp.qq[dfp.qq$cgid %in% dflabel$cgid,]$expected,
                  y = dfp.qq[dfp.qq$cgid %in% dflabel$cgid,]$observed,
                  label = dfp.qq[dfp.qq$cgid %in% dflabel$cgid,]$cgid,
                  size = 3, box.padding = 1.8, point.padding = 0, 
                  force = 100, segment.size  = 0.4, segment.color = "grey50", 
                  direction = "both", max.overlaps = 20)

# get concordance at top plot
tdf.meta <- tdf.meta[order(tdf.meta$p.bf),]
num.dmp <- c(seq(0, 950, 50), seq(1000, nrow(tdf.meta), 1000))
num.validated <- sapply(num.dmp, function(x){
  length(intersect(tdf.meta[c(1:x),1], dmp.study[,2]))})
dfp.cat <- data.frame(num.dmp = num.dmp, 
                      num.validated = num.validated,
                      stringsAsFactors = F)
dfp.cat$color <- ifelse(dfp.cat$num.validated <= 186, 
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
  ylim(0, 200) + xlim(0, 1e4) + xlab("") + ylab("") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.x = element_text(angle = 90))
grobfun <- function(plot){annotation_custom(annotation_custom(ggplotGrob(plot)))}
gg.grob <- dplyr::mutate(grobs = purrr::pmap(gg.cat2, grobfun))
# full plot
gg.cat <- gg.cat1 + 
  annotation_custom(ggplotGrob(gg.cat2), 
                    xmin=1.3e5, xmax=3.5e5, 
                    ymin=5, ymax=150)

#------------
# store plots
#------------
# store individual plots
pdf("figa_pbmc-sexdmp_cat.pdf", width = 4.5, height = 4); gg.cat; dev.off()
pdf("figb_pbmc-sexdmp-grid.pdf", width = 3, height = 3); ggtile; dev.off()
pdf("figc_pbmc-sexdmp-vp.pdf", width = 4, height = 5); gg.vp; dev.off()
pdf("figd_pbmc-sexdmp-qq.pdf", width = 4, height = 4); gg.qq; dev.off()
pdf("fige_pbmc-sexdmp-pt.pdf", width = 4, height = 4); gg.pt; dev.off()