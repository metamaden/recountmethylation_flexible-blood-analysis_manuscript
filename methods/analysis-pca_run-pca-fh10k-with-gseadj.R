#!/usr/bin/env R

# Author: Sean Maden
#
# PCA of feature hashed table (N = 10k features). 

library(data.table)
library(HDF5Array)
# library(gmodels)

#--------------------
# load data -- server
#--------------------
save.dpath <- file.path("home", "metamaden", "bioinfo_appnote", "manuscript_final_results")
# load metadata
md.fname <- "mdf-epic-hm450k-merge_normal-blood_0-0-3.rda"
md <- get(load(file.path(save.dpath, md.fname)))
# get formatted groups
md$group <- ifelse(md$blood.subgroup %in% c("NA", "all", "blood_spot"), "other/NOS",
                    ifelse(grepl("^peripheral.*", md$blood.subgroup), "PBMC", md$blood.subgroup))


#-------------------------------------
# save gseadj data, then get fh tables
#-------------------------------------
# get granges
gr.fname <- "gr-gseadj_h5se_hm450k-epic-merge_0-0-3"
gr.fpath <- newfpath <- file.path(save.dpath, gr.fname)
gr <- HDF5Array::loadHDF5SummarizedExperiment(gr.fpath)

# get bvals
bval <- t(getBeta(gr)); dim(bval) # [1]  12242 442474
# write labels
samp.labels.fpath <- file.path(save.dpath, 
                               "bval-samplabels_blood-groups_2platforms.txt")
write(rownames(bval), sep = " ", file = samp.labels.fpath)

# write dnam bvals
bval.fname <- "bval-gseadj_row-gsm_4-subgroups-2-platforms.csv"
bval.fpath <- file.path(save.dpath, bval.fname)
cnv <- matrix(colnames(bval), nrow = 1)
data.table::fwrite(cnv, file = bval.fpath, sep = ",", append = F, col.names = F, row.names = F)
nrow.chunk <- 2000; indexv <- seq(1, nrow(bval), nrow.chunk)
t1 <- Sys.time()
for(ri in indexv){
  start.index <- ri; end.index <- start.index + nrow.chunk - 1
  end.index <- ifelse(end.index > nrow(bval), nrow(bval), end.index)
  data.table::fwrite(as.matrix(bval[start.index:end.index,]), sep = ",", 
                     file = bval.fpath, append = T, col.names = F, row.names = F)
  message("Finished writing row ", ri, ", time: ", Sys.time()-t1)}
message("done")

# check write file success
# bread <- data.table::fread(bval.fpath, sep = ",", skip = 12240, data.table = F)

# make fh table
# use get_fhtable.py to get feature hashed data

#-----------------------------
# pca, all subgroups -- server
#-----------------------------

# load fh data
fh.fname <- "bval-gseadj-fh10k_all-blood-2-platforms.csv"
fh.fpath <- file.path(save.dpath, fh.fname)
fh <- fread(fh.fpath, sep = ",", header = F, data.table = F)
rownames(fh) <- fh[,1]; fh <- fh[,c(2:10001)]
identical(rownames(md), rownames(fh)) # TRUE

# run pca
pca.dat <- prcomp(fh)
# save pca results
pca.dat.fname <- "pca-results_fh10k-gseadj_2platforms-autodnam.rda"
pca.dat.fpath <- file.path(save.dpath, pca.dat.fname)
save(pca.dat, file = pca.dat.fpath)
message("done!")

# pca dat matrix
mpca <- pca.dat$x
# get anova md
md <- md[rownames(md) %in% rownames(mpca),]
identical(rownames(md), rownames(mpca))
cnv <- c("gse", "predsex", "predage", "group", "platform",
         colnames(md)[grepl("predcell", colnames(md))], 
         "glint.epi.pc1", "glint.epi.pc2")
md.anova <- md[,cnv]
for(c in c(1,2,4,5)){md.anova[,c] <- as.factor(md.anova[,c])}
for(c in c(3,6:13)){md.anova[,c] <- as.numeric(md.anova[,c])}

# run anovas
tdf <- t(mpca[,c(1:10)])
avdf <- do.call(rbind, lapply(seq(nrow(tdf)), function(ri){
  pci <- data.frame(pci = as.numeric(tdf[ri,]))
  rownames(pci) <- colnames(tdf); md.af <- cbind(md.anova, pci)
  stati <- summary(aov(pci ~ gse + predsex + predage + group + platform + predcell.CD8T + 
                         predcell.CD4T + predcell.NK + predcell.Bcell + predcell.Mono + 
                         predcell.Gran + glint.epi.pc1 + glint.epi.pc2, data = md.af))
  varv <- gsub(" ", "", rownames(stati[[1]]))
  lstat <- eval(parse(text = paste0(as.character(stati))))
  mssq.tot <- sum(lstat$`Sum Sq`)
  mssq.percv <- unlist(lapply(lstat$`Sum Sq`, function(mi){100*mi/mssq.tot}))
  dfi <- data.frame(mssq.percv = round(mssq.percv,3), var = varv); dfi$pc <- ri
  return(dfi)}))

# do anovas -- without gse id
avdf.nogse <- do.call(rbind, lapply(seq(nrow(tdf)), function(ri){
  pci <- data.frame(pci = as.numeric(tdf[ri,]))
  rownames(pci) <- colnames(tdf)
  md.af <- cbind(md.anova, pci)
  stati <- summary(aov(pci ~ predsex + predage + group + platform + predcell.CD8T + 
                         predcell.CD4T + predcell.NK + predcell.Bcell + predcell.Mono + 
                         predcell.Gran + glint.epi.pc1 + glint.epi.pc2, data = md.af))
  varv <- gsub(" ", "", rownames(stati[[1]]))
  lstat <- eval(parse(text = paste0(as.character(stati))))
  mssq.tot <- sum(lstat$`Sum Sq`)
  mssq.percv <- unlist(lapply(lstat$`Sum Sq`, function(mi){100*mi/mssq.tot}))
  dfi <- data.frame(mssq.percv = round(mssq.percv,3), var = varv); dfi$pc <- ri
  return(dfi)}))

# do anovas -- without glint epi
avdf.noglint <- do.call(rbind, lapply(seq(nrow(tdf)), function(ri){
  pci <- data.frame(pci = as.numeric(tdf[ri,]))
  rownames(pci) <- colnames(tdf)
  md.af <- cbind(md.anova, pci)
  stati <- summary(aov(pci ~ gse + predsex + predage + group + platform + predcell.CD8T + 
                         predcell.CD4T + predcell.NK + predcell.Bcell + predcell.Mono + 
                         predcell.Gran, data = md.af))
  varv <- gsub(" ", "", rownames(stati[[1]]))
  lstat <- eval(parse(text = paste0(as.character(stati))))
  mssq.tot <- sum(lstat$`Sum Sq`)
  mssq.percv <- unlist(lapply(lstat$`Sum Sq`, function(mi){100*mi/mssq.tot}))
  dfi <- data.frame(mssq.percv = round(mssq.percv,3), var = varv); dfi$pc <- ri
  return(dfi)}))

# save all pca results data
lpca <- list(mpc10 = tdf, pcsd = pca.dat$sdev[1:10], avdf = avdf, 
             avdf.nogse = avdf.nogse, avdf.noglint = avdf.noglint, 
             md = md.anova)
lpca.fname <- "lpca-mpc10-sd-av_fh10k-2plat-autodnam.rda"
lpca.fpath <- file.path(save.dpath, lpca.fname)
save(lpca, file = lpca.fpath)
