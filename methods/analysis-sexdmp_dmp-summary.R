#!/usr/bin/env R

# Author: Sean Maden
#
# Summarize DMP properties from validation of sex DMPs in PBMCs and whole blood.
#
#

#----------
# load data
#----------


# load dmp info
dmpdf <- get(load(file.path("st_sex-dmp_all-dmp-info.rda")))
# format columns
for(c in c(2,3,4,14)){dmpdf[,c] <- as.logical(dmpdf[,c])}
for(c in c(7,8,9)){dmpdf[,c] <- as.numeric(dmpdf[,c])}


# load grant dmps
grant.dmp <- read.csv("grant-2021_wb-sex-dmp.csv")
# bind grant and 
length(intersect(dmpdf$cgid, grant.dmp$Name)) # 272
length(unique(grant.dmp$Name)) # 544
which.int.grantid <- grant.dmp$Name %in% dmpdf$cgid
cg.out <- grant.dmp[!which.int.grantid,]$Name
mdatv <- c(cg.out, rep("NA", length(cg.out)*12))
df.bind <- as.data.frame(matrix(mdatv, ncol = 13))
colnames(df.bind) <- colnames(dmpdf)
dmpdf <- rbind(dmpdf, df.bind)
dmpdf$is.grant.2021.dmp <- ifelse(dmpdf$cgid %in% grant.dmp$Name, TRUE, FALSE)
# format other vars
for(cni in c("is.inoshita.2015.dmp", "is.pbmc.dmp", "is.whole.blood.dmp")){
  dmpdf[,cni] <- ifelse(dmpdf[,cni] == "NA", FALSE, dmpdf[,cni])
}

#------------------------------------------------------------
# get differences by whole blood and pbmc -- compare inoshita
#------------------------------------------------------------
table(dmpdf$is.pbmc.dmp, dmpdf$is.whole.blood.dmp)
#       FALSE TRUE
# FALSE   450  721
# TRUE    721  279

cond <- dmpdf$is.pbmc.dmp==T|dmpdf$is.whole.blood.dmp==T
table(cond)
# FALSE  TRUE 
# 450  1721

which.cond <- which(cond)
ddf <- dmpdf[which.cond,]
table(ddf$is.inoshita.2015.dmp)
# FALSE  TRUE 
# 1609   112
# 112/292 = 0.3835616

# consensus dmps which map to islands
table(ddf[ddf$is.inoshita.2015.dmp,]$relation.to.island)
# Island N_Shore OpenSea S_Shelf S_Shore 
# 54      17      24       4      13
# 54/112 = 0.4821429

cond2 <- as.logical(ddf$is.whole.blood.dmp) & 
  as.logical(ddf$is.inoshita.2015.dmp) &
  !as.logical(ddf$is.pbmc.dmp)
nrow(ddf[cond2,]) # 42

cond2 <- as.logical(ddf$is.pbmc.dmp) & 
  as.logical(ddf$is.inoshita.2015.dmp) &
  !as.logical(ddf$is.whole.blood.dmp)
nrow(ddf[cond2,]) # 17

cond2 <- as.logical(ddf$is.pbmc.dmp) & 
  as.logical(ddf$is.inoshita.2015.dmp) &
  as.logical(ddf$is.whole.blood.dmp)
nrow(ddf[cond2,]) # 53

#------------------------------------------------------------
# get differences by whole blood and pbmc -- compare grants
#------------------------------------------------------------
table(ddf$is.grant.2021.dmp)
# FALSE  TRUE 
# 1471   250
# 250/544 = 0.4595588

table(dmpdf$is.grant.2021.dmp)

#-----------------------------
# get replicated dmp fractions
#-----------------------------
# get fraction of inoshita dmps validated
table(dmpdf$is.inoshita.2015.dmp, dmpdf$is.pbmc.dmp)
#       FALSE TRUE
# FALSE   949  930
# TRUE    222   70
# 70/292 = 0.239726
table(dmpdf$is.inoshita.2015.dmp, dmpdf$is.whole.blood.dmp)
#       FALSE TRUE
# FALSE   974  905
# TRUE    197   95
# 95/292 = 0.32

# get fraction of grant dmps validated
table(dmpdf$is.grant.2021.dmp, dmpdf$is.pbmc.dmp)
#       FALSE TRUE
# FALSE   804  823
# TRUE    367  177
# 177/544 = 0.3253676
table(dmpdf$is.grant.2021.dmp, dmpdf$is.whole.blood.dmp)
#       FALSE TRUE
# FALSE   816  811
# TRUE    355  189
# 189/544 = 0.35

#---------------------------
# dmp all 3 whole blood sets
#---------------------------
cond <- which(dmpdf$is.grant.2021.dmp==T & dmpdf$is.inoshita.2015.dmp==T & 
                dmpdf$is.whole.blood.dmp==T)
nrow(dmpdf[cond,]) # 37
nrow(dmpdf[cond,])/nrow(dmpdf) # 0.017

# dmp overlapping all sets
cond <- which(dmpdf$is.grant.2021.dmp==T & dmpdf$is.inoshita.2015.dmp==T & 
                dmpdf$is.pbmc.dmp==T & dmpdf$is.whole.blood.dmp==T)
nrow(dmpdf[cond,]) # 26
nrow(dmpdf[cond,])/nrow(dmpdf) # 0.01197605

#--------------------------
# compare mean bvals by sex
#--------------------------
# ttest validation results
tdf.wb.fname <- "ttest-df-results-mvalfit_wholeblood-2platforms_inoshita-2015-validate.rda"
tdf.pbmc.fname <- "ttest-df-results-mvalfit_pbmc-2platforms_inoshita-2015-validate.rda"
tdf.wb <- get(load(tdf.wb.fname))
tdf.pbmc <- get(load(tdf.pbmc.fname))

# filter dmps
dff <- dmpdf[as.logical(dmpdf$is.inoshita.2015.dmp),]
tdf.wb <- tdf.wb[tdf.wb$cgid %in% dff[as.logical(dff$is.whole.blood.dmp),]$cgid,]
tdf.pbmc <- tdf.pbmc[tdf.pbmc$cgid %in% dff[as.logical(dff$is.pbmc.dmp),]$cgid,]

# summarize mean bvals
table(tdf.wb$mean.male > tdf.wb$mean.female)
# FALSE  TRUE 
# 77    18 

table(tdf.pbmc$mean.male > tdf.pbmc$mean.female)
# FALSE  TRUE 
# 50    20

# compare directions -- pbmc
dff <- dmpdf[dmpdf$is.pbmc.dmp,]
table(dff$dnam.direction.inoshita.2015.dmp > 0, dff$dnam.direction.pbmc.dmp > 0)
#       FALSE TRUE
# FALSE   683  297
# TRUE      0   20
# (683+20)/(683+20+297) = 0.703
# compare directions -- wb
dff <- dmpdf[dmpdf$is.whole.blood.dmp,]
table(dff$dnam.direction.inoshita.2015.dmp > 0, dff$dnam.direction.whole.blood.dmp > 0)
#       FALSE TRUE
# FALSE   934   48
# TRUE      0   18
# (934+18)/nrow(dff) = 0.952