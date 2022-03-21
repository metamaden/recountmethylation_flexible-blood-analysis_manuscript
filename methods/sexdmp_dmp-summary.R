# load grant dmps
grant.dmp <- read.csv("grant-2021_wb-sex-dmp.csv")

# load dmp info
dmpdf <- get(load(file.path("st_sex-dmp_all-dmp-info.rda")))

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

# get differences by whole blood and pbmc
table(dmpdf$is.pbmc.dmp, dmpdf$is.whole.blood.dmp)
cond.test <- which(dmpdf$is.pbmc.dmp==T|dmpdf$is.whole.blood.dmp==T)
nrow(dmpdf[cond.test,]) # 1739
261/nrow(dmpdf[cond.test,]) # 0.15

# get replicated dmp fractions
table(dmpdf$is.inoshita.2015.dmp, dmpdf$is.pbmc.dmp)
80/292 # 0.27
table(dmpdf$is.inoshita.2015.dmp, dmpdf$is.whole.blood.dmp)
95/292 #  0.32
table(dmpdf$is.grant.2021.dmp, dmpdf$is.pbmc.dmp)
166/544 # 0.31
table(dmpdf$is.grant.2021.dmp, dmpdf$is.whole.blood.dmp)
189/544 # 0.35

# dmp all 3 whole blood sets
cond <- which(dmpdf$is.grant.2021.dmp==T & 
                dmpdf$is.inoshita.2015.dmp==T & 
                dmpdf$is.whole.blood.dmp==T)
nrow(dmpdf[cond,]) # 37
nrow(dmpdf[cond,])/nrow(dmpdf) # 0.017

# dmp overlapping all sets
cond <- which(dmpdf$is.grant.2021.dmp==T & 
                dmpdf$is.inoshita.2015.dmp==T & 
                dmpdf$is.pbmc.dmp==T & 
                dmpdf$is.whole.blood.dmp==T)
nrow(dmpdf[cond,]) # 26
nrow(dmpdf[cond,])/nrow(dmpdf) # 0.012