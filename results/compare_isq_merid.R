#!/usr/bin/env R

# Author: Sean Maden
#
# Compare results from prior work using new cord blood DMP sets.
#

library(ggplot2)
library(gridExtra)

#----------
# load data
#----------
st.fname <- "st_gest-age-dmp_all-dmp-info.rda"
dmpdf <- get(load(file.path(st.fname)))

dmp1.fname <- "dmps1-nocomplications_cord-validation_merid-et-al-2020.csv"
merid1 <- read.csv(dmp1.fname)

#-----------------------------
# compare i-squared properties
#-----------------------------
# compare isq, nocomp
dff <- dmpdf[grepl("no_complications", dmpdf$merid.type),]
dff <- dff[order(match(dff$cgid, merid1$ï..CpGID)),]
identical(dff$cgid, merid1$ï..CpGID)
cor.test(as.numeric(dff$cord.blood.ttests.padjusted), 
         as.numeric(merid1$I_SQUARE), method = "spearman")
# pval = 3.999e-13
# rho = -0.2449905

cor.test(abs(as.numeric(dff$cord.blood.ttest.estimate)), 
         as.numeric(merid1$I_SQUARE), method = "spearman")
# p-value = 8.564e-14
# rho = 0.2517163 

table(dff$is.cord.blood.dmp, merid1$I_SQUARE..50)
#         no  yes
# FALSE 7857  189
# TRUE   723  130

# violin plots
dfp <- data.frame(isq = as.numeric(merid1$I_SQUARE))
dfp$type <- "merid1"
is.cb.dmp <- which(dff$is.cord.blood.dmp==T)
dfp <- rbind(dfp, 
             data.frame(isq = as.numeric(merid1$I_SQUARE[is.cb.dmp]),
                        type = rep("cord_blood", length(is.cb.dmp))))

t.test(dfp[dfp$type=="merid1",]$isq,
       dfp[dfp$type=="cord_blood",]$isq)

ggplot(dfp, aes(x = type, y = isq)) + geom_violin(draw_quantiles = 0.5)

# scatterplot
dff$isq <- as.numeric(merid1$I_SQUARE)

plot(as.numeric(dff$cord.blood.ttest.estimate), dff$isq)

#----------------------------
# compare coefficients, pvals
#----------------------------
dff$coef.re <- merid1$Coefficient..RE
dff$coef.fe <- merid1$Coefficient..FE
dff$pval.fe <- merid1$PVALUE.FE
dff$pval.re <- merid1$PVALUE.RE
dff$neg_log10_pval_fe <- -1*log10(dff$pval.fe)
dff$neg_log10_pval_re <- -1*log10(dff$pval.re)

# get plot data
p1 <- ggplot(dff, aes(x = coef.fe, color = is.cord.blood.dmp)) + 
  geom_density() + theme_bw() + theme(legend.position = "none")

p2 <- ggplot(dff, aes(x = coef.re, color = is.cord.blood.dmp)) + 
  geom_density() + theme_bw() + theme(legend.position = "none")

p3 <- ggplot(dff, aes(x = neg_log10_pval_fe, color = is.cord.blood.dmp)) + 
  geom_density() + theme_bw() + theme(legend.position = "none")

p4 <- ggplot(dff, aes(x = neg_log10_pval_re, color = is.cord.blood.dmp)) + 
  geom_density() + theme_bw() + theme(legend.position = c(0.6, 0.5))

# save new supp fig pdf
pdf("sfig_merid-stats-by-cbdmp.pdf", width = 8, height = 5)
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
dev.off()
