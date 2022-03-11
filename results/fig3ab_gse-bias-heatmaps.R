#!/usr/bin/env R

# Author: Sean Maden
#
# Make heatmaps of GSE bias simulations results.

library(ggplot2); library(data.table)
library(scales); library(gridExtra)

#----------
# load data
#----------
# load data tables
msq.fname <- "msq-gse-bias_all-blood-2-platforms.rda"
msq <- get(load(msq.fname))

mdif.fname <- "mdiff-gse-bias_all-blood-2-platforms.rda"
mdif <- get(load(mdif.fname))

#-----------------
# helper functions
#-----------------
get_medians <- function(mdat, groupv){
  mmed <- do.call(rbind, lapply(groupv, function(typei){
    # get all merged medians
    col.str <- paste0(c(typei, "group"), collapse = "|")
    mff <- mdat[,c(grepl(col.str, colnames(mdat)))]
    statv <- apply(mff[,1:10], 2, median, na.rm = T)
    statv.format <- format(statv, scientific = T, digits = 3)
    med.merge <- c(statv.format, "merged_groups")
    # medians by sample group
    med.groups <- do.call(rbind, lapply(unique(mff$group), function(groupi){
      statv <- apply(mff[mff$group == groupi,1:10], 2, median, na.rm = T)
      statv.format <- format(statv, scientific = T, digits = 3)
      c(statv.format, groupi)}))
    med.all <- rbind(med.merge, med.groups)
    groupv <- med.all[,11]; med.all <- med.all[,c(1:10)]
    variablev <- gsub("_.*", "", colnames(med.all))
    typev <- gsub(".*_", "", colnames(med.all))
    med.all.format <- as.data.frame(t(med.all[,c(1:10)]), 
                                    stringsAsFactors = F)
    colnames(med.all.format) <- groupv
    med.all.format$variable <- varv
    med.all.format$model.type <- typev
    return(med.all.format)
  })); return(mmed)
}

get_variances <- function(mdat, groupv){
  mvar <- do.call(rbind, lapply(groupv, function(typei){
    # get all merged medians
    col.str <- paste0(c(typei, "group"), collapse = "|")
    mff <- mdat[,c(grepl(col.str, colnames(mdat)))]
    statv <- apply(mff[,1:10], 2, var, na.rm = T)
    statv.format <- format(statv, scientific = T, digits = 3)
    var.merge <- c(statv.format, "merged_groups")
    # medians by sample group
    var.groups <- do.call(rbind, lapply(unique(mff$group), function(groupi){
      statv <- apply(mff[mff$group == groupi,1:10], 2, var, na.rm = T)
      statv.format <- format(statv, scientific = T, digits = 3)
      c(statv.format, groupi)}))
    var.all <- rbind(var.merge, var.groups)
    groupv <- var.all[,11]; var.all <- var.all[,c(1:10)]
    variablev <- gsub("_.*", "", colnames(var.all))
    typev <- gsub(".*_", "", colnames(var.all))
    var.all.format <- as.data.frame(t(var.all[,c(1:10)]), stringsAsFactors = F)
    colnames(var.all.format) <- groupv; var.all.format$variable <- varv
    var.all.format$model.type <- typev; return(var.all.format)
  })); return(mvar)
}

get_mstat <- function(mdat, groupv){
  mmed <- get_medians(mdat, groupv); mvar <- get_variances(mdat, groupv)
  mmed$stat <- "median"; mvar$stat <- "variance"; mstat <- rbind(mmed, mvar)
  return(mstat)
}

#----------------------------------------
# fraction explained variances plot data
#----------------------------------------
dfp.fev <- apply(msq[,c(1:39)], 2, function(ci){
  median(as.numeric(ci), na.rm=T)})
# format heatmap data
lvlv <- c("gse", "predsex", "predcell.Mono", "predcell.NK",
          "predcell.CD4T", "predage", "predcell.Bcell", "predcell.CD8T", 
          "predcell.Gran", "platform", "glint.epi.pc2", "glint.epi.pc1", 
          "Residuals")
dfp.fev <- data.frame(var = names(dfp.fev), value = as.numeric(dfp.fev))
dfp.fev$model <- gsub(".*_", "", dfp.fev$var)
dfp.fev$var <- gsub("_.*", "", dfp.fev$var)
dfp.fev$`Median\nFEV` <- as.numeric(dfp.fev$value)
dfp.fev$var <- factor(dfp.fev$var, levels = lvlv)
dfp.fev$model <- factor(dfp.fev$model, 
                        levels = c("unadj", "adj1", "adj2"))
dfp.fev$value.label <- round(100*dfp.fev$value, digits = 2)

#---------------------------------
# median abs differences plot data
#---------------------------------
dfp.mad <- apply(mdif[,c(1:39)], 2, function(ci){
  median(abs(as.numeric(ci)), na.rm = T)})
# format dfp.mad for heatmap
dfp.mad <- data.frame(var = names(dfp.mad), value = as.numeric(dfp.mad))
dfp.mad$difftype <- gsub(".*_", "", dfp.mad$var)
dfp.mad$var <- gsub("_.*", "", dfp.mad$var)
lvlv <- c("gse", "predsex", "predcell.Mono", "predcell.NK",
          "predcell.CD4T", "predage", "predcell.Bcell", "predcell.CD8T", 
          "predcell.Gran", "platform", "glint.epi.pc2", "glint.epi.pc1", 
          "Residuals")
dfp.mad$var <- factor(dfp.mad$var, levels = lvlv)
# make new difflab var
labv <- c("unadj-adj1", "unadj-adj2", "adj1-adj2")
dfp.mad$diff.lab <- ifelse(dfp.mad$difftype == "diff1", labv[1],
                           ifelse(dfp.mad$difftype == "diff2", labv[2],
                                  labv[3]))
dfp.mad$diff.lab <- factor(dfp.mad$diff.lab, levels = labv)

#--------------------------
# ggtile heatmap fev, fig2b
#--------------------------
# get plot object
ylab.str <- paste0(paste0(rep("\n", 3), collapse = ""), 
                   "Model type", collapse = "")
ggtile <- ggplot(dfp.fev, aes(x = var, y = model)) + 
  geom_tile(aes(fill = `Median\nFEV`)) + 
  theme_bw() + geom_text(aes(label = value.label), 
                         color = "white", size = 3) + 
  xlab("Variable") + ylab(ylab.str) +
  scale_fill_gradient(low = "blue", high = "red", name = "Median\nFEV", 
                      labels = scales::label_scientific()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# save new plot
pdf.fname <- "fig2b_heatmap-fev-all_2platforms.pdf"
pdf(pdf.fname, 7.5, 2.5); print(ggtile); dev.off()

#-----------------------------------
# get the heatmap of mads -- fig 2c
#-----------------------------------
# get tile labels
dfp.mad$labels <- round(100*dfp.mad$value, 2)
ggtile <- ggplot(dfp.mad, aes(x = var, y = diff.lab)) +
  geom_tile(aes(fill = value), color = "white") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) +
  scale_fill_gradient(low = "blue", high = "red", name = "MAD", 
                      labels = scales::label_scientific()) +
  ylab("Model difference") + xlab("Variable") +
  geom_text(aes(label = labels), color = "white", size = 3)

# save plot
plot.fname <- "fig2c_ggtile_mads-labeled_gse-bias.pdf"
pdf(plot.fname, 7.5, 2.5); print(ggtile);dev.off()

#---------------------------------------------------------
# get composite heatmaps of each subgroup type -- supp fig
#---------------------------------------------------------
pdf.fname <- "sfig_mad-bytx_gse-bias-results.pdf"
dfp <- mmad[!mmad$tx == "combined",]
ltile <- lapply(unique(dfp$tx), function(txi){
  dfpi <- dfp[dfp$tx == txi,]
  # order vars
  medv <- unlist(lapply(unique(dfpi$var), 
                        function(vari){median(dfpi[dfpi$var == vari,]$mad)}))
  dfpi$var <- factor(dfpi$var, levels = unique(dfpi$var)[order(medv)])
  # get tile labels
  dfpi$labels <- format(dfpi$mad, scientific = T, digits = 3)
  ggtile <- ggplot(dfpi, aes(x=var, y=diff)) +
    geom_tile(aes(fill = mad), color = "white") + theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    scale_fill_gradient(low = "blue", high = "red", name = "MAD", 
                        labels = scales::label_scientific()) +
    ggtitle(txi) + geom_text(aes(label = labels), color = "white", size = 3)
  return(ggtile)})
# make composite plot
pdf(pdf.fname, 8, 8.3)
grid.arrange(ltile[[1]], ltile[[2]], ltile[[3]], ltile[[4]], 
             nrow = 4, bottom = "Variable", 
             left = "Difference type")
dev.off()

#----------------
# get supp tables
#----------------
# format table
# mmad <- dfp.mad
# save table
# save.dpath <- "."
# mmad.fname <- "supptable_mads_gse-bias-sims"
# save(mmad, file = file.path(save.dpath, paste0(mmad.fname,".rda")))
# write.csv(mmad, file = file.path(save.dpath, paste0(mmad.fname,".csv")))
# write.table(mmad, file = file.path(save.dpath, paste0(mmad.fname,".tsv")))