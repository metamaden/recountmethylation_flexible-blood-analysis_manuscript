#!/usr/bin/env R

# Author: Sean Maden
# 
# Summarize sample metadata

#----------
# load data
#----------
md.all.fpath <- "si1_all-md-2platforms.rda"
md.all <- get(load(md.all.fpath))
sit.fpath <- "sitest_results_blood-groups-2platforms.csv"
sit <- data.table::fread(sit.fpath, sep = ",", header = T, data.table = F)
sit <- sit[,c(2:ncol(sit))]
colnames(sit) <- c("gsm", paste0("k=", c(10, 100, 500, 1000)), 
                   "subgroup", "platform")

#-------------
# analyze data
#-------------
sif1 <- sit[!sit$subgroup == "all",]
sif2 <- sit[sit$subgroup == "all" & !sit$gsm %in% sif1$gsm,]
sif <- rbind(sif1, sif2)
dim(sif)

mdf <- md.all[md.all$gsm %in% sif$gsm,]
mdf <- mdf[order(match(mdf$gsm, sif$gsm)),]
identical(mdf$gsm, sif[,1])

mdf$subgroup <- sif$subgroup
dim(mdf)
length(unique(mdf$gse))
table(mdf$platform)

table(sif$subgroup, sif$platform)
table(mdf$gse, sif$subgroup, mdf$platform)

#-------------------------
# make barplot dfp objects
#-------------------------
# make dfp, samples
dfp <- as.data.frame(table(mdf$subgroup))
colnames(dfp) <- c("group", "samples"); dfp$platform <- "Combined"
dfp <- dfp[,c(1,3,2)]
dfp.platform <- as.data.frame(table(mdf$subgroup, mdf$platform))
colnames(dfp.platform) <- c("group", "platform", "samples")
dfp <- rbind(dfp, dfp.platform)
dfp$platform <- as.character(dfp$platform)
dfp[dfp$platform == "hm450k",]$platform <- "HM450k"
dfp[dfp$platform == "epic",]$platform <- "EPIC"
dfp$platform <- factor(dfp$platform, levels = c("Combined", "HM450k", "EPIC"))
dfp$group <- as.character(dfp$group)
dfp[grepl("^periph.*", dfp$group),]$group <- "PBMC"
dfp[grepl("^all.*", dfp$group),]$group <- "other/NOS"
dfp$group <- factor(dfp$group, 
                    levels = c("PBMC", "cord_blood", 
                               "other/NOS", "whole_blood"))
colnames(dfp) <- c("Blood group", "Platform", "Samples")
dfp.samples <- dfp

# make dfp, studies
mdff <- mdf[!duplicated(mdf$gse),]
dfp <- as.data.frame(table(mdff$subgroup))
colnames(dfp) <- c("group", "studies"); dfp$platform <- "Combined"
dfp <- dfp[,c(1,3,2)]
dfp.platform <- as.data.frame(table(mdff$subgroup, mdff$platform))
colnames(dfp.platform) <- c("group", "platform", "studies")
dfp <- rbind(dfp, dfp.platform)
dfp$platform <- as.character(dfp$platform)
dfp[dfp$platform == "hm450k",]$platform <- "HM450k"
dfp[dfp$platform == "epic",]$platform <- "EPIC"
dfp$platform <- factor(dfp$platform, levels = c("Combined", "HM450k", "EPIC"))
dfp$group <- as.character(dfp$group)
dfp[grepl("^periph.*", dfp$group),]$group <- "PBMC"
dfp[grepl("^all.*", dfp$group),]$group <- "other/NOS"
dfp$group <- factor(dfp$group, 
                    levels = c("PBMC", "cord_blood", 
                               "other/NOS", "whole_blood"))
colnames(dfp) <- c("Blood group", "Platform", "Studies")
dfp.studies <- dfp

#------------------------------------------
# make metadata stacked barplots, fig 1a-1b
#------------------------------------------
library(ggplot2)

# define color blind palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# samples, single plot
bp.samples <- ggplot(dfp.samples, aes(x = Platform, y = Samples, fill = `Blood group`)) + 
  geom_bar(stat = "identity") + theme_bw() + scale_fill_manual(values = cbp1) +
  theme(axis.text.x = element_text(angle = 90))
pdf("bp-samples_blood-groups_2platforms.pdf", 3, 2)
print(bp.samples); dev.off()

# studies, single plot
bp.studies <- ggplot(dfp.studies, aes(x = Platform, y = Studies, fill = `Blood group`)) + 
  geom_bar(stat = "identity") + theme_bw() + scale_fill_manual(values = cbp1) +
  theme(axis.text.x = element_text(angle = 90))
pdf("bp-studies_blood-groups_2platforms.pdf", 3, 2)
print(bp.studies); dev.off()

#--------------------------------
# make composite barplots, fig 1a
#--------------------------------
library(gridExtra); library(cowplot)
# devtools::install_github("https://github.com/coolbutuseless/ggpattern")
library(ggpattern)

legend <- get_legend(bp.studies) # store legend

# get fill pattern variable
dfp.samples$fillvar <- dfp.studies$fillvar <- ifelse(dfp$Platform == "Combined", "a", 
                                                     ifelse(dfp$Platform == "HM450k", 
                                                            "b", "c"))

bp.samples.comp <- ggplot(dfp.samples, 
                          aes(x = Platform, y = Samples, fill = `Blood group`)) + 
  geom_bar(stat = "identity") + theme_bw() + xlab("") +
  scale_x_discrete(labels= rep("", 3)) + scale_fill_manual(values = cbp1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 0.95),
        legend.position = "none", plot.margin=unit(c(0.27,0.1,0.1,0.15),"cm"))

bp.studies.comp <- ggplot(dfp.studies, 
                          aes(x = Platform, y = Studies, fill = `Blood group`)) + 
  geom_bar(stat = "identity") + theme_bw() + xlab("Platform") +
  scale_fill_manual(values = cbp1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust =0.95),
        legend.position = "none", plot.margin=unit(c(0.01,0.1,0.25,0.5),"cm"))
  
# define the layout matrix
lm <- matrix(c(rep(c(1,1,1,1,1,2,2,2,2,2,2), 2),3,3,3,3,3,3,3,3,3,3,3), ncol = 3)
pdf.fname <- "fig1a_bp-composite_samples-studies_blood-groups_2platforms.pdf"
pdf(pdf.fname, 4, 4.5)
grid.arrange(bp.samples.comp, bp.studies.comp, legend, 
             # bottom = paste0("Platform", rep(" ", 10)),
             layout_matrix = lm)
dev.off()









