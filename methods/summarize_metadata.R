#!/usr/bin/env R

# Author: Sean Maden
# 
# Summarize sample metadata

library(gridExtra); library(ggpubr); library(ggsci)

#----------
# load data
#----------
# load metadata
md.fname <- "si2_blood-md-2platforms.rda"
md <- get(load(md.fname))
# format groups
md$subgroup <- md$blood_subgroup
md[grepl("^periph.*", md$subgroup),]$subgroup <- "PBMC"
md[grepl("NA|all", md$subgroup),]$subgroup <- "other/NOS"

#-------------
# make table 1
#-------------
t1 <- do.call(rbind, lapply(c("cord_blood", "whole_blood", "PBMC", "all"), function(groupi){
  if(groupi=="all"){mdf <- md} else{mdf <- md[md$subgroup==groupi,]}
  fract.male <- round(nrow(mdf[mdf$predsex == "F",])/nrow(mdf),2)
  mean.age <- round(mean(as.numeric(mdf$predage), na.rm = T), 2)
  sd.age <- round(sd(as.numeric(mdf$predage), na.rm = T),2)
  num.gse <- length(unique(mdf$gse)); num.gsm <- length(unique(mdf$gsm))
  num.gse.hm450k <- length(unique(mdf[mdf$platform=="hm450k",]$gse))
  num.gse.epic <- length(unique(mdf[mdf$platform=="epic",]$gse))
  num.gsm.hm450k <- length(unique(mdf[mdf$platform=="hm450k",]$gsm))
  num.gsm.epic <- length(unique(mdf[mdf$platform=="epic",]$gsm))
  mdat <- c(groupi, num.gse, num.gse.hm450k, num.gse.epic, 
            num.gsm, num.gsm.hm450k, num.gsm.epic, 
            fract.male, mean.age, sd.age)
  return(matrix(mdat, nrow = 1))}))
colnames(t1) <- c("Subgroup", "Num_studies", "Num_studies_HM450K", "Num_studies_EPIC", 
                  "Num_samples", "Num_samples_HM450K", "Num_samples_EPIC", 
                  "Fract_female", "Mean_age", "SD_age")
t1 <- as.data.frame(t(t1), stringsAsFactors = F)
colnames(t1) <- t1[1,]; t1 <- t1[c(2:nrow(t1)),]
t1 <- t1[,c(4,2,1,3)]
# save table 1
write.csv(t1, file = "table1_demo-4groups-2platforms.csv")

#-------------------------
# make barplot dfp objects
#-------------------------
# make dfp, samples
dfp <- as.data.frame(table(md$subgroup))
colnames(dfp) <- c("subgroup", "samples")
dfp$platform <- "Combined"
dfp <- dfp[,c(1,3,2)]
dfp.platform <- as.data.frame(table(md$subgroup, md$platform))
colnames(dfp.platform) <- c("subgroup", "platform", "samples")
dfp <- rbind(dfp, dfp.platform)
dfp$platform <- as.character(dfp$platform)
dfp[dfp$platform == "hm450k",]$platform <- "HM450k"
dfp[dfp$platform == "epic",]$platform <- "EPIC"
dfp$platform <- factor(dfp$platform, levels = c("Combined", "HM450k", "EPIC"))
dfp$subgroup <- as.character(dfp$subgroup)
dfp$subgroup <- factor(dfp$subgroup, 
                    levels = c("PBMC", "cord_blood", 
                               "other/NOS", "whole_blood"))
colnames(dfp) <- c("Subgroup", "Platform", "Samples")
dfp.samples <- dfp[!is.na(dfp[,1]),]

# make dfp, studies
mdff <- md[!duplicated(md$gse),]
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
dfp$group <- factor(dfp$group, 
                    levels = c("PBMC", "cord_blood", 
                               "other/NOS", "whole_blood"))
colnames(dfp) <- c("Subgroup", "Platform", "Studies")
dfp.studies <- dfp[!is.na(dfp[,1]),]

#--------------------------------
# make composite barplots, fig 1a
#--------------------------------
# palette
pal <- c("#CD534CFF", "#EFC000FF", "#3B3B3BFF", "#868686FF")


# get fill pattern variable
dfp.samples$fillvar <- dfp.studies$fillvar <- ifelse(dfp$Platform == "Combined", "a", 
                                                     ifelse(dfp$Platform == "HM450k", 
                                                            "b", "c"))
# samples barplot object
bp.samples.comp <- ggplot(dfp.samples, 
                          aes(x = Platform, y = Samples, fill = Subgroup)) + 
  geom_bar(stat = "identity") + theme_bw() + scale_fill_manual(values = pal) + xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 0.95),
        legend.position = "none", plot.margin=unit(c(0.01,0,0.1,0.1),"cm"))
  
# studies barplot object
bp.studies.comp <- ggplot(dfp.studies, 
                          aes(x = Platform, y = Studies, fill = Subgroup)) + 
  geom_bar(stat = "identity") + theme_bw() + xlab("") + scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust =0.95),
        legend.position = "none", plot.margin=unit(c(0.01,0.1,0.1,0.4),"cm"))
  
# store legend
bp.leg <- ggplot(dfp.studies, aes(x = Platform, y = Studies, fill = Subgroup)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = pal)
legend <- as_ggplot(ggpubr::get_legend(bp.leg))

# define the layout matrix
pdf.fname <- "fig1a_bp-composite_samples-studies_blood-groups_2platforms.pdf"
pdf(pdf.fname, 6, 2)
grid.arrange(bp.samples.comp, bp.studies.comp, legend, 
             layout_matrix = matrix(c(1,2,3), ncol = 3), 
             bottom = paste0("Platform", paste0(rep(" ", 15), 
                                                collapse = "")))
dev.off()









