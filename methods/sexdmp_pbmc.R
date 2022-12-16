#!/usr/bin/env R

# Author: Sean Maden
#
# Validation of Inoshita et al 2015, sex differences in blood using whole blood
# compiled from public HM450K and EPIC data.
#
#

# get dependencies for se objects and plots
library(minfi); library(sva); library(gridExtra); 
library(HDF5Array); library(gaston); library(ggrepel)
library(methyPre)
library(data.table)

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
data(chen_crxcg)
cgfilt <- which(!rownames(gr) %in% chen.crxcg)
gr <- gr[cgfilt,]
dim(gr) # 416045   6866