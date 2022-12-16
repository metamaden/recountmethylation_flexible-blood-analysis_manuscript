#!/usr/bin/env R

#-----
# test
#-----
libv <- c("minfi", "FlowSorted.CordBlood.450k","minfiData")
sapply(libv, library, character.only = T)
# load example data
rg <- get(data(RGsetEx))
# get cell types
cb <- get(load("FlowSorted.CordBlood.450k"))
ctv <- unique(cb$cellTypes)
# get estimates
est <- estimateCellCounts(rg, compositeCellType = "CordBlood", cellTypes = ctv)


# load data
# get rgset data
#rg.fname <- "remethdb_epic-hm850k_h5se_rg_1669220613_0-0-3"
#rg.fpath <- file.path("recount-methylation-files", "compilations", rg.fname)
#rg <- HDF5Array::loadHDF5SummarizedExperiment(rg.fpath)

# get adj gr data
gr.fname <- "gr-gseadj_h5se_hm450k-epic-merge_0-0-3"
gr <- loadHDF5SummarizedExperiment(gr.fname)
dim(gr) # 442474  12242

# get gsm ids
csv.fname <- "table-s1_md-normal-blood.csv"
# csv.fpath <- file.path("recountmethylation_v2_manuscript", "data", "tables", csv.fname)
csv <- read.csv(csv.fname)

# get predictions
# cord blood gsm idv
idv <- as.character(csv[csv$blood.subgroup == "cord_blood",]$gsm)
length(idv) # 1475
gsmv <- as.character(gsub("\\..*", "", colnames(gr)))
length(intersect(gsmv, idv)) # 1475
# get matrix-backed se
gr.sub <- gr[,gr$gsm %in% idv]
dim(gr.sub)
gr.sub.se <- GenomicRatioSet(gr = granges(gr.sub),
                             Beta = as.matrix(getBeta(gr.sub)),
                             annotation = annotation(gr.sub))

# get cell types
cb <- get(data("FlowSorted.CordBlood.450k"))
ctv <- unique(cb$cellTypes)
# get estimates
est <- estimateCellCounts(gr.sub.se, compositeCellType = "CordBlood", cellTypes = ctv)



