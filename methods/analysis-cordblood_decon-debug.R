#!/usr/bin/env R

# Author: Sean Maden
#
# Use an updated deconvolution function to estimate cord blood cell types, 
# including enucleated red blood cells.
#

libv <- c("minfi", "minfiData", "FlowSorted.CordBlood.450k")
sapply(libv, library, character.only = TRUE)

#----------
# load data
#----------
rgcb <- get(data(FlowSorted.CordBlood.450k))
rgSet <- get(data(RGsetEx))

est <- estimateCellCounts(rgSet, cellTypes = unique(rgcb$CellType),
                   compositeCellType = "CordBlood") # works

library(FlowSorted.Blood.EPIC)
est2 <- estimateCellCounts2(rgSet, cellTypes = unique(rgcb$CellType),
                            compositeCellType = "CordBlood")


#url <- "https://github.com/metamaden/minfi/blob/master/R/estimateCellCounts.R"
#download.file(url, destfile = fn)

library(minfi)
fn <- "estimateCellCounts.R"
source(fn)
est <- minfi::estimateCellCounts(rgSet, cellTypes = unique(rgcb$CellType),
                          compositeCellType = "CordBlood") # works

#----------------------------
# unpack estimateCellCounts()
#----------------------------
compositeCellType = "Blood"
processMethod = "auto"
probeSelect = "auto"
cellTypes = unique(dfcb$CellType) # c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran")
referencePlatform = c("IlluminaHumanMethylation450k") 
returnAll = FALSE
meanPlot = FALSE
verbose = TRUE

# run checks
minfi:::.isMatrixBackedOrStop(rgSet, "estimateCellCounts")
minfi:::.isRGOrStop(rgSet)
rgSet <- as(rgSet, "RGChannelSet")
referencePlatform <- match.arg(referencePlatform)
rgPlatform <- sub("IlluminaHumanMethylation", "", 
                  annotation(rgSet)[which(names(annotation(rgSet)) == "array")])
platform <- sub("IlluminaHumanMethylation", "", referencePlatform)


referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
subverbose <- max(as.integer(verbose) - 1L, 0L)
data(list = referencePkg)
referenceRGset <- dfcb # get(referencePkg)


if ((processMethod == "auto") && (compositeCellType %in% c("Blood", "DLPFC"))) {
  processMethod <- "preprocessQuantile"
}
if ((processMethod == "auto") && (!compositeCellType %in% c("Blood", "DLPFC"))) {
  processMethod <- "preprocessNoob"
}

processMethod <- get(processMethod)
if ((probeSelect == "auto") && (compositeCellType == "CordBlood")) {
  probeSelect <- "any"
}
if ((probeSelect == "auto") && (compositeCellType != "CordBlood")) {
  probeSelect <- "both"
}


# bind pheno data
newpd <- DataFrame(sampleNames = c(colnames(rgSet), colnames(referenceRGset)), 
                   studyIndex = rep(x = c("user", "reference"), 
                                    times = c(ncol(rgSet), ncol(referenceRGset))), 
                   stringsAsFactors = FALSE)
referencePd <- colData(referenceRGset)
combinedRGset <- combineArrays(object1 = rgSet, 
                               object2 = referenceRGset, 
                               outType = "IlluminaHumanMethylation450k")
colData(combinedRGset) <- newpd
colnames(combinedRGset) <- newpd$sampleNames
rm(referenceRGset)

# do preprocessNoob for cord blood
# combinedMset <- processMethod(combinedRGset, verbose = subverbose)
combinedMset <- preprocessNoob(combinedRGset)
#compTable <- get(paste0(referencePkg, ".compTable"))
#combinedMset <- combinedMset[which(rownames(combinedMset) %in% rownames(compTable)), ]
#rm(combinedRGset)

# get separate Msets
referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
colData(referenceMset) <- as(referencePd, "DataFrame")
mSet <- combinedMset[, combinedMset$studyIndex == "user"]
colData(mSet) <- as(colData(rgSet), "DataFrame")
rm(combinedMset)
mSet.subj <- mSet

# pick probes for decon
#compData <- minfi:::pickCompProbes(mSet = referenceMset,
#                                   cellTypes = cellTypes,
#                                   compositeCellType = compositeCellType,
#                                   probeSelect = probeSelect)

#------------------------
# unpack pickCompProbes()
#------------------------
mSet = referenceMset
cellTypes = cellTypes
compositeCellType = compositeCellType
probeSelect = probeSelect
numProbes = 50

minfi:::.isMatrixBackedOrStop(mSet)
splitit <- function(x) {split(seq_along(x), x)}
  
# subset bvals
p <- getBeta(mSet)
pd <- as.data.frame(colData(mSet))
if (!is.null(cellTypes)) {
  if (!all(cellTypes %in% pd$CellType)) 
      stop("elements of argument 'cellTypes' is not part of ", 
           "'mSet$CellType'")
  keep <- which(pd$CellType %in% cellTypes)
  pd <- pd[keep, ]
  p <- p[, keep]
}


pd$CellType <- factor(pd$CellType, levels = cellTypes)
ffComp <- rowFtests(p, pd$CellType)
prof <- vapply(X = splitit(pd$CellType), 
               FUN = function(j) rowMeans2(p, cols = j), 
               FUN.VALUE = numeric(nrow(p)))
r <- rowRanges(p)
compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[, 2]))
names(compTable)[1] <- "Fstat"
names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- 
  c("low", "high", "range")
tIndexes <- splitit(pd$CellType)


tstatList <- lapply(tIndexes, function(i) {
  x <- rep(0, ncol(p))
  x[i] <- 1
  return(rowttests(p, factor(x)))
})

if (probeSelect == "any") {
  probeList <- lapply(tstatList, function(x) {
    y <- x[x[, "p.value"] < 1e-08, ]
    yAny <- y[order(abs(y[, "dm"]), decreasing = TRUE), 
    ]
    c(rownames(yAny)[seq(numProbes * 2)])
  })
}
else {
  probeList <- lapply(tstatList, function(x) {
    y <- x[x[, "p.value"] < 1e-08, ]
    yUp <- y[order(y[, "dm"], decreasing = TRUE), ]
    yDown <- y[order(y[, "dm"], decreasing = FALSE), 
    ]
    c(rownames(yUp)[seq_len(numProbes)], rownames(yDown)[seq_len(numProbes)])
  })
}

trainingProbes <- unique(unlist(probeList))
trainingProbes <- trainingProbes[!is.na(trainingProbes)]
p <- p[trainingProbes, ]
pMeans <- colMeans2(p)
names(pMeans) <- pd$CellType
form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse = "+")))
phenoDF <- as.data.frame(model.matrix(~pd$CellType - 1))
colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
if (ncol(phenoDF) == 2) {
  X <- as.matrix(phenoDF)
  coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
} else {
    tmp <- minfi:::validationCellType(Y = p, pheno = phenoDF, modelFix = form)
    coefEsts <- tmp$coefEsts
}

compData <- list(coefEsts = coefEsts, compTable = compTable, sampleMeans = pMeans)

#-----------
# get counts
#-----------
coefs <- compData$coefEsts
rm(referenceMset)

# get counts
mSet <- mSet.subj
counts <- minfi:::projectCellType(getBeta(mSet)[rownames(coefs),], coefs)
rownames(counts) <- colnames(rgSet)




