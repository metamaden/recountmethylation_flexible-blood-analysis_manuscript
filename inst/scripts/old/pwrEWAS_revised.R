#!/usr/bin/env R

# Author: Sean Maden
# Revised functions for pwrEWAS

require(foreach)

#--------
# helpers
#--------

#' getAlphBet
#'
#'
#'
getAlphBet <- function(myMean, myVar){
  alpha <- myMean^2 * ((1-myMean)/myVar - 1/myMean) 
  beta <- alpha * (1/myMean - 1)
  return( list(alpha = alpha, beta = beta) )
}

#' getMeanVar
#'
#'
#'
getMeanVar <- function(myAlpha, myBeta){
  mean <- myAlpha / (myAlpha + myBeta)
  var <- myAlpha * myBeta / ((myAlpha + myBeta)^2 * (myAlpha + myBeta + 1))
  return( list(mean = mean, var = var) )
}

#' beta-value to m-value
#' 
#' 
#' 
beta2Mvalue <- function(beta){ # beta to m-value
  return( log2(beta / (1-beta)) )
}

#' t-test uses UNequal variance
#' 
#' 
#' 
ttestSlow <- function(g1Beta,g2Beta,rCnt,rTx,paired){
  mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
  ttest <- apply(mvals,1,  function (x) 
    stats::t.test(x[seq_len(rCnt)], x[(rCnt+1):(rCnt+rTx)], 
                  paired = FALSE, 
                  var.equal = TRUE))
  temp <- NULL
  temp$pval <- unlist(lapply(ttest, function(x) x$p.value))
  temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
  return(temp)
}


#' t-test uses equal variance
#' 
#' 
#' 
ttestFast <- function(g1Beta,g2Beta,rCnt,rTx){
  mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
  ttest <- genefilter::rowttests(mvals, 
                                 fac = factor(c(rep("g1",rCnt),rep("g2",rTx))), 
                                 tstatOnly = FALSE)  # faster: tstatOnly = T
  temp <- NULL
  temp$pval <- ttest$p.value
  temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
  return(temp)
}


#' Wilcox rank sum test
#' 
#' 
#' 
Wilcox <- function(g1Beta, g2Beta, rCnt,rTx){
  mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
  WRS <- apply(mvals,1, 
               function (x) stats::wilcox.test(x[seq_len(rCnt)] - 
                                                 x[(rCnt+1):(rCnt+rTx)], correct=TRUE))
  temp <- NULL
  temp$pval <- unlist(lapply(WRS, function(x) x$p.value))
  temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
  return(temp)
}

#' limma
#' 
#' 
#' 
limma <- function(g1Beta, g2Beta, rCnt,rTx){
  mvals <- cbind(beta2Mvalue(g1Beta), beta2Mvalue(g2Beta))
  design <- stats::model.matrix(~ c(rep("g1",rCnt),rep("g2",rTx)))
  limmaFit <- limma::lmFit(mvals, design)
  temp <- NULL
  temp$pval <- eBayes(limmaFit)$p.value[,2]
  temp$fdr <- stats::p.adjust(temp$pval, method = "fdr")
  return(temp)
}

#' CPGassoc
#' 
#' 
#' 
CPGassoc <- function(g1Beta, g2Beta, rCnt,rTx){
  assoc <- CpGassoc::cpg.assoc(cbind(g1Beta, g2Beta), 
                               c(rep("g1",rCnt),rep("g2",rTx)))
  temp <- NULL
  temp$pval <- assoc$results$P.value
  temp$fdr <- assoc$results$FDR
  return(temp)
}

#'
#'
#'
#'
getTau <- function(targetDmCpGs, targetDelta, methPara, 
                   detectionLimit, J, CpGonArray, maxCnt = 100){
  out <- NULL; tau <- 1; tauSteps <- 1; lookForTau <- TRUE;cnt<-0;
  while(cnt < maxCnt & lookForTau){
    percentile <- NULL
    for(i in seq_len(100)){ # simulate deltas for J CpG's (num sim CpG's later)
      cpgIdx4Tau <- sample(x = seq_len(CpGonArray), size = J, replace = TRUE) 
      delta <- truncnorm::rtruncnorm(1, mean = 0, sd = tau, 
                                     a=0.5 - methPara$mu[cpgIdx4Tau] - 
                                       sqrt(0.25-methPara$var[cpgIdx4Tau]), 
                                     b=0.5 - methPara$mu[cpgIdx4Tau] + 
                                       sqrt(0.25-methPara$var[cpgIdx4Tau]))
      percentile[i] <- stats::quantile(abs(delta),0.9999, na.rm = TRUE)}
    if(mean(percentile) < targetDelta - 0.5*detectionLimit & tau >= 1){
      tau <- tau + 1
    } else if(mean(percentile) < targetDelta - 0.5*detectionLimit & tau < 1){
      tauSteps <- 0.5 * tauSteps; tau <- tau + tauSteps
    } else if(mean(percentile) > targetDelta + 0.5*detectionLimit){
      tauSteps <- 0.5 * tauSteps; tau <- tau - tauSteps
    } else {lookForTau <- FALSE};cnt <- cnt + 1
    if(cnt == maxCnt){stop("Max iterations reached with finding Tau. ",
                           "Try increasing value of `maxCnt` argument")}
    }
  truelyDMperc <- mean(abs(delta) > detectionLimit, na.rm = TRUE)
  out$tau <- tau; targetK <- round(1/truelyDMperc * targetDmCpGs)
  out$K <- ifelse(targetK > J, J, targetK); return(out)
}

#'
#'
#'
#'
getK <- function(targetDmCpGs, methPara, detectionLimit, J, CpGonArray, tau){
  cpgIdx4Tau <- sample(x = seq_len(CpGonArray), size = J, replace = TRUE) # pick J random CpG's to be changed in mean meth
  delta <- truncnorm::rtruncnorm(1, mean = 0, sd = tau, 
                                 a=0.5 - methPara$mu[cpgIdx4Tau] - sqrt(0.25-methPara$var[cpgIdx4Tau]), 
                                 b=0.5 - methPara$mu[cpgIdx4Tau] + sqrt(0.25-methPara$var[cpgIdx4Tau]))
  
  truelyDMperc <- mean(abs(delta) > detectionLimit)
  targetK <- round(1/truelyDMperc * targetDmCpGs)
  K <- ifelse(targetK > J, J, targetK)
  return(K)
}

#' combine_tau
#'
#' @param listA
#' @param listB
#' @return 
#' @export
combine_tau <- function(listA, listB){
  if (is.null(listA)){return(listB)}
    
  if (!methods::is(listA[["power"]], "array") & 
      !methods::is(listA[["power"]], "matrix")){
    listA[["power"]] <- matrix(listA[["power"]])}
    
  if (!methods::is(listB[["power"]], "array") & 
      !methods::is(listB[["power"]], "matrix")){
    listB[["power"]] <- matrix(listB[["power"]])} 
    
  if (!methods::is(listA[["metric"]]$marTypeI, "array") & 
      !methods::is(listA[["metric"]]$marTypeI, "matrix")){
    listA[["metric"]]$marTypeI <- matrix(listA[["metric"]]$marTypeI)} 
    
  if (!methods::is(listB[["metric"]]$marTypeI, "array") & 
      !methods::is(listB[["metric"]]$marTypeI, "matrix")){
    listB[["metric"]]$marTypeI <- matrix(listB[["metric"]]$marTypeI)} 
    
  if (!methods::is(listA[["metric"]]$classicalPower, "array") & 
      !methods::is(listA[["metric"]]$classicalPower, "matrix")){
    listA[["metric"]]$classicalPower <- matrix(
      listA[["metric"]]$classicalPower)} 
    
  if (!methods::is(listB[["metric"]]$classicalPower, "array") & 
      !methods::is(listB[["metric"]]$classicalPower, "matrix")){
    listB[["metric"]]$classicalPower <- matrix(
      listB[["metric"]]$classicalPower)} 
    
  if (!methods::is(listA[["metric"]]$FDR, "array") & 
      !methods::is(listA[["metric"]]$FDR, "matrix")){
    listA[["metric"]]$FDR <- matrix(listA[["metric"]]$FDR)} 
    
  if (!methods::is(listB[["metric"]]$FDR, "array") & 
      !methods::is(listB[["metric"]]$FDR, "matrix")){
    listB[["metric"]]$FDR <- matrix(listB[["metric"]]$FDR)} 
    
  if (!methods::is(listA[["metric"]]$FDC, "array") & 
      !methods::is(listA[["metric"]]$FDC, "matrix")){
    listA[["metric"]]$FDC <- matrix(listA[["metric"]]$FDC)} 
   
  if (!methods::is(listB[["metric"]]$FDC, "array") & 
      !methods::is(listB[["metric"]]$FDC, "matrix")){
    listB[["metric"]]$FDC <- matrix(listB[["metric"]]$FDC)} 
    
  if (!methods::is(listA[["metric"]]$probTP, "array") & 
      !methods::is(listA[["metric"]]$probTP, "matrix")){
    listA[["metric"]]$probTP <- matrix(listA[["metric"]]$probTP)} 
    
  if (!methods::is(listB[["metric"]]$probTP, "array") & 
      !methods::is(listB[["metric"]]$probTP, "matrix")){
    listB[["metric"]]$probTP <- matrix(listB[["metric"]]$probTP)} 
    
  if (!methods::is(listA[["delta"]], "list")){
    listA[["delta"]] <- list(listA[["delta"]])} 
    
  returnList <- list()
  returnList[["power"]] <- abind::abind(listA[["power"]], 
                                        listB[["power"]], along = 3)
  returnList[["delta"]] <- listA[["delta"]]
  returnList[["delta"]][[length(listA[["delta"]]) + 1]] <- listB[["delta"]]
  returnList[["metric"]]$marTypeI <- abind::abind(listA[["metric"]]$marTypeI, 
                                                  listB[["metric"]]$marTypeI, 
                                                  along = 3)
  returnList[["metric"]]$classicalPower <- abind::abind(listA[["metric"]]$classicalPower, 
                                                        listB[["metric"]]$classicalPower, 
                                                        along = 3)
  returnList[["metric"]]$FDR <- abind::abind(listA[["metric"]]$FDR, 
                                             listB[["metric"]]$FDR, along = 3)
  returnList[["metric"]]$FDC <- abind::abind(listA[["metric"]]$FDC, 
                                             listB[["metric"]]$FDC, along = 3)
  returnList[["metric"]]$probTP <- abind::abind(listA[["metric"]]$probTP, 
                                                listB[["metric"]]$probTP, 
                                                along = 3)
  return(returnList)
}


#' combine_totSampleSizes
#' 
#' @param listA
#' @param listB
#' @return 
#' @export
combine_totSampleSizes <- function(listA, listB) {
  if (is.null(listA)){return(listB)} 
  returnList <- list()
  returnList[["power"]] <- cbind(listA[["power"]], listB[["power"]])
  returnList[["delta"]] <- cbind(listA[["delta"]], listB[["delta"]])
  returnList[["metric"]]$marTypeI <- cbind(listA[["metric"]]$marTypeI, 
                                           listB[["metric"]]$marTypeI)
  returnList[["metric"]]$classicalPower <- cbind(
    listA[["metric"]]$classicalPower, listB[["metric"]]$classicalPower)
  returnList[["metric"]]$FDR <- cbind(listA[["metric"]]$FDR, 
                                      listB[["metric"]]$FDR)
  returnList[["metric"]]$FDC <- cbind(listA[["metric"]]$FDC, 
                                      listB[["metric"]]$FDC)
  returnList[["metric"]]$probTP <- cbind(listA[["metric"]]$probTP, 
                                         listB[["metric"]]$probTP)
  return(returnList)
}

#' Progress bar function
#'
#' @param pb
#' @param n
#' @return 
#' @export
progress <- function(n){return(utils::setTxtProgressBar(pb, n))} 



#' pwrEWAS
#'
#' Main power function. Wraps parallelization script for simulations. Starts
#' by loading the 2-column table of tissue data (beta-value means and 
#' variances).
#'
#' @param tissueType Path to stats table containing columns for Beta-value mean 
#' ("mu") and variance ("var").
#' @param minTotSampleSize
#' @param maxTotSampleSize
#' @param SampleSizeSteps
#' @param NcntPer
#' @param targetDelta
#' @param deltaSD
#' @param J
#' @param targetDmCpGs
#' @param detectionLimit
#' @param DMmethod
#' @param FDRcritVal
#' @param core
#' @param sims
#' @return
#' @export 
pwrEWAS_itable <- function (tissueType, minTotSampleSize = 10, maxTotSampleSize = 20, 
                            SampleSizeSteps = 10, NcntPer = 0.5, maxCnt.tau = 1000,
                            targetDelta = NULL, deltaSD = NULL, J = 1e+05, 
                            targetDmCpGs = 10, detectionLimit = 0.01, 
                            DMmethod = c("limma", "t-test (unequal var)",
                                         "t-test (equal var)", "Wilcox rank sum", "CPGassoc"),
                            FDRcritVal = 0.05, core = 1, sims = 50)
{
  DMmethod <- match.arg(DMmethod)
  if (!is.null(targetDelta) & !is.null(deltaSD)){
    stop("Please specify only one: 'targetDelta' or 'deltaSD'")}
  # methPara <- get(load(tissueType)) # pwrEWAS.data:::loadDataset(tissueType)
  methPara <- tissueType
  CpGonArray <- length(methPara$mu); output <- NULL
  totSampleSizes <- seq(minTotSampleSize, maxTotSampleSize, SampleSizeSteps)
  cl <- parallel::makeCluster(core, setup_strategy = "sequential")
  doSNOW::registerDoSNOW(cl)
  if(is.null(deltaSD)){
    cat(paste("[", Sys.time(), "] ", "Finding tau...", sep = ""))
    K <- NULL; tau <- NULL
    for(d in seq_along(targetDelta)){
      myTau <- getTau(targetDmCpGs, targetDelta[d], methPara, 
                      detectionLimit, J, CpGonArray, maxCnt = maxCnt.tau)
      tau[d] <- myTau$tau; K[d] <- myTau$K}
    cat(paste("done", " [", Sys.time(), "]\n", sep = ""))
    print(paste("The following taus were chosen: ", 
                paste(tau, collapse = ", "), sep = ""))
  } else{
    tau <- deltaSD; K <- NULL
    for(d in seq_along(tau)){
      K[d] <- getK(targetDmCpGs,methPara,detectionLimit,
                            J,CpGonArray,tau)}}
  startTime <- Sys.time()
  cat(paste("[", startTime, "] ", "Running simulation\n", sep = ""))
  iterations <- length(totSampleSizes) * length(tau)
  pb <- utils::txtProgressBar(min = 0, max = iterations, style = 3)
  opts <- list(progress = progress); Ntot <- NULL
  multiThreadOut <- foreach(d = seq_along(tau), .combine = combine_tau, 
                            .packages = c("truncnorm", "limma", 
                                          "CpGassoc", "genefilter"),
                            .export = c("getAlphBet","getMeanVar","beta2Mvalue",
                                        "limma", "ttestSlow", "ttestFast", 
                                        "Wilcox", "CPGassoc")) %:% 
    foreach(Ntot = totSampleSizes, .combine = combine_totSampleSizes, 
            .options.snow = opts) %dopar% {
              utils::setTxtProgressBar(pb = pb, value = (d - 1) * length(totSampleSizes) + 
                                         which(Ntot == totSampleSizes))
              Ncnt <- round(Ntot * NcntPer)
              Ntx <- Ntot - Ncnt; marPower <- NULL; deltaSim <- NULL; marTypeI <- NULL
              FDR <- NULL; classicalPower <- NULL; FDC <- NULL; probTP <- NULL
              for (sim in seq_len(sims)) {
                cpgIdx <- sample(x = seq_len(CpGonArray), size = J, replace = TRUE)
                cpgIdxName <- paste(seq_len(J), "_", rownames(methPara)[cpgIdx], sep = "")
                changedCpgsIdx <- sample(x = cpgIdx, size = K[d])
                changedCpgsIdxName <- cpgIdxName[match(changedCpgsIdx, cpgIdx)]
                # handle delta maths with possible NaN's
                mean.a <- methPara$mu[changedCpgsIdx]
                sqrt.a <- sqrt(0.25 - methPara$var[changedCpgsIdx])
                which.na.a <- sqrt.a == "NaN"
                mean.b <- methPara$mu[changedCpgsIdx]
                sqrt.b <- sqrt(0.25 - methPara$var[changedCpgsIdx])
                which.na.b <- sqrt.b == "NaN"
                delta.a <- 0.5 - mean.a[!which.na.a] - sqrt.a[!which.na.a]
                delta.b <- 0.5 - mean.b[!which.na.b] + sqrt.b[!which.na.b]
                delta <- truncnorm::rtruncnorm(1, mean = 0, sd = as.numeric(tau[d]), 
                                               a = delta.a, b = delta.b)
                deltaSim <- c(deltaSim, delta)
                meaningfulDM <- (abs(delta) >= detectionLimit)
                meaningfulDMName <- changedCpgsIdxName[meaningfulDM]
                muToBeSimuUNchanged <- methPara$mu[cpgIdx]
                muToBeSimuChanged <- methPara$mu[cpgIdx]
                muToBeSimuChanged[match(changedCpgsIdx, cpgIdx)] <- muToBeSimuChanged[
                  match(changedCpgsIdx, cpgIdx)] + delta
                params_unchanged <- getAlphBet(myMean = muToBeSimuUNchanged, 
                                               myVar = methPara$var[cpgIdx])
                alpha_unchanged <- params_unchanged$alpha
                beta_unchanged <- params_unchanged$beta
                params_changed <- getAlphBet(myMean = muToBeSimuChanged, 
                                             myVar = methPara$var[cpgIdx])
                alpha_changed <- params_changed$alpha; beta_changed <- params_changed$beta
                g1Beta <- NULL; g2Beta <- NULL
                g1Beta <- matrix(stats::rbeta(J * Ncnt, 
                                              rep(alpha_unchanged, each = Ncnt), 
                                              rep(beta_unchanged, each = Ncnt)), 
                                 ncol = Ncnt, byrow = TRUE)
                g2Beta <- matrix(stats::rbeta(J * Ntx, 
                                              rep(alpha_changed, each = Ntx), 
                                              rep(beta_changed, each = Ntx)), 
                                 ncol = Ntx, byrow = TRUE)
                g1Beta[g1Beta == 1] <- max(g1Beta[g1Beta != 1])
                g2Beta[g2Beta == 1] <- max(g2Beta[g2Beta != 1])
                g1Beta[g1Beta == 0] <- min(g1Beta[g1Beta != 0])
                g2Beta[g2Beta == 0] <- min(g2Beta[g2Beta != 0])
                #rownames(g1Beta) <- rownames(g2Beta) <- 
                #  paste(seq_len(J), "_", names(alpha_unchanged), sep = "")
                rownames(g1Beta) <- rownames(g2Beta) <- 
                  paste(seq_len(J), "_", as.character(alpha_unchanged), sep = "")
                DMtest <- switch(DMmethod,
                                 `t-test (unequal var)` = ttestSlow(g1Beta, g2Beta, 
                                                                    Ncnt, Ntx, 
                                                                    paired = FALSE), 
                                 `t-test (equal var)` = ttestFast(g1Beta, g2Beta, 
                                                                  Ncnt, Ntx), 
                                 CPGassoc = CPGassoc(g1Beta, g2Beta, Ncnt, Ntx), 
                                 `Wilcox rank sum` = Wilcox(g1Beta, g2Beta, 
                                                            Ncnt, Ntx), 
                                 limma = limma(g1Beta, g2Beta, Ncnt, Ntx), 
                                 stop("Test not found"))
                notDM <- cpgIdxName[!(cpgIdxName %in% changedCpgsIdxName)]
                DM_negligible <- changedCpgsIdxName[!(changedCpgsIdxName %in% 
                                                        meaningfulDMName)]
                DM_meaningful <- changedCpgsIdxName[changedCpgsIdxName %in% 
                                                      meaningfulDMName]
                FP <- intersect(cpgIdxName[DMtest$fdr < FDRcritVal], notDM)
                NP <- intersect(cpgIdxName[DMtest$fdr < FDRcritVal], DM_negligible)
                TP <- intersect(cpgIdxName[DMtest$fdr < FDRcritVal], DM_meaningful)
                detectedCpGs <- cpgIdxName[DMtest$fdr < FDRcritVal]
                TN <- intersect(cpgIdxName[!(DMtest$fdr < FDRcritVal)], notDM)
                NN <- intersect(cpgIdxName[!(DMtest$fdr < FDRcritVal)], DM_negligible)
                FN <- intersect(cpgIdxName[!(DMtest$fdr < FDRcritVal)], DM_meaningful)
                marPower[sim] <- ifelse(length(DM_meaningful) > 0, 
                                        length(TP)/length(DM_meaningful), NA)
                marTypeI[sim] <- ifelse(length(notDM) > 0, length(FP)/length(notDM), NA)
                FDR[sim] <- ifelse(length(detectedCpGs) > 0, 
                                   length(FP)/length(detectedCpGs), NA)
                FDC[sim] <- ifelse(length(TP) > 0, 
                                   (length(FP))/length(TP), NA)
                classicalPower[sim] <- (length(NP) + length(TP))/(length(DM_negligible) + 
                                                                    length(DM_meaningful))
                probTP[sim] <- ifelse(length(TP) > 0, 1, 0)}
              outSim <- list(); outSim[["power"]] <- marPower
              outSim[["delta"]] <- deltaSim; outSim[["metric"]]$marTypeI <- marTypeI
              outSim[["metric"]]$FDR <- FDR
              outSim[["metric"]]$classicalPower <- classicalPower; outSim[["metric"]]$FDC <- FDC
              outSim[["metric"]]$probTP <- probTP; outSim}
  close(pb); parallel::stopCluster(cl)
  cat(paste("[", startTime, "] Running simulation ... done [",
            Sys.time(), "]\n", sep = ""))
  if (is.null(targetDelta)){targetDelta <- tau}
  if (length(targetDelta) == 1 & length(totSampleSizes) == 1){
    output$meanPower <- matrix(mean(multiThreadOut[["power"]], 
                                    na.rm = TRUE))} 
  if (length(targetDelta) == 1 & length(totSampleSizes) > 1){
    output$meanPower <- matrix(apply(multiThreadOut[["power"]], 
                                     2, mean, na.rm = TRUE), ncol = 1)}
  if (length(targetDelta) > 1){
    output$meanPower <- apply(multiThreadOut[["power"]], 
                              c(2, 3), mean, na.rm = TRUE)}
  rownames(output$meanPower) <- totSampleSizes
  colnames(output$meanPower) <- targetDelta
  if (length(targetDelta) == 1){
    output$powerArray <- array(data = multiThreadOut[["power"]], 
                               dim = c(sims, length(totSampleSizes), 
                                       length(targetDelta)))}
  if (length(targetDelta) > 1 & length(totSampleSizes) == 1){
    output$powerArray <- multiThreadOut[["power"]]}
  if (length(targetDelta) > 1){output$powerArray <- multiThreadOut[["power"]]}
  dimnames(output$powerArray) <- list(seq_len(sims),totSampleSizes,targetDelta)
  if (length(targetDelta) == 1 & length(totSampleSizes) == 1) {
    output$deltaArray <- list(matrix(multiThreadOut[["delta"]]))}
  if (length(targetDelta) == 1 & length(totSampleSizes) > 1){
    output$deltaArray <- list(multiThreadOut[["delta"]])}
  if (length(targetDelta) > 1 & length(totSampleSizes) == 1){
    output$deltaArray <- lapply(multiThreadOut[["delta"]], as.matrix)}
  if (length(targetDelta) > 1 & length(totSampleSizes) > 1){
    output$deltaArray <- multiThreadOut[["delta"]]} 
  names(output$deltaArray) <- targetDelta
  for (d in seq_along(targetDelta)) {
    colnames(output$deltaArray[[d]]) <- totSampleSizes}
  if (length(targetDelta) == 1 & length(totSampleSizes) == 1){
    output$metric$marTypeI <- matrix(mean(multiThreadOut[["metric"]]$marTypeI, 
                                          na.rm = TRUE))}
  if (length(targetDelta) == 1 & length(totSampleSizes) > 1){
    output$metric$marTypeI <- matrix(apply(multiThreadOut[["metric"]]$marTypeI, 
                                           2, mean, na.rm = TRUE), ncol = 1)}
  if (length(targetDelta) > 1){
    output$metric$marTypeI <- apply(multiThreadOut[["metric"]]$marTypeI, 
                                    c(2, 3), mean, na.rm = TRUE)}
   
  rownames(output$metric$marTypeI) <- totSampleSizes
  colnames(output$metric$marTypeI) <- targetDelta
  if (length(targetDelta) == 1 & length(totSampleSizes) == 1){
    output$metric$classicalPower <- matrix(
      mean(multiThreadOut[["metric"]]$classicalPower, na.rm = TRUE))} 
  if (length(targetDelta) == 1 & length(totSampleSizes) > 1){
    output$metric$classicalPower <- matrix(apply(
      multiThreadOut[["metric"]]$classicalPower, 2, mean, na.rm = TRUE), ncol = 1)}
    
  if (length(targetDelta) > 1){
    output$metric$classicalPower <- apply(multiThreadOut[["metric"]]$classicalPower, 
                                          c(2, 3), mean, na.rm = TRUE)}
    
  rownames(output$metric$classicalPower) <- totSampleSizes
  colnames(output$metric$classicalPower) <- targetDelta
  if (length(targetDelta) == 1 & length(totSampleSizes) == 1){
    output$metric$FDR <- matrix(mean(multiThreadOut[["metric"]]$FDR, 
                                     na.rm = TRUE))} 
  if (length(targetDelta) == 1 & length(totSampleSizes) > 1){
    output$metric$FDR <- matrix(apply(multiThreadOut[["metric"]]$FDR, 
                                      2, mean, na.rm = TRUE), ncol = 1)} 
  if (length(targetDelta) > 1){
    output$metric$FDR <- apply(multiThreadOut[["metric"]]$FDR, 
                               c(2, 3), mean, na.rm = TRUE)} 
    
  rownames(output$metric$FDR) <- totSampleSizes
  colnames(output$metric$FDR) <- targetDelta
  if (length(targetDelta) == 1 & length(totSampleSizes) == 1){
    output$metric$FDC <- matrix(mean(multiThreadOut[["metric"]]$FDC, 
                                     na.rm = TRUE))} 
  if (length(targetDelta) == 1 & length(totSampleSizes) > 1){
    output$metric$FDC <- matrix(apply(multiThreadOut[["metric"]]$FDC, 
                                      2, mean, na.rm = TRUE), ncol = 1)} 
    
  if (length(targetDelta) > 1){
    output$metric$FDC <- apply(multiThreadOut[["metric"]]$FDC, 
                               c(2, 3), mean, na.rm = TRUE)} 
  rownames(output$metric$FDC) <- totSampleSizes
  colnames(output$metric$FDC) <- targetDelta
  if (length(targetDelta) == 1 & length(totSampleSizes) == 1){
    output$metric$probTP <- matrix(mean(multiThreadOut[["metric"]]$probTP, 
                                        na.rm = TRUE))} 
  if (length(targetDelta) == 1 & length(totSampleSizes) > 1){
    output$metric$probTP <- matrix(apply(multiThreadOut[["metric"]]$probTP, 
                                         2, mean, na.rm = TRUE), ncol = 1)}
  if (length(targetDelta) > 1){
    output$metric$probTP <- apply(multiThreadOut[["metric"]]$probTP, 
                                  c(2, 3), mean, na.rm = TRUE)}
  rownames(output$metric$probTP) <- totSampleSizes
  colnames(output$metric$probTP) <- targetDelta; return(output)
}
