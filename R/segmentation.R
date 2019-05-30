# file:   segmentation.R
# author: Gavin Ha, Ph.D.
#         Justin Rhoades
#               Dana-Farber Cancer Institute
#               Broad Institute
# contact: <gavinha@broadinstitute.org>
# ULP-WGS website: http://www.broadinstitute.org/~gavinha/ULP-WGS/
# HMMcopy website: http://compbio.bccrc.ca/software/hmmcopy/ and https://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html
# date:   Oct 26, 2016
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

HMMsegment <- function(x, validInd = NULL, dataType = "copy", param = NULL, 
    chrTrain = c(1:22), maxiter = 50, estimateNormal = TRUE, estimatePloidy = TRUE, 
    estimatePrecision = TRUE, estimateSubclone = TRUE, estimateTransition = TRUE,
    estimateInitDist = TRUE, logTransform = FALSE, verbose = TRUE) {
  chr <- as.factor(seqnames(x[[1]]))
	# setup columns for multiple samples #
	dataMat <- as.matrix(as.data.frame(lapply(x, function(y) { mcols(x[[1]])[, dataType] })))
	
	# normalize by median and log data #
	if (logTransform){
    dataMat <- apply(dataMat, 2, function(x){ log(x / median(x, na.rm = TRUE)) })
	}else{
	  dataMat <- log(2^dataMat)
	}
	## update variable x with loge instead of log2
  for (i in 1:length(x)){
    mcols(x[[i]])[, dataType] <- dataMat[, i]
  }
  if (!is.null(chrTrain)) {
		chrInd <- chr %in% chrTrain
  }else{
  	chrInd <- !logical(length(chr))
  }
  if (!is.null(validInd)){
    chrInd <- chrInd & validInd
  }  

	if (is.null(param)){
		param <- getDefaultParameters(dataMat[chrInd])
	}
	#if (param$n_0 == 0){
	#	param$n_0 <- .Machine$double.eps
	#}
	####### RUN EM ##########
  convergedParams <- runEM(dataMat, chr, chrInd, param, maxiter, 
      verbose, estimateNormal = estimateNormal, estimatePloidy = estimatePloidy, 
      estimateSubclone = estimateSubclone, estimatePrecision = estimatePrecision, 
      estimateTransition = estimateTransition, estimateInitDist = estimateInitDist)
  # Calculate likelihood using converged params
 # S <- param$numberSamples
 # K <- length(param$ct)
 # KS <- K ^ S
 # py <- matrix(0, KS, nrow(dataMat))
 # iter <- convergedParams$iter
  # lambdasKS <- as.matrix(expand.grid(as.data.frame(convergedParams$lambda[, , iter])))
  # for (ks in 1:KS) {
  #   probs <- tdistPDF(dataMat, convergedParams$mus[ks, , iter], lambdasKS[ks, ], param$nu)
  #   py[ks, ] <- apply(probs, 1, prod) # multiply across samples for each data point to get joint likelihood.
  # }
  # 
  viterbiResults <- runViterbi(convergedParams, chr)
  
  # setup columns for multiple samples #
  segs <- segmentData(x, validInd, viterbiResults$states, convergedParams)
  #output$segs <- processSegments(output$segs, chr, start(x), end(x), x$DataToUse)
  names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  #if (c(0) %in% param$ct){ #if state 0 HOMD is IN params#
  	#names <- c("HOMD", names)
  	# shift states to start at 2 (HETD)
    #tmp <- lapply(segs, function(x){ x$state <- x$state + 1; x})
    #viterbiResults$states <- as.numeric(viterbiResults$states) + 1
	#}
	### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
  cnaList <- list()
  S <- length(x)
  for (s in 1:S){
    id <- names(x)[s]
    copyNumber <- param$jointCNstates[viterbiResults$state, s]
    subclone.status <- param$jointSCstatus[viterbiResults$state, s]
  	cnaList[[id]] <- data.frame(cbind(sample = as.character(id), 
                  chr = as.character(seqnames(x[[s]])),	
                  start = start(x[[s]]), end = end(x[[s]]), 
                  copy.number = copyNumber,
                  event = names[copyNumber + 1], 
                  logR = round(log2(exp(dataMat[,s])), digits = 4),
                  subclone.status = as.numeric(subclone.status)
  	))
  
    cnaList[[id]] <- transform(cnaList[[id]], 
                              start = as.integer(as.character(start)),
                              end = as.integer(as.character(end)), 
                              copy.number = as.numeric(copy.number),
                              logR = as.numeric(as.character(logR)),
                              subclone.status = as.numeric(subclone.status))
  
  	## order by chromosome ##
  	chrOrder <- unique(chr) #c(1:22,"X","Y")
  	cnaList[[id]] <- cnaList[[id]][order(match(cnaList[[id]][, "chr"],chrOrder)),]
  	## remove MT chr ##
    cnaList[[id]] <- cnaList[[id]][cnaList[[id]][,"chr"] %in% chrOrder, ]
    
    ## segment mean loge -> log2
    #segs[[s]]$median.logR <- log2(exp(segs[[s]]$median.logR))
    segs[[s]]$median <- log2(exp(segs[[s]]$median))
    ## add subclone status
    segs[[s]]$subclone.status <-  param$jointSCstatus[segs[[s]]$state, s]
  }	
  convergedParams$segs <- segs
  return(list(cna = cnaList, results = convergedParams, viterbiResults = viterbiResults))
}

getTransitionMatrix <- function(K, e, strength){
  A <- matrix(0, K, K)
  for (j in 1:K) {
    A[j, ] <- (1 - e[1]) / (K - 1)
    A[j, j] <- e[1]
  }
  A <- normalize(A)
  A_prior <- A
  dirPrior <- A * strength[1]
  return(list(A=A, dirPrior=dirPrior))
}

getDefaultParameters <- function(x, maxCN = 5, ct.sc = NULL, ploidy = 2, e = 0.9999999, e.sameState = 10, strength = 10000000, includeHOMD = FALSE){
  if (includeHOMD){
    ct <- 0:maxCN
  }else{
    ct <- 1:maxCN
  }
	param <- list(
		strength = strength, e = e,
		ct = c(ct, ct.sc),
		ct.sc.status = c(rep(FALSE, length(ct)), rep(TRUE, length(ct.sc))),
		phi_0 = 2, alphaPhi = 4, betaPhi = 1.5,
		n_0 = 0.5, alphaN = 2, betaN = 2,
		sp_0 = 0.5, alphaSp = 2, betaSp = 2,
		lambda = as.matrix(rep(100, length(ct)+length(ct.sc)), ncol=1),
		nu = 2.1,
		kappa = rep(75, length(ct)), 
		alphaLambda = 5
	)
	K <- length(param$ct)
  ## initialize hyperparameters for precision using observed data ##
	if (!is.null(dim(x))){ # multiple samples (columns)
    param$numberSamples <- ncol(x)
    #betaLambdaVal <- ((apply(x, 2, function(x){ sd(diff(x), na.rm=TRUE) }) / sqrt(length(param$ct))) ^ 2)
    betaLambdaVal <- ((apply(x, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)   
	}else{ # only 1 sample
	  param$numberSamples <- 1
	  betaLambdaVal <- ((sd(x, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
	}
	param$betaLambda <- matrix(betaLambdaVal, ncol = param$numberSamples, nrow = length(param$ct), byrow = TRUE)
  param$alphaLambda <- rep(param$alphaLambda, K)
  
	# increase prior precision for -1, 0, 1 copies at ploidy
	#param$lambda[param$ct %in% c(1,2,3)] <- 1000 # HETD, NEUT, GAIN
	#param$lambda[param$ct == 4] <- 100 
	#param$lambda[which.max(param$ct)] <- 50 #highest CN
	#param$lambda[param$ct == 0] <- 1 #HOMD
	S <- param$numberSamples
	logR.var <- 1 / ((apply(x, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
	if (!is.null(dim(x))){ # multiple samples (columns)
		param$lambda <- matrix(logR.var, nrow=K, ncol=S, byrow=T, dimnames=list(c(),colnames(x)))
	}else{ # only 1 sample    
		#logR.var <- 1 / ((sd(x, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
    param$lambda <- matrix(logR.var, length(param$ct))
    param$lambda[param$ct %in% c(2)] <- logR.var 
    param$lambda[param$ct %in% c(1,3)] <- logR.var 
    param$lambda[param$ct >= 4] <- logR.var / 5
    param$lambda[param$ct == max(param$ct)] <- logR.var / 15
    param$lambda[param$ct.sc.status] <- logR.var / 10
  }
  # define joint copy number states #
  param$jointCNstates <- expand.grid(rep(list(param$ct), S))
  param$jointSCstatus <- expand.grid(rep(list(param$ct.sc.status), S))
  colnames(param$jointCNstates) <- paste0("Sample.", 1:param$numberSamples)
  colnames(param$jointSCstatus) <- paste0("Sample.", 1:param$numberSamples)
  
	# Initialize transition matrix to the prior
	txn <- getTransitionMatrix(K ^ S, e, strength)
  ## set higher transition probs for same CN states across samples ##
  # joint states where at least "tol" fraction of samples with the same CN state
	#apply(param$jointCNstates, 1, function(x){ sum(duplicated(as.numeric(x))) > 0 })
  cnStateDiff <- apply(param$jointCNstates, 1, function(x){ (abs(max(x) - min(x)))})
  if (e.sameState > 0 & S > 1){
		txn$A[, cnStateDiff == 0] <- txn$A[, cnStateDiff == 0] * e.sameState * K 
		txn$A[, cnStateDiff >= 3] <- txn$A[, cnStateDiff >=3]  / e.sameState / K
	}
  for (i in 1:nrow(txn$A)){
    for (j in 1:ncol(txn$A)){
      if (i == j){
        txn$A[i, j] <- e
      }
    }
  }
  txn$A <- normalize(txn$A)
	param$A <- txn$A
	param$dirPrior <- txn$A * strength[1] 
  param$A[, param$ct.sc.status] <- param$A[, param$ct.sc.status] / 10
  param$A <- normalize(param$A)
  param$dirPrior[, param$ct.sc.status] <- param$dirPrior[, param$ct.sc.status] / 10
  
  if (includeHOMD){
    K <- length(param$ct)
    param$A[1, 2:K] <- param$A[1, 2:K] * 1e-5; param$A[2:K, 1] <- param$A[2:K, 1] * 1e-5;
    param$A[1, 1] <- param$A[1, 1] * 1e-5
    param$A <- normalize(param$A); param$dirPrior <- param$A * param$strength
  }

  param$kappa <- rep(75, K ^ S)
  param$kappa[cnStateDiff == 0] <- param$kappa[cnStateDiff == 0] + 125
	param$kappa[cnStateDiff >=3] <- param$kappa[cnStateDiff >=3] - 50
	param$kappa[which(rowSums(param$jointCNstates==2) == S)] <- 800
  
  return(param)
}



segmentData <- function(dataGR, validInd, states, convergedParams){
  if (sum(convergedParams$param$ct == 0) ==0){
  	includeHOMD <- FALSE
  }else{
  	includeHOMD <- TRUE
  }
  if (!includeHOMD){
    names <- c("HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  }else{
    names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  }
  states <- states[validInd]
  S <- length(dataGR)
  jointStates <- convergedParams$param$jointCNstates
  jointSCstatus <- convergedParams$param$jointSCstatus
  colNames <- c("seqnames", "start", "end", "copy")
  segList <- list()
  for (i in 1:S){
  	id <- names(dataGR)[i]
  	dataIn <- dataGR[[i]][validInd, ]
    rleResults <- t(sapply(runValue(seqnames(dataIn)), function(x){
      ind <- as.character(seqnames(dataIn)) == x
      r <- rle(states[ind])
    }))
    rleLengths <- unlist(rleResults[, "lengths"])
    rleValues <- unlist(rleResults[, "values"])
    sampleDF <- as.data.frame(dataIn)
    numSegs <- length(rleLengths)
    segs <- as.data.frame(matrix(NA, ncol = 7, nrow = numSegs, 
                   dimnames = list(c(), c("chr", "start", "end", "state", "event", "median", "copy.number"))))
    prevInd <- 0
    for (j in 1:numSegs){
      start <- prevInd + 1
      end <- prevInd + rleLengths[j]
      segDF <- sampleDF[start:end, colNames]
      prevInd <- end
      numR <- nrow(segDF)
      segs[j, "chr"] <- as.character(segDF[1, "seqnames"])
      segs[j, "start"] <- segDF[1, "start"]
      segs[j, "state"] <- rleValues[j]
      segs[j, "copy.number"] <- jointStates[rleValues[j], i]
      if (segDF[1, "seqnames"] == segDF[numR, "seqnames"]){
        segs[j, "end"] <- segDF[numR, "end"]
        segs[j, "median"] <- round(median(segDF$copy, na.rm = TRUE), digits = 6)
        if (includeHOMD){
        	segs[j, "event"] <- names[segs[j, "copy.number"] + 1]
        }else{
        	segs[j, "event"] <- names[segs[j, "copy.number"]]
        }
      }else{ # segDF contains 2 different chromosomes
        print(j)
      }                                      
    }
    segList[[id]] <- segs
  }
  return(segList)
}
    

runViterbi <- function(convergedParams, chr){
  message("runViterbi: Segmenting and classifying")
  chrs <- levels(chr)
  chrsI <- vector('list', length(chrs))
  # initialise the chromosome index and the init state distributions
  for(i in 1:length(chrs)) {
    chrsI[i] <- list(which(chr == chrs[i]))
  }
  segs <- vector('list', length(chrs))
  py <- convergedParams$py
  N <- ncol(py)
  Z <- rep(0, N)
  convergeIter <- convergedParams$iter
  piG <- convergedParams$pi[, convergeIter]
  A <- convergedParams$A


  for(c in 1:length(chrsI)) {
    I <- chrsI[[c]]
    output <- .Call("viterbi", log(piG), log(A), log(py[, I]), PACKAGE = "HMMcopy")
    Z[I] <- output$path
    segs[[c]] <- output$seg
  }
  return(list(segs=segs, states=Z))
}

# Normalize a given array to sum to 1
normalize <- function(A) {
	vectorNormalize <- function(x){ x / (sum(x) + (sum(x) == 0)) }
	if (length(dim(A)) < 2){
  	M <- vectorNormalize(A)
  }else{
  	M <- t(apply(A, 1, vectorNormalize))
  }
  return(M);
}


# processSegments <- function(seg, chr, start, end, copy) {
#   segment <- data.frame()
#   chromosomes <- levels(chr)
#   for (i in 1:length(chromosomes)) {
#     seg_length = dim(seg[[i]])[1]
#     chr_name <- rep(chromosomes[i], seg_length)
#     chr_index <- which(chr == chromosomes[i])
#     chr_start <- start[chr_index][seg[[i]][, 1]]
#     chr_stop <- end[chr_index][seg[[i]][, 2]]
#     chr_state <- seg[[i]][, 3]
#     chr_median <- rep(0, seg_length)
#     for(j in 1:seg_length) {
#       chr_median[j] <-
#         median(na.rm = TRUE, log2(exp(copy[chr_index][seg[[i]][j, 1]:seg[[i]][j, 2]])))
#     }
#     segment <- rbind(segment, cbind(chr = chr_name,
#       start = as.numeric(chr_start), end = chr_stop, state = chr_state,
#       median = chr_median))
#   }
#   segment <- transform(segment, start = as.numeric(as.character(start)),
#     end = as.numeric(as.character(end)), as.numeric(as.character(state)),
#     median = as.numeric(as.character(median)))
#   return(segment)
# }
