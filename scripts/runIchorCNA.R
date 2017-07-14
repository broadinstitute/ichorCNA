# file:   ichorCNA.R
# author: Gavin Ha, Ph.D.
#               Dana-Farber Cancer Institute
#               Broad Institute
# contact: <gavinha@broadinstitute.org>
#         Justin Rhoades
#               Broad Institute
# contact: <rhoades@broadinstitute.org>

# ichorCNA: https://github.com/broadinstitute/ichorCNA
# HMMcopy website: http://compbio.bccrc.ca/software/hmmcopy/ and https://www.bioconductor.org/packages/release/bioc/html/HMMcopy.html
# date:   March 29, 2017
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This script is the main script to run the HMM.

library(optparse)

option_list <- list(
  make_option(c("--WIG"), type = "character", default=NULL, help = "Path to tumor WIG file."),
  make_option(c("--NORMWIG"), type = "character", default=NULL, help = "Path to normal WIG file."),
  make_option(c("--gcWig"), type = "character", default=NULL, help = "Path to GC-content WIG file."),
  make_option(c("--mapWig"), type = "character", default=NULL, help = "Path to mappability score WIG file."),
  make_option(c("--normalPanel"), type="character", default=NULL, help="Median corrected depth from panel of normals"),
  make_option(c("--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions."),
  make_option(c("--id"), type = "character", help = "Patient ID."),
  #make_option(c("--gcOnly"), type="logical", default=FALSE, help = "TRUE if only correct for GC content bias"),
  make_option(c("--centromere"), type="character", default=NULL, help = "File containing Centromere locations"),
  make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove"),
  make_option(c("--normal"), type="character", default="0.5", help = "Initial normal contamination"),
  make_option(c("--scStates"), type="character", default="NULL", help = "Subclonal states to consider"),
  make_option(c("--coverage"), type="numeric", default=NULL, help = "PICARD sequencing coverage"),
  make_option(c("--lambda"), type="character", default="NULL", help="Initial Student's t precision; must contain 4 values (e.g. c(1500,1500,1500,1500))"),
  make_option(c("--lambdaScaleHyperParam"), type="numeric", default=3, help="Hyperparameter (scale) for Gamma prior on Student's-t precision."),
  #	make_option(c("--kappa"), type="character", default=50, help="Initial state distribution"),
  make_option(c("--ploidy"), type="character", default="2", help = "Initial tumour ploidy"),
  make_option(c("--maxCN"), type="numeric", default=7, help = "Total clonal CN states"),
  make_option(c("--estimateNormal"), type="logical", default=TRUE, help = "Estimate normal."),
  make_option(c("--estimateScPrevalence"), type="logical", default=TRUE, help = "Estimate subclonal prevalence."),
  make_option(c("--estimatePloidy"), type="logical", default=TRUE, help = "Estimate tumour ploidy."),
  make_option(c("--maxFracCNASubclone"), type="numeric", default=0.7, help="Exclude solutions with fraction of subclonal events greater than this value."),
  make_option(c("--maxFracGenomeSubclone"), type="numeric", default=0.5, help="Exclude solutions with subclonal genome fraction greater than this value."),
  make_option(c("--minSegmentBins"), type="numeric", default=50, help="Minimum number of bins for largest segment threshold required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
  make_option(c("--altFracThreshold"), type="numeric", default=0.05, help="Minimum proportion of bins altered required to estimate tumor fraction; if below this threshold, then will be assigned zero tumor fraction."),
  make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases"),
  make_option(c("--chrTrain"), type="character", default="c(1:22)", help = "Specify chromosomes to estimate params."),
  make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze."),
  make_option(c("--diploidChrX"), type="logical", default=FALSE, help = "Assume 2 copies of chrX"),
  make_option(c("--normalizeMaleX"), type="logical", default=FALSE, help = "If male, then normalize chrX by median"),
  make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male."),
  make_option(c("--includeHOMD"), type="logical", default=FALSE, help="If false, then exclude HOMD state"),
  make_option(c("--txnE"), type="numeric", default=0.9999999, help = "Self-transition probability"),
  make_option(c("--txnStrength"), type="numeric", default=10000000, help = "Transition pseudo-counts"),
  	make_option(c("--plotFileType"), type="character", default="pdf", help = "File format for output plots"),
	make_option(c("--plotYLim"), type="character", default="c(-2,2)", help = "ylim to use for chromosome plots"),
  make_option(c("--outDir"), type="character", default="./", help = "Output Directory"),
  make_option(c("--libdir"), type = "character", default=NULL, help = "Script library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)
options(scipen=0, stringsAsFactors=F)

library(HMMcopy)
library(ichorCNA)
options(stringsAsFactors=FALSE)
options(bitmapType='cairo')

patientID <- opt$id
tumour_file <- opt$WIG
normal_file <- opt$NORMWIG
gcWig <- opt$gcWig
mapWig <- opt$mapWig
normal_panel <- opt$normalPanel
gcOnly <- as.logical(opt$gcOnly)  # {TRUE, FALSE}
exons.bed <- opt$exons.bed  # "0" if none specified
centromere <- opt$centromere
flankLength <- opt$rmCentromereFlankLength
normal <- eval(parse(text = opt$normal))
scStates <- eval(parse(text = opt$scStates))
lambda <- eval(parse(text = opt$lambda))
lambdaScaleHyperParam <- opt$lambdaScaleHyperParam
estimateNormal <- opt$estimateNormal
estimatePloidy <- opt$estimatePloidy
estimateScPrevalence <- opt$estimateScPrevalence
maxFracCNASubclone <- opt$maxFracCNASubclone
maxFracGenomeSubclone <- opt$maxFracGenomeSubclone
minSegmentBins <- opt$minSegmentBins
altFracThreshold <- opt$altFracThreshold
ploidy <- eval(parse(text = opt$ploidy))
coverage <- opt$coverage
maxCN <- opt$maxCN
txnE <- opt$txnE
txnStrength <- opt$txnStrength
diploidChrX <- as.logical(opt$diploidChrX)
normalizeMaleX <- as.logical(opt$normalizeMaleX)
includeHOMD <- as.logical(opt$includeHOMD)
fracReadsInChrYForMale <- opt$fracReadsInChrYForMale
outDir <- opt$outDir
libdir <- opt$libdir
plotFileType <- opt$plotFileType
plotYLim <- eval(parse(text=opt$plotYLim))
gender <- NULL
outImage <- paste0(outDir,"/", patientID,".RData")
diploidChrXThreshold <- 0.5
chrs <- eval(parse(text = opt$chrs))
#chrs <- c(chrs, "Y")
chrTrain <- eval(parse(text=opt$chrTrain))
chrNormalize <- eval(parse(text=opt$chrNormalize))

if (!is.null(libdir)){
	source(paste0(libdir,"/utils.R"))
	source(paste0(libdir,"/segmentation.R"))
	source(paste0(libdir,"/EM.R"))
	source(paste0(libdir,"/output.R"))
	source(paste0(libdir,"/plotting.R"))
}
if (substr(tumour_file,nchar(tumour_file)-2,nchar(tumour_file)) == "wig") {
  wigFiles <- data.frame(cbind(patientID, tumour_file))
} else {
  wigFiles <- read.delim(tumour_file, header=F, as.is=T)
}
save.image()
## FILTER BY EXONS IF PROVIDED ##
## add gc and map to RangedData object ##
if (is.null(exons.bed) || exons.bed == "None" || exons.bed == "NULL"){
  targetedSequences <- NULL
}else{
  targetedSequences <- read.delim(exons.bed, header=T, sep="\t")  
}

## load PoN
if (is.null(normal_panel) || normal_panel == "None" || normal_panel == "NULL"){
	normal_panel <- NULL
}

if (is.null(centromere) || centromere == "None" || centromere == "NULL"){ # no centromere file provided
	centromere <- system.file("extdata", "GRCh37.p13_centromere_UCSC-gapTable.txt", 
			package = "ichorCNA")
}
centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
save.image(outImage)

## LOAD IN WIG FILES ##
numSamples <- nrow(wigFiles)
tumour_counts <- list()
tumour_copy <- list()
for (i in 1:numSamples) {
  id <- wigFiles[i,1]
  ## create output directories for each sample ##
  dir.create(paste0(outDir, "/", id, "/"), recursive = TRUE)
  ### LOAD TUMOUR AND NORMAL FILES ###
  message("Loading tumour file:", wigFiles[i,1])
  tumour_reads <- wigToRangedData(wigFiles[i,2])
  
  ## LOAD GC/MAP WIG FILES ###
	# find the bin size and load corresponding wig files #
	binSize <- as.data.frame(tumour_reads[1,])$width 
	if (is.null(gcWig) || gcWig == "None" || gcWig == "NULL"){
		gcWig <- system.file("extdata", paste0("gc_hg19_", binSize / 1000, "kb.wig"), 
				package = "ichorCNA")
	}
	if (is.null(mapWig) || mapWig == "None" || mapWig == "NULL"){
		mapWig <- system.file("extdata", paste0("map_hg19_", binSize / 1000, "kb.wig"), 
				package = "ichorCNA")
	}
	message("Reading GC and mappability files")
	gc <- wigToRangedData(gcWig)
	map <- wigToRangedData(mapWig)
  message("Correcting Tumour")
  
  counts <- loadReadCountsFromWig(tumour_reads, chrs = chrs, gc = gc, map = map, 
                                       centromere = centromere, flankLength = flankLength, 
                                       targetedSequences = targetedSequences, 
                                       chrNormalize = chrNormalize, mapScoreThres = 0.9)
  tumour_copy[[id]] <- counts$counts #as(counts$counts, "GRanges")
  gender <- counts$gender
 	## load in normal file if provided 
 	if (!is.null(normal_file)){
		message("Loading normal file:", normal_file)
		normal_reads <- wigToRangedData(normal_file)
		message("Correcting Normal")
		counts <- loadReadCountsFromWig(normal_reads, chrs=chrs, gc=gc, map=map, 
				centromere=centromere, flankLength = flankLength, targetedSequences=targetedSequences,
				chrNormalize = chrNormalize, mapScoreThres = 0.9)
		normal_copy <- counts$counts #as(counts$counts, "GRanges")
		gender.normal <- counts$gender
	}else{
	  normal_copy <- NULL
	}

	### DETERMINE GENDER ###
	## if normal file not given, use chrY, else use chrX
	message("Determining gender...", appendLF = FALSE)
	gender.mismatch <- FALSE
	if (!is.null(normal_copy)){
	  if (gender$gender != gender.normal$gender){ #use tumour # use normal if given
		# check if normal is same gender as tumour
		  gender.mismatch <- TRUE
		}
	}
	message("Gender ", gender$gender)

  ## NORMALIZE GENOME-WIDE BY MATCHED NORMAL OR NORMAL PANEL (MEDIAN) ##
  tumour_copy[[id]] <- normalizeByPanelOrMatchedNormal(tumour_copy[[id]], chrs = c(1:22, "X", "Y"), 
      normal_panel = normal_panel, normal_copy = normal_copy, 
      gender = gender$gender, normalizeMaleX = normalizeMaleX)
	
	### OUTPUT FILE ###
	### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
	outMat <- as.data.frame(tumour_copy[[id]])
	outMat <- outMat[,c(1,2,3,12)]
	colnames(outMat) <- c("chr","start","end","log2_TNratio_corrected")
	outFile <- paste0(outDir,"/",id,".correctedDepth.txt")
	message(paste("Outputting to:", outFile))
	write.table(outMat, file=outFile, row.names=F, col.names=T, quote=F, sep="\t")

} ## end of for each sample


chrInd <- space(tumour_copy[[1]]) %in% chrTrain
## get positions that are valid
valid <- tumour_copy[[1]]$valid
if (length(tumour_copy) >= 2) {
  for (i in 2:length(tumour_copy)){ 
    valid <- valid & tumour_copy[[i]]$valid 
  } 
}

### RUN HMM ###
## store the results for different normal and ploidy solutions ##
ptmTotalSolutions <- proc.time() # start total timer
results <- list()
loglik <- as.data.frame(matrix(NA, nrow = length(normal) * length(ploidy), ncol = 7, 
                 dimnames = list(c(), c("init", "n_est", "phi_est", "BIC", 
                 												"Frac_genome_subclonal", "Frac_CNA_subclonal", "loglik"))))
counter <- 1
compNames <- rep(NA, nrow(loglik))
#### restart for purity and ploidy values ####
for (n in normal){
  for (p in ploidy){
    logR <- as.data.frame(lapply(tumour_copy, "[[", "copy")) # NEED TO EXCLUDE CHR X #
    param <- getDefaultParameters(logR[valid & chrInd, , drop=F], maxCN = maxCN, includeHOMD = includeHOMD, 
                ct.sc=scStates, ploidy = floor(p), e=txnE, e.same = 50, strength=txnStrength)
    param$phi_0 <- rep(p, numSamples)
    param$n_0 <- rep(n, numSamples)
    
    ############################################
    ######## CUSTOM PARAMETER SETTINGS #########
    ############################################
    # 0.1x cfDNA #
    if (is.null(lambda)){
			logR.var <- 1 / ((apply(logR, 2, sd, na.rm = TRUE) / sqrt(length(param$ct))) ^ 2)
			param$lambda <- rep(logR.var, length(param$ct))
			param$lambda[param$ct %in% c(2)] <- logR.var 
			param$lambda[param$ct %in% c(1,3)] <- logR.var 
			param$lambda[param$ct >= 4] <- logR.var / 5
			param$lambda[param$ct == max(param$ct)] <- logR.var / 15
			param$lambda[param$ct.sc.status] <- logR.var / 10
    }else{
			param$lambda[param$ct %in% c(2)] <- lambda[2]
			param$lambda[param$ct %in% c(1)] <- lambda[1]
			param$lambda[param$ct %in% c(3)] <- lambda[3]
			param$lambda[param$ct >= 4] <- lambda[4]
			param$lambda[param$ct == max(param$ct)] <- lambda[2] / 15
			param$lambda[param$ct.sc.status] <- lambda[2] / 10
		}
		param$alphaLambda <- rep(lambdaScaleHyperParam, length(param$ct))  
    # 1x bulk tumors #
    #param$lambda[param$ct %in% c(2)] <- 2000
    #param$lambda[param$ct %in% c(1)] <- 1750
    #param$lambda[param$ct %in% c(3)] <- 1750
    #param$lambda[param$ct >= 4] <- 1500
    #param$lambda[param$ct == max(param$ct)] <- 1000 / 25
		#param$lambda[param$ct.sc.status] <- 1000 / 75
		#param$alphaLambda[param$ct.sc.status] <- 4
		#param$alphaLambda[param$ct %in% c(1,3)] <- 5
		#param$alphaLambda[param$ct %in% c(2)] <- 5
		#param$alphaLambda[param$ct == max(param$ct)] <- 4
				
		#############################################
		################ RUN HMM ####################
		#############################################
    hmmResults.cor <- HMMsegment(tumour_copy, valid, dataType = "copy", 
                                 param = param, chrTrain = chrTrain, maxiter = 50,
                                 estimateNormal = estimateNormal, estimatePloidy = estimatePloidy,
                                 estimateSubclone = estimateScPrevalence, verbose = TRUE)
    
    for (s in 1:numSamples){
			iter <- hmmResults.cor$results$iter
			id <- names(hmmResults.cor$cna)[s]
	
			## convert full diploid solution (of chrs to train) to have 1.0 normal or 0.0 purity
			## check if there is an altered segment that has at least a minimum # of bins
			segsS <- hmmResults.cor$results$segs[[s]]
			segsS <- segsS[segsS$chr %in% chrTrain, ]
			segAltInd <- which(segsS$event != "NEUT")
			maxBinLength = -Inf
			if (sum(segAltInd) > 0){
				maxInd <- which.max(segsS$end[segAltInd] - segsS$start[segAltInd] + 1)
				maxSegRD <- RangedData(space=segsS$chr[segAltInd[maxInd]], 
									ranges=IRanges(start=segsS$start[segAltInd[maxInd]], end=segsS$end[segAltInd[maxInd]]))
				hits <- findOverlaps(query=maxSegRD, subject=tumour_copy[[s]][valid, ])
				maxBinLength <- length(subjectHits(hits))
			}
			## check if there are proportion of total bins altered 
			# if segment size smaller than minSegmentBins, but altFrac > altFracThreshold, then still estimate TF
			cnaS <- hmmResults.cor$cna[[s]]
			altInd <- cnaS[cnaS$chr %in% chrTrain, "event"] == "NEUT"
			altFrac <- sum(!altInd, na.rm=TRUE) / length(altInd)
			if ((maxBinLength <= minSegmentBins) & (altFrac <= altFracThreshold)){
				hmmResults.cor$results$n[s, iter] <- 1.0
    	}
      ## plot solution ##
      outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide_", "n", n, "-p", p)
      plotGWSolution(hmmResults.cor, s=s, outPlotFile=outPlotFile, plotFileType=plotFileType, 
                     plotYLim=plotYLim, estimateScPrevalence=estimateScPrevalence, main=id)
    }
    iter <- hmmResults.cor$results$iter
    results[[counter]] <- hmmResults.cor
    loglik[counter, "loglik"] <- signif(hmmResults.cor$results$loglik[iter], digits = 4)
    subClonalBinCount <- unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$subclone.status) }))
    fracGenomeSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ nrow(x) }))
    fracAltSub <- subClonalBinCount / unlist(lapply(hmmResults.cor$cna, function(x){ sum(x$copy.number != 2) }))
    fracAltSub <- lapply(fracAltSub, function(x){if (is.na(x)){0}else{x}})
    loglik[counter, "Frac_genome_subclonal"] <- paste0(signif(fracGenomeSub, digits=2), collapse=",")
    loglik[counter, "Frac_CNA_subclonal"] <- paste0(signif(as.numeric(fracAltSub), digits=2), collapse=",")
    loglik[counter, "init"] <- paste0("n", n, "-p", p)
    loglik[counter, "n_est"] <- paste(signif(hmmResults.cor$results$n[, iter], digits = 2), collapse = ",")
    loglik[counter, "phi_est"] <- paste(signif(hmmResults.cor$results$phi[, iter], digits = 4), collapse = ",")

    counter <- counter + 1
  }
}
## get total time for all solutions ##
elapsedTimeSolutions <- proc.time() - ptmTotalSolutions
message("Total ULP-WGS HMM Runtime: ", format(elapsedTimeSolutions[3] / 60, digits = 2), " min.")

### SAVE R IMAGE ###
save.image(outImage)
#save(tumour_copy, results, loglik, file=paste0(outDir,"/",id,".RData"))

### SELECT SOLUTION WITH LARGEST LIKELIHOOD ###
if (estimateScPrevalence){ ## sort but excluding solutions with too large % subclonal 
	fracInd <- which(loglik[, "Frac_CNA_subclonal"] <= maxFracCNASubclone & 
						 		   loglik[, "Frac_genome_subclonal"] <= maxFracGenomeSubclone)
	if (length(fracInd) > 0){ ## if there is a solution satisfying % subclonal
		ind <- fracInd[order(loglik[fracInd, "loglik"], decreasing=TRUE)]
	}else{ # otherwise just take largest likelihood
		ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
	}
}else{#sort by likelihood only
  ind <- order(as.numeric(loglik[, "loglik"]), decreasing=TRUE) 
}

hmmResults.cor <- results[[ind[1]]]
hmmResults.cor$results$loglik <- as.data.frame(loglik)
hmmResults.cor$results$gender <- gender$gender
hmmResults.cor$results$chrYCov <- gender$chrYCovRatio
hmmResults.cor$results$chrXMedian <- gender$chrXMedian
hmmResults.cor$results$coverage <- coverage

outputHMM(cna = hmmResults.cor$cna, segs = hmmResults.cor$results$segs, 
                      results = hmmResults.cor$results, patientID = patientID, outDir=outDir)
outFile <- paste0(outDir, "/", patientID, ".params.txt")
outputParametersToFile(hmmResults.cor, file = outFile)

## plot solutions for all samples 
plotSolutions(hmmResults.cor, tumour_copy, chrs, outDir, numSamples=numSamples,
              plotFileType=plotFileType, plotYLim=plotYLim, 
              estimateScPrevalence=estimateScPrevalence, maxCN=maxCN)

