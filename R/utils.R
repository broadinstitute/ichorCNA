# file:   utils.R
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

####################################
##### FUNCTION TO FILTER CHRS ######
####################################
keepChr <- function(tumour_reads, chr = c(1:22,"X","Y")){	
	tumour_reads <- tumour_reads[space(tumour_reads) %in% chr, ]
	tumour_reads <- as.data.frame(tumour_reads)
	tumour_reads$space <- droplevels(tumour_reads$space)
	tumour_reads$space <- factor(tumour_reads$space,levels=chr)
	tumour_reads <- as(tumour_reads,"RangedData")
	return(tumour_reads)
}

filterEmptyChr <- function(gr){
	require(plyr)
	ind <- daply(as.data.frame(gr), .variables = "seqnames", .fun = function(x){
	  rowInd <- apply(x[, 6:ncol(x), drop = FALSE], 1, function(y){
	    sum(is.na(y)) == length(y)
	  })
	  sum(rowInd) == nrow(x)
	})	
	return(keepSeqlevels(gr, value = names(which(!ind))))
}


##################################################
##### FUNCTION TO FILTER CENTROMERE REGIONS ######
##################################################
excludeCentromere <- function(x, centromere, flankLength = 0, genomeStyle = "NCBI"){
	require(GenomeInfoDb)
	colnames(centromere)[1:3] <- c("space","start","end")
	centromere$space <- setGenomeStyle(centromere$space, genomeStyle)
	centromere$start <- centromere$start - flankLength
	centromere$end <- centromere$end + flankLength
	centromere$space <- factor(centromere$space, levels = centromere$space)
	centromere <- as(centromere, "RangedData")
	
	hits <- findOverlaps(query = x, subject = centromere)
	ind <- queryHits(hits)
	message("Removed ", length(ind), " bins near centromeres.")
	return(x[-ind, ])
}

##################################################
##### FUNCTION TO USE NCBI CHROMOSOME NAMES ######
##################################################
setGenomeStyle <- function(x, genomeStyle = "NCBI", species = "Homo_sapiens"){
        require(GenomeInfoDb)
        #chrs <- genomeStyles(species)[c("NCBI","UCSC")]
        if (!genomeStyle %in% seqlevelsStyle(as.character(x))){
        x <- suppressWarnings(mapSeqlevels(as.character(x),
                                        genomeStyle, drop = FALSE)[1,])
    }

    autoSexMChr <- extractSeqlevelsByGroup(species = species,
                                style = genomeStyle, group = "all")
    x <- x[x %in% autoSexMChr]
    return(x)
}

loadReadCountsFromWig <- function(counts, chrs = c(1:22, "X", "Y"), gc = NULL, map = NULL, centromere = NULL, flankLength = 100000, targetedSequences = NULL, genomeStyle = "NCBI", applyCorrection = TRUE, mapScoreThres = 0.9, chrNormalize = c(1:22, "X", "Y"), fracReadsInChrYForMale = 0.002, useChrY = TRUE){
	require(HMMcopy)
	require(GenomeInfoDb)
	counts.raw <- counts
	names(counts) <- setGenomeStyle(names(counts), genomeStyle)
	counts <- keepChr(counts, chrs)
	if (!is.null(gc)){ 
		names(gc) <- setGenomeStyle(names(gc), genomeStyle)
		counts$gc <- keepChr(gc, chrs)$value
	}
	if (!is.null(map)){ 
		names(map) <- setGenomeStyle(names(map), genomeStyle)
		counts$map <- keepChr(map, chrs)$value
	}
	colnames(counts)[1] <- c("reads")
	
	# remove centromeres
	if (!is.null(centromere)){ 
		counts <- excludeCentromere(counts, centromere, flankLength = flankLength)
	}
	# keep targeted sequences
	if (!is.null(targetedSequences)){
		countsExons <- filterByTargetedSequences(counts, targetedSequences)
		counts <- counts[countsExons$ix,]
	}
	gender <- NULL
	if (applyCorrection){
	## correct read counts ##
    counts <- correctReadCounts(counts, chrNormalize = chrNormalize)
    if (!is.null(map)) {
      ## filter bins by mappability
      counts <- filterByMappabilityScore(counts, map=map, mapScoreThres = mapScoreThres)
    }
    ## get gender ##
    gender <- getGender(counts.raw, counts, gc, map, fracReadsInChrYForMale = fracReadsInChrYForMale, useChrY = useChrY,
                        centromere=centromere, flankLength=flankLength, targetedSequences = targetedSequences)
    }
  return(list(counts = counts, gender = gender))
}

filterByMappabilityScore <- function(counts, map, mapScoreThres = 0.9){
	message("Filtering low uniqueness regions with mappability score < ", mapScoreThres)
	counts <- counts[counts$map >= mapScoreThres, ]
	return(counts)
}

filterByTargetedSequences <- function(counts, targetedSequences){
 ### for targeted sequencing (e.g.  exome capture),
    ### ignore bins with 0 for both tumour and normal
    ### targetedSequence = RangedData (IRanges object)
    ### containing list of targeted regions to consider;
    ### 3 columns: chr, start, end
	message("Analyzing targeted regions...")
	targetIR <- RangedData(ranges = IRanges(start = targetedSequences[, 2], 
				end = targetedSequences[, 3]), space = targetedSequences[, 1])
				
	hits <- findOverlaps(query = counts, subject = targetIR)
	keepInd <- unique(queryHits(hits))    

	return(list(counts=counts, ix=keepInd))
}

selectFemaleChrXSolution <- function(){
	
}

##################################################
### FUNCTION TO DETERMINE GENDER #################
##################################################
getGender <- function(rawReads, normReads, gc, map, fracReadsInChrYForMale = 0.002, useChrY = TRUE,
											centromere=NULL, flankLength=1e5, targetedSequences=NULL){
	chrXInd <- normReads$space == "X"
	if (sum(chrXInd) > 1){ ## if no X 
		chrXMedian <- median(normReads[chrXInd, ]$copy, na.rm = TRUE)
		# proportion of reads in chrY #
		tumY <- loadReadCountsFromWig(rawReads, chrs="Y", gc=gc, map=map, applyCorrection = FALSE,
				centromere=centromere, flankLength=flankLength, targetedSequences=targetedSequences)$counts
		chrYCov <- sum(tumY$reads) / sum(rawReads$value)
		if (chrXMedian < -0.5){
			if (useChrY && (chrYCov < fracReadsInChrYForMale)){ #trumps chrX if using chrY
					gender <- "female"  
			}else{
				gender <- "male" # satisfies decreased chrX log ratio and/or increased chrY coverage
			}
		}else{
			gender <- "female" # chrX is provided but does not satisfies male critera
		}
	}else{
		gender <- "unknown" # chrX is not provided
		chrYCov <- NA
		chrXMedian <- NULL
	}
	return(list(gender=gender, chrYCovRatio=chrYCov, chrXMedian=chrXMedian))
}
	
	
normalizeByPanelOrMatchedNormal <- function(tumour_copy, chrs = c(1:22, "X", "Y"), 
      normal_panel = NULL, normal_copy = NULL, gender = "female", normalizeMaleX = FALSE){
 	### COMPUTE LOG RATIO FROM MATCHED NORMAL OR PANEL AND HANDLE CHRX ###
	## NO PANEL
	# matched normal but NO panel, then just normalize by matched normal (WES)
	## WHY DO WE NOT NORMALIZE BY NORMAL WITH PANEL? ##
	chrXInd <- tumour_copy$space == "X"
	chrXMedian <- median(tumour_copy[chrXInd, ]$copy, na.rm = TRUE)
	if (!is.null(normal_copy) && is.null(normal_panel)){
			message("Normalizing Tumour by Normal")
			tumour_copy$copy <- tumour_copy$copy - normal_copy$copy
			rm(normal_copy)
	}
	# matched normal and panel and male, then compute normalized chrX median (WES)
	if (!is.null(normal_copy) && !is.null(normal_panel) && gender=="male"){
			message("Normalizing by matched normal for ChrX")
			chrX.MNnorm <- tumour_copy$copy[chrXInd] - normal_copy$copy[chrXInd]
			chrXMedian.MNnorm <- median(chrX.MNnorm, na.rm = TRUE)
	}
	# if male, then just normalize chrX to median (ULP and WES)
	if (is.null(normal_copy) && gender=="male" && !gender.mismatch && normalizeMaleX){
			tumour_copy$copy[chrXInd] <- tumour_copy$copy[chrXInd] - chrXMedian
	}
	# PANEL, then normalize by panel instead of matched normal (ULP and WES)
	if (!is.null(normal_panel)){
		panel <- readRDS(normal_panel) ## load in IRanges object
		panel <- keepChr(panel, chr = chrs)
        # intersect bins in sample and panel
        hits <- findOverlaps(tumour_copy, panel, type="equal")
        tumour_copy <- tumour_copy[queryHits(hits),]
        panel <- panel[subjectHits(hits),]
        # subtract out panel median
		tumour_copy$copy <- tumour_copy$copy - panel$Median
		# if male, then shift chrX by +chrXMedian.MNnorm
		if (gender == "male" && exists("chrXMedian.MNnorm")){
			tumour_copy$copy[chrXInd] <- tumour_copy$copy[chrXInd] + chrXMedian.MNnorm
		}
	}
	return(tumour_copy)
}

##################################################
###### FUNCTION TO CORRECT GC/MAP BIASES ########
##################################################
correctReadCounts <- function(x, chrNormalize = c(1:22), mappability = 0.9, samplesize = 50000,
    verbose = TRUE) {
  if (length(x$reads) == 0 | length(x$gc) == 0) {
    stop("Missing one of required columns: reads, gc")
  }
	chrInd <- space(x) %in% chrNormalize
  if(verbose) { message("Applying filter on data...") }
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- quantile(x$reads[x$valid & chrInd], prob = c(0, 1 - routlier), na.rm = TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid & chrInd], prob = c(doutlier, 1 - doutlier), na.rm = TRUE)
  if (length(x$map) != 0) {
    x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
      x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  } else {
    x$ideal[!x$valid | x$reads <= range[1] |
      x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE
  }

  if (verbose) { message("Correcting for GC bias...") }
  set <- which(x$ideal & chrInd)
  select <- sample(set, min(length(set), samplesize))
  rough = loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final = loess(predict(rough, i) ~ i, span = 0.3)
  x$cor.gc <- x$reads / predict(final, x$gc)

  if (length(x$map) != 0) {
    if (verbose) { message("Correcting for mappability bias...") }
    coutlier <- 0.01
    range <- quantile(x$cor.gc[which(x$valid & chrInd)], prob = c(0, 1 - coutlier), na.rm = TRUE)
    set <- which(x$cor.gc < range[2] & chrInd)
    select <- sample(set, min(length(set), samplesize))
    final = approxfun(lowess(x$map[select], x$cor.gc[select]))
    x$cor.map <- x$cor.gc / final(x$map)
  } else {
    x$cor.map <- x$cor.gc
  }
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  return(x)
}



computeBIC <- function(params){
  iter <- params$iter
  KS <- nrow(params$rho) # num states
  N <- ncol(params$rho) # num data points
  NP <- nrow(params$n) + nrow(params$phi) # normal + ploidy
  L <- 1 # precision (lambda)
  numParams <- KS * (KS - 1) + KS * (L + NP) - 1
  b <- -2 * params$loglik[iter] + numParams * log(N)
  return(b)
} 


#############################################################
## function to compute power from ULP-WGS purity/ploidy #####
#############################################################
# current assumptions: power for clonal heterozygous mutations

