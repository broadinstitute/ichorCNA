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
# updated for GRanges #
keepChr <- function(tumour_reads, chrs = c(1:22,"X","Y")){	
	tumour_reads <- keepSeqlevels(tumour_reads, chrs, pruning.mode="coarse")
	sortSeqlevels(tumour_reads)
	return(sort(tumour_reads))
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

####################################
##### FUNCTION GET SEQINFO ######
####################################
getSeqInfo <- function(genomeBuild = "hg19", genomeStyle = "NCBI"){
	bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
	if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
		seqinfo <- Seqinfo(genome=genomeBuild)
	} else {
		seqinfo <- seqinfo(get(bsg))
	}
	seqlevelsStyle(seqinfo) <- genomeStyle
	seqinfo <- keepSeqlevels(seqinfo, value = chrs)
	#seqinfo <- cbind(seqnames = seqnames(seqinfo), as.data.frame(seqinfo))
	return(seqinfo)	
}

##################################################
##### FUNCTION TO FILTER CENTROMERE REGIONS ######
##################################################
excludeCentromere <- function(x, centromere, flankLength = 0, genomeStyle = "NCBI"){
	require(GenomeInfoDb)
	colnames(centromere)[1:3] <- c("seqnames","start","end")
	centromere$start <- centromere$start - flankLength
	centromere$end <- centromere$end + flankLength
	centromere <- as(centromere, "GRanges")
	seqlevelsStyle(centromere) <- genomeStyle
	centromere <- sort(centromere)	
	hits <- findOverlaps(query = x, subject = centromere)
	ind <- queryHits(hits)
	message("Removed ", length(ind), " bins near centromeres.")
	if (length(ind) > 0){
		x <- x[-ind, ]
	}
	return(x)
}

##################################################
##### FUNCTION TO USE NCBI CHROMOSOME NAMES ######
##################################################
## deprecated ##
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

wigToGRanges <- function(wigfile, verbose = TRUE){
  if (verbose) { message(paste("Slurping:", wigfile)) }
  input <- readLines(wigfile, warn = FALSE)
  breaks <- c(grep("fixedStep", input), length(input) + 1)
  temp <- NULL
  span <- NULL
  for (i in 1:(length(breaks) - 1)) {
    data_range <- (breaks[i] + 1):(breaks[i + 1] - 1)
    track_info <- input[breaks[i]]
    if (verbose) { message(paste("Parsing:", track_info)) }
    tokens <- strsplit(
      sub("fixedStep chrom=(\\S+) start=(\\d+) step=(\\d+) span=(\\d+)",
          "\\1 \\2 \\3 \\4", track_info, perl = TRUE), " ")[[1]]
    span <- as.integer(tokens[4])
    chr <- rep.int(tokens[1], length(data_range))
    pos <- seq(from = as.integer(tokens[2]), by = as.integer(tokens[3]),
               length.out = length(data_range))
    val <- as.numeric(input[data_range])
    temp <- c(temp, list(data.frame(chr, pos, val)))
  }
  if (verbose) { message("Sorting by decreasing chromosome size") }
  lengths <- as.integer(lapply(temp, nrow))
  temp <- temp[order(lengths, decreasing = TRUE)]
  temp = do.call("rbind", temp)
  output <- GenomicRanges::GRanges(ranges = IRanges(start = temp$pos, width = span),
                       seqnames = temp$chr, value = temp$val)
  return(output)
}


loadReadCountsFromWig <- function(counts, chrs = c(1:22, "X", "Y"), gc = NULL, map = NULL, centromere = NULL, flankLength = 100000, targetedSequences = NULL, genomeStyle = "NCBI", applyCorrection = TRUE, mapScoreThres = 0.9, chrNormalize = c(1:22, "X", "Y"), fracReadsInChrYForMale = 0.002, chrXMedianForMale = -0.5, useChrY = TRUE){
	require(HMMcopy)
	require(GenomeInfoDb)
	seqlevelsStyle(counts) <- genomeStyle
	counts.raw <- counts	
	counts <- keepChr(counts, chrs)
	
	if (!is.null(gc)){ 
		seqlevelsStyle(gc) <- genomeStyle
		counts$gc <- keepChr(gc, chrs)$value
	}
	if (!is.null(map)){ 
		seqlevelsStyle(map) <- genomeStyle
		counts$map <- keepChr(map, chrs)$value
	}
	colnames(values(counts))[1] <- c("reads")
	
	# remove centromeres
	if (!is.null(centromere)){ 
		counts <- excludeCentromere(counts, centromere, flankLength = flankLength, genomeStyle=genomeStyle)
	}
	# keep targeted sequences
	if (!is.null(targetedSequences)){
		colnames(targetedSequences)[1:3] <- c("chr", "start", "end")
		targetedSequences.GR <- as(targetedSequences, "GRanges")
		seqlevelsStyle(targetedSequences.GR) <- genomeStyle
		countsExons <- filterByTargetedSequences(counts, targetedSequences.GR)
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
		gender <- getGender(counts.raw, counts, gc, map, fracReadsInChrYForMale = fracReadsInChrYForMale, 
							chrXMedianForMale = chrXMedianForMale, useChrY = useChrY,
							centromere=centromere, flankLength=flankLength, targetedSequences = targetedSequences,
							genomeStyle = genomeStyle)
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
    ### targetedSequence = GRanges object
    ### containing list of targeted regions to consider;
    ### 3 columns: chr, start, end
					
	hits <- findOverlaps(query = counts, subject = targetedSequences)
	keepInd <- unique(queryHits(hits))    

	return(list(counts=counts, ix=keepInd))
}

selectFemaleChrXSolution <- function(){
	
}

##################################################
### FUNCTION TO DETERMINE GENDER #################
##################################################
getGender <- function(rawReads, normReads, gc, map, fracReadsInChrYForMale = 0.002, chrXMedianForMale = -0.5, useChrY = TRUE,
					  centromere=NULL, flankLength=1e5, targetedSequences=NULL, genomeStyle="NCBI"){
	chrXStr <- grep("X", runValue(seqnames(normReads)), value = TRUE)
	chrYStr <- grep("Y", runValue(seqnames(rawReads)), value = TRUE)
	chrXInd <- as.character(seqnames(normReads)) == chrXStr
	if (sum(chrXInd) > 1){ ## if no X 
		chrXMedian <- median(normReads[chrXInd, ]$copy, na.rm = TRUE)
		# proportion of reads in chrY #
		tumY <- loadReadCountsFromWig(rawReads, chrs=chrYStr, genomeStyle=genomeStyle,
				gc=gc, map=map, applyCorrection = FALSE, centromere=centromere, flankLength=flankLength, 
				targetedSequences=targetedSequences)$counts
		chrYCov <- sum(tumY$reads) / sum(rawReads$value)
		if (chrXMedian < chrXMedianForMale){
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
    genomeStyle <- seqlevelsStyle(tumour_copy)[1]
    seqlevelsStyle(chrs) <- genomeStyle
 	### COMPUTE LOG RATIO FROM MATCHED NORMAL OR PANEL AND HANDLE CHRX ###
	## NO PANEL
	# matched normal but NO panel, then just normalize by matched normal (WES)
	## WHY DO WE NOT NORMALIZE BY NORMAL WITH PANEL? ##
	chrXInd <- grep("X", as.character(seqnames(tumour_copy)))
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
		## load in IRanges object, then convert to GRanges
		panel <- readRDS(normal_panel)
		seqlevelsStyle(panel) <- genomeStyle
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
  chrInd <- as.character(seqnames(x)) %in% chrNormalize
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

## Recompute integer CN for high-level amplifications ##
## compute logR-corrected copy number ##
correctIntegerCN <- function(cn, segs, callColName = "event", 
		purity, ploidy, cellPrev, maxCNtoCorrect.autosomes = NULL, 
		maxCNtoCorrect.X = NULL, correctHOMD = TRUE, minPurityToCorrect = 0.2, gender = "male", chrs = c(1:22, "X")){
	names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP", rep("HLAMP", 1000))

	## set up chromosome style
	autosomeStr <- grep("X|Y", chrs, value=TRUE, invert=TRUE)
	chrXStr <- grep("X", chrs, value=TRUE)
	
	if (is.null(maxCNtoCorrect.autosomes)){
		maxCNtoCorrect.autosomes <- max(segs[segs$chr %in% autosomeStr, "copy.number"], na.rm = TRUE)
	}
	if (is.null(maxCNtoCorrect.X) & gender == "female" & length(chrXStr) > 0){
		maxCNtoCorrect.X <- max(segs[segs$chr == chrXStr, "copy.number"], na.rm=TRUE)
	}
	## correct log ratio and compute corrected CN
	cellPrev.seg <- rep(1, nrow(segs))
	cellPrev.seg[as.logical(segs$subclone.status)] <- cellPrev
	segs$logR_Copy_Number <- logRbasedCN(segs[["median"]], purity, ploidy, cellPrev.seg, cn=2)
	if (gender == "male" & length(chrXStr) > 0){ ## analyze chrX separately
		ind.cnChrX <- which(segs$chr == chrXStr)
		segs$logR_Copy_Number[ind.cnChrX] <- logRbasedCN(segs[["median"]][ind.cnChrX], purity, ploidy, cellPrev.seg[ind.cnChrX], cn=1)
	}

	## assign copy number to use - Corrected_Copy_Number
	# 1) ame ichorCNA calls for autosomes - no change in copy number
	segs$Corrected_Copy_Number <- as.integer(segs$copy.number)
	segs$Corrected_Call <- segs[[callColName]]

	ind.change <- c()
	if (purity >= minPurityToCorrect){
		# 2) ichorCNA calls adjusted for >= copies - HLAMP
		# perform on all chromosomes
		ind.cn <- which(segs$copy.number >= maxCNtoCorrect.autosomes | 
						(segs$logR_Copy_Number >= maxCNtoCorrect.autosomes * 1.2 & !is.infinite(segs$logR_Copy_Number)))
		segs$Corrected_Copy_Number[ind.cn] <- as.integer(round(segs$logR_Copy_Number[ind.cn]))
		segs$Corrected_Call[ind.cn] <- names[segs$Corrected_Copy_Number[ind.cn] + 1]
		ind.change <- c(ind.change, ind.cn)
		
		# 3) ichorCNA calls adjust for HOMD
		if (correctHOMD){
			ind.cn <- which(segs$chr %in% chrs & 
				(segs$copy.number == 0 | segs$logR_Copy_Number == 1/2^6))
			segs$Corrected_Copy_Number[ind.cn] <- as.integer(round(segs$logR_Copy_Number[ind.cn]))
			segs$Corrected_Call[ind.cn] <- names[segs$Corrected_Copy_Number[ind.cn] + 1]
			ind.change <- c(ind.change, ind.cn)
		}
		# 4) Re-adjust chrX copy number for males (females already handled above)
		if (gender == "male" & length(chrXStr) > 0){
			ind.cn <- which(segs$chr == chrXStr & 
				(segs$copy.number >= maxCNtoCorrect.X | segs$logR_Copy_Number >= maxCNtoCorrect.X * 1.2))
			segs$Corrected_Copy_Number[ind.cn] <- as.integer(round(segs$logR_Copy_Number[ind.cn]))
			segs$Corrected_Call[ind.cn] <- names[segs$Corrected_Copy_Number[ind.cn] + 2]
			ind.change <- c(ind.change, ind.cn)
		}
	}
	## adjust the bin level data ##
	# 1) assign the original calls
	cn$Corrected_Copy_Number <- as.integer(cn$copy.number)
	cn$Corrected_Call <- cn[[callColName]]
	cellPrev.cn <- rep(1, nrow(cn))
	cellPrev.cn[as.logical(cn$subclone.status)] <- cellPrev
	cn$logR_Copy_Number <- logRbasedCN(cn[["logR"]], purity, ploidy, cellPrev.cn, cn=2)
	if (gender == "male" & length(chrXStr) > 0){ ## analyze chrX separately
		ind.cnChrX <- which(cn$chr == chrXStr)
		cn$logR_Copy_Number[ind.cnChrX] <- logRbasedCN(cn[["logR"]][ind.cnChrX], purity, ploidy, cellPrev.cn[ind.cnChrX], cn=1)
	}
	if (purity >= minPurityToCorrect){
		# 2) correct bins overlapping adjusted segs
		ind.change <- unique(ind.change)
		ind.overlapSegs <- c()
		if (length(ind.change) > 0){		
			cn.gr <- as(cn, "GRanges")
			segs.gr <- as(segs, "GRanges")
			hits <- findOverlaps(query = cn.gr, subject = segs.gr[ind.change])
			cn$Corrected_Copy_Number[queryHits(hits)] <- segs$Corrected_Copy_Number[ind.change][subjectHits(hits)]
			cn$Corrected_Call[queryHits(hits)] <- segs$Corrected_Call[ind.change][subjectHits(hits)]
			ind.overlapSegs <- queryHits(hits)
		}
		# 3) correct bins that are missed as high level amplifications
		ind.hlamp <- which(cn$copy.number >= maxCNtoCorrect.autosomes | 
	 					(cn$logR_Copy_Number >= maxCNtoCorrect.autosomes * 1.2 & !is.infinite(cn$logR_Copy_Number)))
		ind.cn <- unique(ind.hlamp, ind.overlapSegs)
	 	cn$Corrected_Copy_Number[ind.cn] <- as.integer(round(cn$logR_Copy_Number[ind.cn]))
	 	cn$Corrected_Call[ind.cn] <- names[cn$Corrected_Copy_Number[ind.cn] + 1]
	 }
	 
	return(list(cn = cn, segs = segs))
}


## compute copy number using corrected log ratio ##
logRbasedCN <- function(x, purity, ploidyT, cellPrev=NA, cn = 2){
	if (length(cellPrev) == 1 && is.na(cellPrev)){
		cellPrev <- 1
	}else{ #if cellPrev is a vector
		cellPrev[is.na(cellPrev)] <- 1
	}
	ct <- (2^x 
		* (cn * (1 - purity) + purity * ploidyT * (cn / 2)) 
		- (cn * (1 - purity)) 
		- (cn * purity * (1 - cellPrev))) 
	ct <- ct / (purity * cellPrev)
	ct <- sapply(ct, max, 1/2^6)
	return(ct)
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

