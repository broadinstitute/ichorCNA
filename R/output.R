# file:   output.R
# author: Gavin Ha, Ph.D.
#		Fred Hutchinson Cancer Research Center
# contact: <gha@fredhutch.org>
# author: Justin Rhoades
#       Dana-Farber Cancer Institute
#       Broad Institute
# https://github.com/broadinstitute/ichorCNA
# date:  July 12, 2019
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.

##################################################
###### FUNCTION TO GET OUTPUT HMM RESULTS ########
##################################################
outputHMM <- function(cna, segs, results, patientID = NULL, outDir = "."){
  names <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))

  S <- results$param$numberSamples
  
  segout <- NULL
  shuffle <- NULL
  if (is.null(patientID)){
    patientID <- names(cna)[1]
  }
  ### PRINT OUT SEGMENTS FOR EACH SAMPLE TO SEPARATE FILES ###
  for (s in 1:S){
    id <- names(cna)[s]
    bin_size <- as.numeric(cna[[s]][1,"end"]) - as.numeric(cna[[s]][1,"start"]) + 1
    markers <- (segs[[s]]$end - segs[[s]]$start + 1) / bin_size  
    segTmp <- cbind(sample = as.character(id), segs[[s]][, 1:3],
                    event = names[segs[[s]]$copy.number + 1],
                    copy.number = segs[[s]]$copy.number, 
                    bins = markers, median = segs[[s]]$median,
                    subclone.status=segs[[s]]$subclone.status,
                    segs[[s]][, 9:ncol(segs[[s]])])
    segout <- rbind(segout, segTmp[, 1:8])
    ### Re-ordering the columns for output ###
    shuffleTmp <- segTmp[, c(1, 2, 3, 4, 7, 8, 6, 5, 9:ncol(segTmp))]
    colnames(shuffleTmp)[1:9] <- c("ID", "chrom", "start", "end", 
                              "num.mark", "seg.median.logR", "copy.number", "call", "subclone.status")
    shuffle <- rbind(shuffle, shuffleTmp) 
  }
  
  shu_out = paste(outDir,"/", patientID,".seg.txt",sep="")
  seg_out = paste(outDir,"/", patientID,".seg",sep="")
  message("Writing segments to ", seg_out)
  write.table(shuffle, file = shu_out, quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(segout, file = seg_out, quote = FALSE, sep = "\t", row.names = FALSE)
  ### PRINT OUT BIN-LEVEL DATA FOR ALL SAMPLES IN SAME FILE ##
  cnaout <- cbind(cna[[1]][, -c(1)])
  colnames(cnaout)[4:ncol(cnaout)] <- paste0(names(cna)[1], ".", colnames(cnaout)[4:ncol(cnaout)])
  if (S >= 2) {
    for (s in 2:S){
      id <- names(cna)[s]
      cnaTmp <- cna[[s]][, -c(1)]
      colnames(cnaTmp)[4:ncol(cnaTmp)] <- paste0(id, ".", colnames(cnaTmp)[4:ncol(cnaTmp)])
      cnaout <- merge(cnaout, cnaTmp, by = c("chr", "start", "end"))
    }
  }
  cnaout$chr <- factor(cnaout$chr, levels = unique(cna[[1]]$chr))
  cnaout <- cnaout[order(cnaout[, 1], cnaout[, 2], cnaout[, 3]), ]
  cna_out = paste(outDir,"/", patientID,".cna.seg",sep="")
  message("Outputting to bin-level results to ", cna_out)
  write.table(cnaout, file = cna_out, quote = FALSE, sep = "\t", row.names = FALSE)
}

outputParametersToFile <- function(hmmResults, file){
  S <- hmmResults$results$param$numberSamples
  x <- hmmResults$results
  i <- x$iter
  fc <- file(file, "w+")
  outMat <- as.data.frame(cbind(`Tumor Fraction` = signif(1 - x$n[, i], digits = 4), Ploidy = signif(x$phi[, i], digits = 4)))
  outMat <- cbind(Sample = names(hmmResults$cna), outMat)
  rownames(outMat) <- names(hmmResults$cna)
  write.table(outMat, file = fc, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table("\n", file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  for (s in 1:S){
    id <- names(hmmResults$cna)[s]
    write.table(id, file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(paste0("Gender:\t", x$gender), file = fc, col.names = FALSE, 
                row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(paste0("Tumor Fraction:\t", signif(1 - x$n[s, i], digits = 4)), file = fc, col.names = FALSE, 
                row.names = FALSE, quote = FALSE, sep = "\t")
    ploidy <- (1 - x$n[i]) * x$phi[i] + x$n[i] * 2
    write.table(paste0("Ploidy:\t", signif(x$phi[s, i], digits = 4)), file = fc, col.names = FALSE, 
                row.names = FALSE, quote = FALSE, sep = "\t")
   	subcloneGenomeFrac <- sum(hmmResults.cor$cna[[s]]$subclone.status) / nrow(hmmResults.cor$cna[[s]])
   	subcloneCNAFrac <- sum(hmmResults.cor$cna[[s]]$subclone.status) / sum(hmmResults.cor$cna[[s]]$copy.num != 2)
   	scFrac <- 1 - x$sp[s, i]
   	if (subcloneGenomeFrac == 0){
   		scFrac <- NA
   	}else{
   		scFrac <- signif(scFrac, digits = 4)
   	}
    #if (sum(x$param$ct.sc.status) != 0){    	
      write.table(paste0("Subclone Fraction:\t", scFrac), file = fc, col.names = FALSE, 
                  row.names = FALSE, quote = FALSE, sep = "\t")
      write.table(paste0("Fraction Genome Subclonal:\t", signif(subcloneGenomeFrac, digits = 4)), file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
      write.table(paste0("Fraction CNA Subclonal:\t", signif(subcloneCNAFrac, digits = 4)), file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    #}
    if (is.null(x$coverage)){
    	coverage <- NA
    }else{
    	coverage <- signif(coverage, digits = 4)
    }
		write.table(paste0("Coverage:\t", coverage), file = fc, col.names = FALSE, 
								row.names = FALSE, quote = FALSE, sep = "\t")
	
    if (!is.null(x$chrYCov)){
      write.table(paste0("ChrY coverage fraction:\t", signif(x$chrYCov[s], digits = 4)), file = fc, col.names = FALSE, 
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
    if (!is.null(x$chrXMedian)){
      write.table(paste0("ChrX median log ratio:\t", signif(x$chrXMedian[s], digits = 4)), file = fc, col.names = FALSE, 
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
    write.table(paste0("Student's t mean: ", paste0(signif(x$mus[,s,i], digits = 2), collapse = ", ")), 
                file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(paste0("Student's t precision: ", paste0(signif(x$lambdas[,s,i], digits = 2), collapse = ", ")), file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(paste0("Gamma Rate Init:\t", signif(hmmResults.cor$results$param$betaLambda[1], digits=2)), file=fc, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
    write.table(paste0("GC-Map correction MAD:\t", format(mad(diff(2^as.numeric(hmmResults$cna[[s]][,"logR"])), na.rm=T), digits=4)), file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    write.table("\n", file = fc, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  write.table(x$loglik, file = fc, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  close(fc)
  invisible()
}
