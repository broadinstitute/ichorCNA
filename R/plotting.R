# file:   plotting.R
# author: Gavin Ha, Ph.D.
#         Fred Hutchinson Cancer Research Center
# contact: <gha@fredhutch.org>
# 
# author: Justin Rhoades, Broad Institute
#
# ichorCNA website: https://github.com/broadinstitute/ichorCNA
# date:   August 12, 2019
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This file contains R functions for plotting.

## plot solutions for all samples
plotSolutions <- function(hmmResults.cor, tumour_copy, logR.column = "logR", call.column = "event", 
					      chrs, outDir, plotSegs = TRUE,
                          numSamples=1, plotFileType="pdf", plotYLim=c(-2,2), seqinfo = NULL,
                          estimateScPrevalence=FALSE, maxCN){
  ## for each sample ##
  for (s in 1:numSamples){
    iter <- hmmResults.cor$results$iter
    id <- names(hmmResults.cor$cna)[s]
    ploidyEst <- hmmResults.cor$results$phi[s, iter]
    normEst <- hmmResults.cor$results$n[s, iter]
    purityEst <- 1 - normEst
    ploidyAll <- (1 - normEst) * ploidyEst + normEst * 2

    outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide")
    plotGWSolution(hmmResults.cor, s, outPlotFile, seqinfo, 
    			   logR.column = logR.column, call.column = call.column, 
    			   plotFileType=plotFileType, plotYLim=plotYLim,
                   plotSegs=plotSegs, estimateScPrevalence=estimateScPrevalence, main=id)

    ### PLOT THE LOG RATIO DATA ALONG WITH COLOUR-CODING FOR PREDICTED CALLS ###
    for (i in chrs){
      ## PLOT CNA BY CHROMOSOME ##
      outPlot <- paste0(outDir,"/",id,"/",id,"_CNA_chr",i)
      if (plotFileType == "png"){ 
        outPlot <- paste0(outPlot, ".png")
        png(outPlot,width=15,height=5,units="in",res=300)
      }else{
        outPlot <- paste0(outPlot, ".pdf")
        pdf(outPlot,width=15,height=5)
      }			
      par(mfrow=c(1,1))
      plotCNlogRByChr(hmmResults.cor$cna[[s]], segs=hmmResults.cor$results$segs[[s]], chr=i,
                      ploidy = ploidyAll, plotSegs=plotSegs, seqinfo = seqinfo,
                      logR.column = logR.column, call.column = call.column, 
                      cytoBand=T, yrange=plotYLim, cex=0.75, spacing=8)	
      dev.off()
    }

    ### PLOT THE CORRECTION COMPARISONS ###
    outPlotFile <- paste0(outDir,"/",id,"/",id,"_correct")
    if (plotFileType == "png"){ 
      outPlotFile <- paste0(outPlotFile, ".png")
      png(outPlotFile,width=10,height=12,units="in",res=300)
    }else{
      outPlotFile <- paste0(outPlotFile, ".pdf")
      pdf(outPlotFile,width=10,height=12)
    }
    plotCorrectionGenomeWide(tumour_copy[[s]], seqinfo = seqinfo, pch = ".")
    dev.off()

    ### PLOT THE BIAS ###
    outPlotFile <- paste0(outDir,"/",id,"/",id,"_bias")
    if (plotFileType == "png"){ 
      outPlotFile <- paste0(outPlotFile, ".png")
      png(outPlotFile,width=7,height=7,units="in",res=300)
    }else{
      outPlotFile <- paste0(outPlotFile, ".pdf")
      pdf(outPlotFile,width=7,height=7)
    }
    try(plotBias(tumour_copy[[s]], pch = 20, cex = 0.5), silent=TRUE)
    dev.off()

    ### PLOT TPDF ##
    outPlotFile <- paste0(outDir,"/",id,"/",id,"_tpdf.pdf")
    pdf(outPlotFile)
    plotParam(mus = unique(hmmResults.cor$results$mus[, s, iter]), 
              lambdas = hmmResults.cor$results$lambdas[, s, iter], 
              subclone = hmmResults.cor$results$param$ct.sc.status,
              nu = hmmResults.cor$results$param$nu, copy.states = 1:maxCN)
    dev.off()

  }
}


plotGWSolution <- function(hmmResults.cor, s, outPlotFile, plotFileType="pdf", 
						   logR.column = "logR", call.column = "event",
						   seqinfo = NULL, plotSegs = TRUE, 
                           plotYLim=c(-2,2), estimateScPrevalence, main,
                           turnDevOn=TRUE, turnDevOff=TRUE){
    ## plot genome wide figures for each solution ##
    iter <- hmmResults.cor$results$iter
    ploidyEst <- hmmResults.cor$results$phi[s, iter]
    normEst <- hmmResults.cor$results$n[s, iter]
    purityEst <- 1 - normEst
    ploidyAll <- (1 - normEst) * ploidyEst + normEst * 2
    subclone <- 1 - hmmResults.cor$results$sp[s, iter]
    #outPlotFile <- paste0(outDir, "/", id, "/", id, "_genomeWide")
    if (turnDevOn){
      if (plotFileType == "png"){ 
          outPlotFile <- paste0(outPlotFile, ".png")
          png(outPlotFile,width=20,height=6,units="in",res=300)
      }else{
          outPlotFile <- paste0(outPlotFile, ".pdf")
          pdf(outPlotFile,width=20,height=6)
      }	
    }
    if (plotSegs){
    	segsToUse <- hmmResults.cor$results$segs[[s]]
    }else{
    	segsToUse <- NULL
    }
    plotCNlogRByChr(hmmResults.cor$cna[[s]], segs = segsToUse, plotSegs=plotSegs, seqinfo=seqinfo,
                    param = hmmResults.cor$results$param, chr=NULL,
                    logR.column = logR.column, call.column = call.column,  
                    ploidy = ploidyAll, cytoBand=T, yrange=plotYLim, main=main)  #ylim for plot
    annotStr <- paste0("Tumor Fraction: ", signif(purityEst, digits=4), ", Ploidy: ", signif(ploidyEst, digits=3))
    if (!is.null(coverage)){
      annotStr <- paste0(annotStr, ", Coverage: ", signif(coverage, digits=2))
    }
    mtext(line=-1, annotStr, cex=1.5)
    if (estimateScPrevalence){
      sampleBins <- hmmResults.cor$cna[[s]]
      #sampleBins <- sampleBins[sampleBins$chr %in% chrTrain, ]
      subBinCount <- sum(sampleBins$subclone.status) 
      fracGenomeSub <- subBinCount / nrow(sampleBins)
      fracCNAsub <- subBinCount / sum(sampleBins$copy.number != 2)
      if (fracGenomeSub > 0){
        annotSubStr <- paste0("Subclone Fraction: ", signif(subclone, digits=3), 
                              ", Frac. Genome Subclonal: ", format(round(fracGenomeSub, 2), nsmall = 2),
                              ", Frac. CNA Subclonal: ", format(round(fracCNAsub, 2), nsmall = 2))
        mtext(line=-2, annotSubStr, cex=1.5)
      }
    }
    if (turnDevOff){
        dev.off()
    }
}


#data is the output format of HMMcopy (*.cna.txt)
#cytoBand = {T, F}
#alphaVal = [0,1]
#geneAnnot is a dataframe with 4 columns: geneSymbol, chr, start, stop
#spacing is the distance between each track
plotCNlogRByChr <- function(dataIn, segs, param = NULL, logR.column = "logR", call.column = "event", plotSegs = TRUE, seqinfo=NULL, chr=NULL, ploidy = NULL, geneAnnot=NULL, yrange=c(-4,6), xlim=NULL, xaxt = "n", cex = 0.5, gene.cex = 0.5, plot.title = NULL, spacing=4, cytoBand=T, alphaVal=1, main){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  alphaSubcloneVal <- ceiling(alphaVal / 2 * 255); class(alphaVal) = "hexmode"
  cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 26))
  subcloneCol <- c("#00FF00")
  cnCol <- paste(cnCol,alphaVal,sep="")
  names(cnCol) <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
#  segCol <- cnCol
#  ## add in colors for subclone if param provided
#  if (!is.null(param)){
#    ind <- ((which.max(param$ct) + 1) : length(param$ct)) + 1
#    cnCol[ind] <- paste0(cnCol[ind], alphaSubcloneVal / 2)
#    segCol[ind] <- "#00FF00"
#  }
  # adjust for ploidy #
  if (!is.null(ploidy)){
    dataIn[, logR.column] <- as.numeric(dataIn[, logR.column]) + log2(ploidy / 2)
    
    if (!is.null(segs)){
      segs[, "median"] <- segs[, "median"] + log2(ploidy / 2)
    }
  }
  
  
  if (!is.null(chr)){
    for (i in chr){
      
      dataByChr <- dataIn[dataIn[,"chr"]==as.character(i),]
      
      #plot the data
      #if (outfile!=""){ pdf(outfile,width=10,height=6) }
      par(mar=c(spacing,8,4,2))
      #par(xpd=NA)
      coord <- (as.numeric(dataByChr[,"end"]) + as.numeric(dataByChr[,"start"]))/2
      if (is.null(xlim)){
        xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"start"]))
        xaxt <- "n"
      }
      if (is.null(plot.title)){
        plot.title <- paste("Chromosome ",i,sep="")
      }
      ## plot logR for bins ##
      plot(coord,as.numeric(dataByChr[, logR.column]),col=cnCol[dataByChr[, call.column]],
           pch=16, ylim=yrange,
           xlim=xlim, xaxt = xaxt, xlab="",ylab="Copy Number (log2 ratio)",
           cex.lab=1.5,cex.axis=1.5, cex=cex,las=1)
      title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)
      ## plot centre line ##
      lines(c(1,as.numeric(dataByChr[dim(dataByChr)[1],3])),rep(0,2),type="l",col="grey",lwd=0.75)
      if (!is.null(segs) & plotSegs){
        segsByChr <- segs[segs[,"chr"]==as.character(i),,drop=FALSE]
        ind <- segsByChr$subclone.status == FALSE
        apply(segsByChr[ind, ], 1, function(x){
          lines(x[c("start","end")], rep(x["median"], 2), col = cnCol[x[call.column]], lwd = 3)
          invisible()
        })
        if (sum(!ind) > 0){
          apply(segsByChr[!ind, ], 1, function(x){
            lines(x[c("start","end")], rep(x["median"], 2), col = subcloneCol, lwd = 3)
            invisible()
          })
        }
      }
      
      if (!is.null(geneAnnot)){
        #par(xpd=F)
        colnames(geneAnnot) <- c("Gene","Chr","Start","Stop")
        geneAnnot <- geneAnnot[geneAnnot[,"Chr"]==chr,]
        for (g in 1:dim(geneAnnot)[1]){
          print(geneAnnot[g,"Gene"])
          abline(v=as.numeric(geneAnnot[g,"Start"]),col="black",lty=3,xpd=F)
          abline(v=as.numeric(geneAnnot[g,"Stop"]),col="black",lty=3,xpd=F)			
          atP <- (as.numeric(geneAnnot[g,"Stop"]) - as.numeric(geneAnnot[g,"Start"]))/2 + as.numeric(geneAnnot[g,"Start"])
          if (atP < dataByChr[1,"start"]){ atP <- dataByChr[1,"start"] }
          else if (atP > dataByChr[dim(dataByChr)[1],"start"]){ atP <- dataByChr[dim(dataByChr)[1],"start"] }
          mtext(geneAnnot[g,"Gene"],side=3,line=0,at=atP,cex=gene.cex)
          
        }
      }
    }
  }else{  #plot for all chromosomes
    par(mar=c(spacing,8,2,2))
    #midpt <- (as.numeric(dataIn[,"end"]) + as.numeric(dataIn[,"start"]))/2
    #coord <- getGenomeWidePositions(dataIn[,"chr"],midpt)
    coord <- getGenomeWidePositions(dataIn[,"chr"],dataIn[,"end"], seqinfo)
    plot(coord$posns,as.numeric(dataIn[, logR.column]),
         col=cnCol[as.character(dataIn[,call.column])],pch=16,xaxt="n", ylim=yrange,
         xlim=c(1,as.numeric(coord$posns[length(coord$posns)])),
         xlab="",ylab="Copy Number (log2 ratio)",
         cex.lab=1.5,cex.axis=1.5,cex=0.5,las=1,bty="n",
         #main=dataIn[1,"sample"])
         main=main)
    #plot segments
    
    if (plotSegs){      
      coordEnd <- getGenomeWidePositions(segs[, "chr"], segs[, "end"], seqinfo)
      coordStart <- coordEnd$posns - (segs[, "end"] - segs[, "start"] + 1)
      xlim <- as.numeric(c(1, coordEnd$posns[length(coordEnd$posns)]))
      #col <- cnCol[as.numeric(segs[, "state"] + 1)]
      col <- cnCol[segs[, call.column]]
      #write.table(segs, "~/Documents/multisample/segs_debug.seg", quote=F, sep="\t", row.names=F)  ## debug
      value <- as.numeric(segs[, "median"])
      sc.status <- as.logical(segs[, "subclone.status"])
      mat <- as.data.frame(cbind(coordStart, coordEnd$posns, value, sc.status, col))
      rownames(mat) <- 1:nrow(mat)
      ## clonal CN
      ind <- mat$sc.status == FALSE
      apply(mat[ind, ], 1, function(x){
        lines(x[1:2], rep(x[3], 2), col = x[5], lwd = 3)
        invisible()
      })
      ## subclonal CN
      if (sum(!ind) > 0){
        apply(mat[!ind, ], 1, function(x){
          lines(x[1:2], rep(x[3], 2), col = subcloneCol, lwd = 3)
          invisible()
        })
      }
    }
    lines(as.numeric(c(1,coord$posns[length(coord$posns)])),rep(0,2),type="l",col="grey",lwd=2)
    if (plotSegs){
   	  plotChrLines(dataIn[,"chr"],coordEnd$chrBkpt,yrange)
   	}else{
   		plotChrLines(dataIn[,"chr"],coord$chrBkpt,yrange)
   	}
  }
}


plotCorrectionGenomeWide <- function(correctOutput, seqinfo = NULL, ...) {
  
  midpt <- (start(correctOutput) + end(correctOutput))/2
  coord <- getGenomeWidePositions(seqnames(correctOutput),midpt, seqinfo)
  
  
  par(mfrow = c(3, 1))
  #correctOutput <- correctOutput[paste(chr)]
  pos <- coord$posns
  from <- min(coord$chrBkpt)
  to <- max(coord$chrBkpt)
  
  copy <- correctOutput$reads / median(correctOutput$reads, na.rm = TRUE)
  top <- quantile(copy, 0.99)
  bot <- quantile(copy, 0.01)
  
  set <- which(correctOutput$valid & pos >= from & pos <= to &
                 copy >= bot & copy <= top)
  
  y <- copy[set]
  m <- signif(mad(y, na.rm = TRUE), digits = 3)
  r <- c(min(y, na.rm = TRUE), max(y, na.rm = TRUE))
  plot(pos[set], y, ylab = "Estimated Copy", xaxt="n",
       main = paste("Uncorrected Readcount, MAD = ", m), ylim = r, ...)
  plotChrLines(as.vector(unique(seqnames(correctOutput))),coord$chrBkpt,yrange=r+c(-0.5,0.5))
  m <- signif(mad(correctOutput$cor.gc[set], na.rm = TRUE), digits = 3)
  plot(pos[set], correctOutput$cor.gc[set], xaxt="n",
       ylab = "Estimated Copy",
       main = paste("CG-corrected Readcount, MAD = ", m), ylim = r, ...)
  plotChrLines(as.vector(unique(seqnames(correctOutput))),coord$chrBkpt,yrange=r+c(-0.5,0.5))
  m <- signif(mad(correctOutput$cor.map[set], na.rm = TRUE), digits = 3)
  plot(pos[set], correctOutput$cor.map[set], xaxt="n", ylab = "Estimated Copy",
       main = paste("Mappability and GC-corrected Readcount, MAD = ", m),
       ylim = r, ...)
  plotChrLines(as.vector(unique(seqnames(correctOutput))),coord$chrBkpt,yrange=r+c(-0.5,0.5))
}

##################################################
### HELPER FUNCTION TO GET GENOME-WIDE COORDS ####
##################################################
plotChrLines <- function(chrs,chrBkpt,yrange){
  #plot vertical chromosome lines
  for (j in 1:length(chrBkpt)){
    lines(rep(chrBkpt[j],2),yrange,type="l",lty=2,col="black",lwd=0.75)
  }
  numLines <- length(chrBkpt)
  mid <- (chrBkpt[1:(numLines-1)]+chrBkpt[2:numLines])/2
  #chrs <- mapSeqlevels(as.vector(chrs), style = "NCBI")
  if (seqlevelsStyle(chrs)[1] != "NCBI"){
  	seqlevelsStyle(chrs) <- "NCBI"
  }
  chrs <- sortSeqlevels(unique(chrs))
  axis(side=1,at=mid,labels=c(chrs),cex.axis=1.5,tick=FALSE)
}

# getGenomeWidePositions <- function(chrs,posns){  
#   #create genome coordinate scaffold
#   positions <- as.numeric(posns)
#   chrsNum <- unique(chrs)
#   chrBkpt <- rep(0,length(chrsNum)+1)
#   for (i in 2:length(chrsNum)){
#     chrInd <- which(chrs==chrsNum[i])
#     prevChrPos <- positions[chrInd[1]-1]      
#     chrBkpt[i] = prevChrPos
#     positions[chrInd] = positions[chrInd] + prevChrPos
#   }
#   chrBkpt[i+1] <- positions[length(positions)]
#   return(list(posns=positions,chrBkpt=chrBkpt))
# }

getGenomeWidePositions <- function(chrs, posns, seqinfo = NULL) {
    # create genome coordinate scaffold
    positions <- as.numeric(posns)
    chrs <- as.character(chrs)
    chrsNum <- unique(chrs)
    chrBkpt <- rep(0, length(chrsNum) + 1)
    prevChrPos <- 0
    i <- 1
    if (length(chrsNum) > 1){
		for (i in 2:length(chrsNum)) {
			chrInd <- which(chrs == chrsNum[i])
			if (!is.null(seqinfo)){
				prevChrPos <- seqlengths(seqinfo)[i-1] + prevChrPos
			}else{
				prevChrPos <- positions[chrInd[1] - 1]
			}
			chrBkpt[i] = prevChrPos
			positions[chrInd] = positions[chrInd] + prevChrPos
		}
	}
    chrBkpt[i + 1] <- positions[length(positions)]
    return(list(posns = positions, chrBkpt = chrBkpt))
}

#stateCols <- function() {
#  return(c("#74C476", "#238B45", "#00008B", "#A50F15", "#DE2D26", "#FB6A4A", "#FB6A4A", "#FB6A4A", "#FB6A4A", "#FB6A4A", "#FB6A4A", "#FB6A4A"))
#}

# plotSegments <- function(correctOutput, segmentOutput,
#                          chr = space(correctOutput)[1], ...){
#   if (is.null(segmentOutput$segs)) {
#     warning("Processed segments now found, automatically processing")
#     segmentOutput$segs <- processSegments(segments$segs,
#                                           space(correctOutput), start(correctOutput), end(correctOutput),
#                                           correctOutput$copy)
#   }
#   
#   segs <- segmentOutput$segs
#   correctOutput$state <- segmentOutput$state
#   cols <- stateCols()
#   range <- quantile(correctOutput$copy, na.rm = TRUE, prob = c(0.01, 0.99))
#   
#   a <- correctOutput[as.character(chr)]
#   b <- segs[segs$chr == chr, ]
#   plot(start(a), a$copy,
#        col = cols[as.numeric(as.character(a$state))], ylim = range, ...)
#   for (k in 1:nrow(b)){
#     lines(c(b$start[k], b$end[k]), rep(b$median[k], 2), lwd = 3,
#           col = "green")
#   }
# }

####  TODO: SUBCLONAL STATES #####
plotParam <- function(mus, lambdas, nu, subclone = NULL, copy.states = 0:6, ...) {
  #cols <- stateCols()
  cols <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 10))
  cols <- cols[copy.states + 1]
  #last <- ncol(mus)
  #lastmu <- mus[, last]
  #lastlambda <- lambdas[, last]
  
  domain <- (max(mus) - min(mus)) / 2
  left <- min(mus) - domain
  right <- max(mus) + domain
  x = seq(left * 10, right * 10, by = 0.01);
  height = 0
  
  for(state in 1:length(mus)) {
    #y = tdistPDF(x, param$mu[state,1], param$lambda[state,1], param$param$nu);
    z = tdistPDF(x, mus[state], lambdas[state], nu);
    height <- max(height, z)
  }
  
  plot(c(left, right), c(0, height), type = "n", xlab = "Normalized Log2 Ratios",
       ylab = "Student-t Density", ...);
  
  for(state in 1:length(mus)) {
    #y = tdistPDF(x, param$mu[state,1], param$lambda[state,1], param$param$nu);
    z = tdistPDF(x, mus[state], lambdas[state], nu);
    #lines(x, y, col = cols[state], lwd = 2, lty = 3);
    lines(x, z, col = cols[state], lwd = 3);
  }
}
