# file:   plotting.R
# author: Gavin Ha, Ph.D.
#         Fred Hutchinson Cancer Research Center
# contact: <gha@fredhutch.org>
# # website: https://GavinHaLab.org
#
# author: Justin Rhoades, Broad Institute
#
# ichorCNA website: https://github.com/broadinstitute/ichorCNA
# date:   January 6, 2020
# description: Hidden Markov model (HMM) to analyze Ultra-low pass whole genome sequencing (ULP-WGS) data.
# This file contains R functions for plotting.

## plot solutions for all samples
plotSolutions <- function(hmmResults.cor, tumour_copy, chrs, outDir, counts,
                          logR.column = "logR", call.column = "event", likModel = "t",
                          plotSegs = TRUE, numSamples=1, plotFileType="pdf", 
                          plotYLim=c(-2,2), seqinfo = NULL,
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

    ### PLOT THE GENOME-WIDE CORRECTION COMPARISONS ###
    outPlotFile <- paste0(outDir,"/",id,"/",id,"_genomeWideCorrection")
    if (plotFileType == "png"){ 
      outPlotFile <- paste0(outPlotFile, ".png")
      png(outPlotFile,width=10,height=12,units="in",res=300)
    }else{
      outPlotFile <- paste0(outPlotFile, ".pdf")
      pdf(outPlotFile,width=10,height=12)
    }
    plotCorrectionGenomeWide(tumour_copy[[s]], seqinfo = seqinfo, pch = ".", xlab = "Chromosomes")
    dev.off()
    
    ### PLOT THE CORRECTION COMPARISONS BY CHR ###
    for (i in chrs){
      outPlotFile <- paste0(outDir,"/",id,"/",id,"_correction_chr", i)
      if (plotFileType == "png"){ 
        outPlotFile <- paste0(outPlotFile, ".png")
        png(outPlotFile,width=10,height=12,units="in",res=300)
      }else{
        outPlotFile <- paste0(outPlotFile, ".pdf")
        pdf(outPlotFile,width=10,height=12)
      }
      plotCorrectionGenomeWide(tumour_copy[[s]], seqinfo = seqinfo, chr = i, 
                               cex = 3, pch = ".", xlab = paste0("Chr", i))
      dev.off()
    }
    
    ### PLOT FIT COMPARISON ### only if replication correction was performed
    if (!is.null(counts[[s]]$counts$cor.rep) && length(counts[[s]]$counts) > 1e4){
      for (i in chrs){
        outPlotFile <- paste0(outDir, "/", id, "/", id, "_fit_chr", i)
        if (plotFileType == "png"){ 
          outPlotFile <- paste0(outPlotFile, ".png")
          png(outPlotFile,width=10,height=12,units="in",res=300)
        }else{
          outPlotFile <- paste0(outPlotFile, ".pdf")
          pdf(outPlotFile,width=10,height=12)
        }
        plotFitCompareByChr(counts[[s]], chr = i, covar = "repTime", covarName = "Replication Timing",
                                        before = "cor.map", beforeName = "Mappability-Corrected",
                                        after = "cor.rep", afterName = "Replication-Timing-Corrected")
        dev.off()
      }
    }

    ### PLOT THE BIAS ###
    outPlotFile <- paste0(outDir,"/",id,"/",id,"_bias")
    if (plotFileType == "png"){ 
      outPlotFile <- paste0(outPlotFile, ".png")
      png(outPlotFile,width=7,height=7,units="in",res=300)
    }else{
      outPlotFile <- paste0(outPlotFile, ".pdf")
      pdf(outPlotFile,width=7,height=7)
    }
    par(mfrow = c(3, 2))
    try(plotCovarBias(counts[[s]], covar = "gc", before = "reads", 
                    after = "cor.gc", fit = "gc.fit", xlab = "GC Content",
                    pch = 20, cex = 0.5, mfrow = NULL),
    silent = TRUE)
    try(plotCovarBias(counts[[s]], covar = "map", before = "cor.gc", 
                      after = "cor.map", fit = NULL, xlab = "Mappability Score",
                      pch = 20, cex = 0.5, mfrow = NULL),
    silent = TRUE)
    try(plotCovarBias(counts[[s]], covar = "repTime", before = "cor.map", 
                      after = "cor.rep", fit = "rep.fit", xlab = "Replication Timing", 
                      pch = 20, cex = 0.5, mfrow = NULL),
    silent = TRUE)
    dev.off()

    ### PLOT TPDF ##
    outPlotFile <- paste0(outDir,"/",id,"/",id,"_tpdf.pdf")
    pdf(outPlotFile)
    plotParam(mus = unique(hmmResults.cor$results$mus[, s, iter]), 
              lambdas = hmmResults.cor$results$lambdas[, s, iter], 
              vars = hmmResults.cor$results$vars[, s],
              jointStates = hmmResults.cor$results$param$jointStates[, s],
              likModel = likModel,
              subclone = hmmResults.cor$results$param$ct.sc.status,
              nu = hmmResults.cor$results$param$nu, copy.states = hmmResults.cor$results$param$ct)
    dev.off()
  }
}


plotGWSolution <- function(hmmResults.cor, s, outPlotFile, plotFileType="pdf", 
						   logR.column = "logR", call.column = "event",
						   seqinfo = NULL, plotSegs = TRUE, 
               cex = 0.5, cex.axis = 1.5, cex.lab=1.5, cex.text=1.5,
               plotYLim=c(-2,2), estimateScPrevalence, main=NULL, spacing=4,
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
    if (par()$mfrow[1] > 1){ # more than one plot in figure
      lines <- c(0.25, -1.5, -3.25)
    }else{       
      lines <- c(0.25, -1, -2.25)
    }
    #par(oma=c(0, 0, 2, 0))
    if (plotSegs){
    	segsToUse <- hmmResults.cor$results$segs[[s]]
    }else{
    	segsToUse <- NULL
    }
    plotCNlogRByChr(hmmResults.cor$cna[[s]], segs = segsToUse, plotSegs=plotSegs, seqinfo=seqinfo,
                    param = hmmResults.cor$results$param, chr=NULL,
                    logR.column = logR.column, call.column = call.column,  
                    cex = cex, cex.axis = cex.axis, cex.lab=cex.lab,
                    ploidy = ploidyAll, cytoBand=T, yrange=plotYLim, spacing=spacing, main=NULL)  #ylim for plot
    annotStr <- paste0("Tumor Fraction: ", signif(purityEst, digits=4), ", Ploidy: ", signif(ploidyEst, digits=3))
    if (!is.null(coverage)){
      annotStr <- paste0(annotStr, ", Coverage: ", signif(coverage, digits=2))
    }
    mtext(line=lines[1], main, cex=cex.text)
    mtext(line=lines[2], annotStr, cex=cex.text)
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
        mtext(line=lines[3], annotSubStr, cex=cex.text)
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
plotCNlogRByChr <- function(dataIn, segs, param = NULL, logR.column = "logR", 
  call.column = "event", plotSegs = TRUE, seqinfo=NULL, chr=NULL, ploidy = NULL, 
  geneAnnot=NULL, yrange=c(-4,6), xlim=NULL, xaxt = "n", cex = 0.5, cex.axis = 1.5, cex.lab=1.5, gene.cex = 0.5, 
  plot.title = NULL, spacing=4, cytoBand=T, alphaVal=1, main){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  alphaSubcloneVal <- ceiling(alphaVal / 2 * 255); class(alphaVal) = "hexmode"
  cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 1001))
  subcloneCol <- c("#00FF00")
  cnCol <- paste(cnCol,alphaVal,sep="")
  names(cnCol) <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0("HLAMP", 2:1000))
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
      par(mar=c(spacing,8,spacing,2))
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
           cex.lab=cex.lab,cex.axis=cex.axis, cex=cex,las=1)
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
         cex.lab=cex.lab,cex.axis=cex.axis,cex.main=1.5,cex=cex,las=1,bty="n",
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


plotCorrectionGenomeWide <- function(correctOutput, seqinfo = NULL, chr = NULL, ...) {
  
  if (is.null(chr)){
    midpt <- (start(correctOutput) + end(correctOutput))/2
    coord <- getGenomeWidePositions(seqnames(correctOutput),midpt, seqinfo)
  }else {
    correctOutput <- correctOutput[seqnames(correctOutput) == chr]
    midpt <- (start(correctOutput) + end(correctOutput))/2
    coord <- getGenomeWidePositions(seqnames(correctOutput),midpt, seqinfo)
  }
  
  mapCor <- !is.null(correctOutput$cor.map)
  repTimeCor <- !is.null(correctOutput$cor.rep)
  rows <- 2 + sum(c(mapCor, repTimeCor)) # must include uncorrected and GC panels
  par(mfrow = c(rows, 1))
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
  # plot Uncorrected read ratios (tum:norm)
  m <- signif(mad(diff(y), na.rm = TRUE), digits = 3)
  r <- c(min(y, na.rm = TRUE), max(y, na.rm = TRUE))
  plot(pos[set], y, ylab = "Estimated Copy Ratio", xaxt="n",
       main = paste("Uncorrected Readcount, MAD = ", m), ylim = r, ...)
  if (is.null(chr)){
    plotChrLines(as.vector(unique(seqnames(correctOutput))),coord$chrBkpt,yrange=r+c(-0.5,0.5))
  }
  # plot GC corrected read counts
  m <- signif(mad(diff(correctOutput$cor.gc[set]), na.rm = TRUE), digits = 3)
  plot(pos[set], correctOutput$cor.gc[set], xaxt="n",
       ylab = "Estimated Copy Ratio",
       main = paste("GC-corrected Readcount, MAD = ", m), ylim = r, ...)
  if (is.null(chr)){
    plotChrLines(as.vector(unique(seqnames(correctOutput))),coord$chrBkpt,yrange=r+c(-0.5,0.5))
  }
  # plot GC corrected + mappability corrected read counts
  if (mapCor){
    m <- signif(mad(diff(correctOutput$cor.map[set]), na.rm = TRUE), digits = 3)
    plot(pos[set], correctOutput$cor.map[set], xaxt="n", ylab = "Estimated Copy Ratio",
         main = paste("GC-corrected, mappability-corrected Readcount, MAD = ", m),
         ylim = r, ...)
    if (is.null(chr)){
      plotChrLines(as.vector(unique(seqnames(correctOutput))),coord$chrBkpt,yrange=r+c(-0.5,0.5))
    }
  }
  # plot GC corrected + mappability corrected + replication-timing corrected read counts
  if (repTimeCor){
    m <- signif(mad(diff(correctOutput$cor.rep[set]), na.rm = TRUE), digits = 3)
    plot(pos[set], correctOutput$cor.rep[set], xaxt="n", ylab = "Estimated Copy Ratio",
         main = paste("GC-corrected, mappability-corrected, replication-timing-corrected Readcount, MAD = ", m),
         ylim = r, ...)
    if (is.null(chr)){
      plotChrLines(as.vector(unique(seqnames(correctOutput))),coord$chrBkpt,yrange=r+c(-0.5,0.5))
    }
  }
}

##################################################
### HELPER FUNCTION TO GET GENOME-WIDE COORDS ####
##################################################
plotChrLines <- function(chrs,chrBkpt,yrange, cex.axis=1.5){
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
  axis(side=1,at=mid,labels=c(chrs),cex.axis=cex.axis,tick=FALSE)
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

plotParam <- function(mus, lambdas, nu, vars = NULL, likModel = "t", 
                      jointStates = NULL, subclone = NULL, copy.states = 0:6, ...) {
  #cols <- stateCols()
  cols <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 100))
  cols <- cols[copy.states + 1]
  
  domain <- (max(mus) - min(mus)) / 2
  left <- min(mus) - domain
  right <- max(mus) + domain
  x = seq(min(-0.2, left * 2), max(0.2, right * 2), by = 0.01);
  height = 0
  K <- length(mus)
  z <- vector('list', K)
  for(state in 1:K) {
    if (likModel == "t"){
      z[[state]] <- tdistPDF(x, mus[state], lambdas[state], nu);
      ylab <- "Student-t Density"
    }else if (grepl("gauss", likModel, ignore.case = TRUE)){
      if (is.null(var) || is.null(jointStates)){
        stop("plotParam: likModel is Gaussian but var and/or jointStates are not provided.")
      }
      varsToUse <- mean(vars[jointStates == state])
      z[[state]] <- normalpdf(x, mus[state], varsToUse)
      ylab <- "Gaussian Density"
    }
  }
  height <- max(unlist(z))
  
  plot(c(min(x), max(x)), c(0, height), type = "n", xlab = "Normalized Log2 Ratios",
       ylab = ylab, ...)
  
  for (state in 1:K){
    if (subclone[state]){ # subclonal
      lines(x, z[[state]], col = cols[state], lwd = 3, lty = 2)
    }else{ #not subclonal
      lines(x, z[[state]], col = cols[state], lwd = 3)
    }
  }
  
}

plotCovarBias <- function(correctOutput, covar = "gc", 
                          before = "reads", after = "cor.gc", fit = "gc.fit", 
                          points = 10000, xlab = "GC content", 
                          mfrow = c(1,2), ...){
  if (!is.null(mfrow)){
    par(mfrow = mfrow)
  }
  coutlier = 0.001
  counts <- as.data.frame(correctOutput$counts)
  
  # select points to show (before)
  set <- which(counts$ideal)
  #range <- quantile(counts[[before]][counts$ideal], 
  #                  prob = c(0, 1 - coutlier), na.rm = TRUE)
  #valid <- which(counts[[before]] >= range[1] & counts[[before]] <= range[2])
  #set <- intersect(counts$ideal, valid)
  #select.1 <- sample(valid, min(length(valid), points))
  select.1 <- sample(set, min(length(set), points))
  # plot before correction
  plot(counts[[covar]][select.1], counts[[before]][select.1], 
       col = densCols(counts[[covar]][select.1], counts[[before]][select.1]), 
       ylab = "Uncorrected", xlab = xlab, 
       ...)
  # plot curve fit line
  if (!is.null(fit)){ # want to look at fit
    if (!is.null(correctOutput[[fit]])){ # actually have the fit object/model as list element
      domain <- c(min(counts[[covar]], na.rm = TRUE), max(counts[[covar]], na.rm = TRUE))
      fit.covar <- correctOutput[[fit]]
      i <- seq(domain[1], domain[2], by = 0.001)
      y <- predict(fit.covar, i)
      lines(i, y, type="l", lwd=1, col="red")
    }
  }
  
  # select points to show (after)
  range <- quantile(counts[[after]][counts$ideal], 
                    prob = c(0, 1 - coutlier), na.rm = TRUE)
  valid <- which(counts[[after]] >= range[1] & counts[[after]] <= range[2])
  select.2 <- intersect(valid, select.1)
  # plot after correction
  plot(counts[[covar]][select.2], counts[[after]][select.2], 
       col = densCols(counts[[covar]][select.2], counts[[after]][select.2]), 
       ylab = "Corrected", xlab = xlab, 
       ...)
}

######### INCOMPLETE ############
## y is the character string for the column to fit the data
## x is the GRanges object
data.fit <- function(x, y){
  ind <- !is.na(values(x)[[y]])
  x <- x[ind]
  midpt <- start(ranges(x)) + (end(ranges(x)) - start(ranges(x)))
  fit <- loess(values(x)[[y]] ~ midpt, span = 0.03)
  fit.predict <- loess(predict(fit, midpt) ~ midpt, span = 0.3)
  y.hat <- predict(fit, midpt)
  return(list(fit = fit, y.hat = y.hat, midpt = midpt))
}

plotFitCompareByChr <- function(x, chr, covar = "repTime", covarName = "Replication Timing",
                                before = "cor.map", beforeName = "Mappability-Corrected",
                                after = "cor.rep", afterName = "Replication-Timing-Corrected"){
  counts.chr <- x$counts[seqnames(x$counts) == chr]
  #counts.chr <- counts[[1]]$counts[seqnames(counts[[1]]$counts) == chr]
  
  fit.map <- data.fit(counts.chr, y = before)
  midpt <- fit.map$midpt
  domain <- c(min(midpt, na.rm = TRUE), max(midpt, na.rm = TRUE))
  xlim <- c(domain[1], domain[2])
  ylim <- c(min(c(counts.chr$cor.map, counts.chr$repTime), na.rm = TRUE),
            max(c(counts.chr$cor.map, counts.chr$repTime), na.rm = TRUE))
  
  par(mfrow = c(6,1), mar = c(4,6,1,3))
  
  # before rep time, after mappability
  plot(midpt, fit.map$y.hat, type="l", lwd=1, col="red", xlab = paste0("Chromosome ", chr), 
       ylab = paste0(beforeName, "\nTumor Read Counts"))
  # after rep time
  fit.rep <- data.fit(counts.chr, y = after)
  plot(fit.rep$midpt, fit.rep$y.hat, type="l", lwd=1, col="red", xlab = paste0("Chromosome ", chr), 
       ylab = paste0(afterName, "\nTumor Read Counts"))
  if (length(values(counts.chr)[[paste0(before, ".normal")]]) > 0){
    fit.map.normal <- data.fit(counts.chr, y = paste0(before, ".normal"))
    plot(fit.map.normal$midpt, fit.map.normal$y.hat, type="l", lwd=1, col="red", xlab = paste0("Chromosome ", chr), 
         ylab = paste0(beforeName, "\nNormal Read Counts"))
    fit.rep.normal <- data.fit(counts.chr, y = paste0(after, ".normal"))
    plot(fit.rep.normal$midpt, fit.rep.normal$y.hat, type="l", lwd=1, col="red", xlab = paste0("Chromosome ", chr), 
         ylab = paste0(afterName, "\nNormal Read Counts"))
    fit.copy <- data.fit(counts.chr, y = "copy")
    plot(fit.copy$midpt, 2^fit.copy$y.hat, type="l", lwd=1, col="blue", xlab = paste0("Chromosome ", chr), 
         ylab = "Tumour:Normal Read Counts")
  }
  
  # covariance (e.g. gc, repTime)
  fit.covar <- data.fit(counts.chr, y = covar)
  covar.midpt <- start(ranges(counts.chr)) + (end(ranges(counts.chr)) - start(ranges(counts.chr)))
  plot(covar.midpt, values(counts.chr)[[covar]], type = "l", xlab = paste0("Chromosome ", chr), 
       ylab = covarName)
  
}
