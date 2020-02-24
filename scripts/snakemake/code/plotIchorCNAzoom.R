#' plotIchorCNAzoom.R
#' author: Gavin Ha 
#' Fred Hutchinson Cancer Research Center
#' contact: <gha@fredhutch.org>
#' date:  October 8, 2019
#' description: Generate plots of copy number from On- and Off-target ichorCNA results


library(optparse)
option_list <- list(
  make_option(c("--id"), type = "character", help = "Sample ID"),
  #make_option(c("--plot_funcs"), type = "character", help = "Path to file containing plotting R functions to source."),
  make_option(c("--cnFile"), type="character", help = "Path to ichorCNA cna.seg output file."),
  make_option(c("--paramFile"), type="character", help = "Path to ichorCNA params.txt output file."),
  make_option(c("--geneFile"), type="character", default=NULL, help = "Path to file containing list of genes with chr, start, end coordinates."),
  make_option(c("--genomeBuild"), type="character", default="hg19", help = "Genome build: hg19 or hg38. Default [%default]"),
  make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
  make_option(c("--zoom"), type = "logical", default = FALSE, help = "Zoom plot; if TRUE, then requires --chrs --start --end to be set. [Default: %default]"),
  make_option(c("--yaxis"), type = "character", default = "integer", help = "Data type to plot for y-axis (\"integer\" copy number or \"logratio\"). [Default: %default]"),
  make_option(c("--chrs"), type = "character", default = "c(1:22, 'X')", help = "Chromosomes to plot; string [Default: %default"),
  make_option(c("--start"), type = "integer", default = NULL, help = "Start coordinate for zoom plots"),
  make_option(c("--end"), type = "integer", default = NULL, help = "End coordinate for zoom plots"),
  make_option(c("--plotYlim"), type = "character", default = "c(-2,2)", help = "Y limits for plotting log ratio. [Default: %default]."),
  make_option(c("--plotSize"), type = "character", default = "c(5,3)", help = "width and height in inches. [Default: %default]."),
 # make_option(c("--plotFormat"), type = "character", default = "png", help = "File format of plot. E.g. pdf or png. [Default: %default]."),
  make_option(c("--outPlotFile"), type="character", help="Path to output figure file.")
 # make_option(c("--outDir"), type="character", help="Path to output directory.")
)
parseobj <- OptionParser(option_list=option_list, usage = "usage: Rscript %prog [options]")
opt <- parse_args(parseobj)
print(opt)

library(data.table)
library(GenomicRanges)
library(stringr)
#library(ggplot2)
#library(reshape2)
#library(diagram)
#library(igraph)
#library(tools)
#library(SNPchip)
#library(foreach)
#library(VariantAnnotation)
#library(doMC)
#library(quantsmooth)

#source(opt$plot_funcs)

args <- commandArgs(TRUE)
options(stringsAsFactors=F, scipen=999, bitmapType = "cairo", width=175, useDingbats = FALSE)

id <- opt$id
cnFile <- opt$cnFile
paramFile <- opt$paramFile
chrStr <- as.character(eval(parse(text = "c(opt$chrs)")))
startPos <- opt$start
endPos <- opt$end
zoom <- opt$zoom
ylim <- eval(parse(text = opt$plotYlim))
geneList <- opt$geneFile
outPlotFile <- opt$outPlotFile
plotFormat <- tools::file_ext(outPlotFile)
#outDir <- opt$outDir
outImage <- gsub(plotFormat, "RData", outPlotFile)
genomeBuild <- opt$genomeBuild
genomeStyle <- opt$genomeStyle
yaxis <- opt$yaxis

plotSize <- eval(parse(text=opt$plotSize))
width <- plotSize[1]  #6 8 
height <- plotSize[2]  #3 3.5 #4 
spacing <- 3

loess.span <- 0.1
numCores <- 6
coverageThres <- 0.1
minPurity <- 0.10
diffLogRThres <- 0.5
madThres <- 0.2
numPlotCols <- 3
ylimMin <- ylim[1]
ylimSV <- ylim
cex.txt <- 0.75
plotSegs <- FALSE
buffer <- 1e4
offset.factor <- 1.15 # sv drawn outside of plot
pch <- 20
if (zoom){
  startTitle <- paste0(format(round(startPos/1e6,2), nsmall=2))
  endTitle <- paste0(format(round(endPos/1e6,2),nsmall=2), "Mb")
  xlim <- c(startPos, endPos)
  cex <- 0.75
  showCytoBand <- FALSE
  xaxt <- "s"
  plotAtCentre <- FALSE
  cnColor <- FALSE
  plotIdio <- FALSE
}else{
  xlim <- NULL
  cex <- 0.25
  showCytoBand <- TRUE
  xaxt <- "n"
  #yaxis <- "logratio"
  plotAtCentre <- FALSE
  cnColor <- FALSE
  plotIdio <- TRUE
  plotSegs <- FALSE
  spacing <- 6
}
if (plotFormat == "png"){
  width <- width * 300
  height <- height * 300
  cex <- 1
  pch <- 19
  cnColor <- TRUE
  cex.txt <- 1.5
}
if (!cnColor){
  #cnCol <- rep("#000000", 30)
  cnCol <- rep("black", 30)
}else{
  cnCol <- NULL
}
if (chrStr == "0" || is.null(chrStr) || chrStr == "None"){
  chrStr <- as.character(c(1:22, "X"))
}
bsg <- paste0("BSgenome.Hsapiens.UCSC.", genomeBuild)
if (!require(bsg, character.only=TRUE, quietly=TRUE, warn.conflicts=FALSE)) {
  seqinfo <- Seqinfo(genome=genomeBuild)
} else {
  seqinfo <- seqinfo(get(bsg))
}
seqlevelsStyle(seqinfo) <- genomeStyle
seqlevelsStyle(chrStr) <- genomeStyle
seqinfo <- keepSeqlevels(seqinfo, value = chrStr)


if (!is.null(geneList) && geneList != "None"){
  genes <- read.delim(geneList, header=F, as.is=T)
  seqlevelsStyle(genes$V2) <- genomeStyle
}else{
  genes <- NULL
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

message("Analyzing ", id)
ulp <- fread(cnFile)
colnames(ulp) <- c("Chr", "Start", "End", "copy.number", "event", "LogRatio", "subclone.status", "Corrected_Copy_Number", "Corrected_Call", "logR_Copy_Number")
params <- read.delim(paramFile, header=T, as.is=T)
purity <- as.numeric(params[1, 2]) 
ploidyT <- as.numeric(params[1, 3])
normCN <- 2
ploidyS <- purity * ploidyT + (1-purity) * normCN
if (yaxis == "integer"){
	#ulp[Chr!="X", LogRatio := log2(logRbasedCN(LogRatio, purity, ploidyT, cn=2))]
	#ulp[Chr=="X", LogRatio := log2(logRbasedCN(LogRatio, purity, ploidyT, cn=1))]
	colName <- "logR_Copy_Number"
}else{
  ulp[, LogRatio := LogRatio + log2(ploidyS / 2)]
	#segs$LogRatio <- segs$Median_logR
	#segs$LogRatio <- segs$LogRatio + log2(ploidyS / 2)
	colName <- "LogRatio"
}

##################################################
########## FUNCTION TO PLOT ZOOM RESULTS #########
##################################################
## format of dataIn is output from /home/unix/gavinha/software/code/git/scripts/titan/analysis/combineTITAN-ichor.R
plotTitanIchorCNA <- function(dataIn, param = NULL, colName = "LogRatio", callColName="Corrected_Call", segs=NULL, chr=NULL, purity = NULL, ploidyT = NULL, geneAnnot=NULL, yrange=c(-4,6), yaxis = "logRatio", xlim=NULL, xaxt = "n", cex = 0.5, gene.cex = 0.5, pch = 16, plot.title = NULL, cnCol = NULL, spacing=4, cytoBand=T, alphaVal=1, plot.new = TRUE, ...){
  #color coding
  alphaVal <- ceiling(alphaVal * 255); class(alphaVal) = "hexmode"
  alphaSubcloneVal <- ceiling(alphaVal / 2 * 255); class(alphaVal) = "hexmode"
  subcloneCol <- c("#00FF00")
  if (is.null(cnCol)){
    cnCol <- c("#00FF00","#006400","#0000FF","#8B0000",rep("#FF0000", 26))
    cnCol <- c(cnCol, "HET"="#0000FF", "DLOH"="#006400", "NLOH"="#0000FF", "ALOH"="#FF0000", "ASCNA"="#FF0000", "BCNA"="#FF0000", "UBCNA"="#FF0000")
  }else{
    cnCol.col <- as.character(cnCol[1])
    cnCol <- c(cnCol, "HET"=cnCol.col, "DLOH"=cnCol.col, "NLOH"=cnCol.col, "ALOH"=cnCol.col, "ASCNA"=cnCol.col, "BCNA"=cnCol.col, "UBCNA"=cnCol.col)
  }
  names(cnCol)[1:30] <- c("HOMD","HETD","NEUT","GAIN","AMP","HLAMP",paste0(rep("HLAMP", 8), 2:25))
  #cnCol <- paste(cnCol,alphaVal,sep="")
  # adjust for ploidy #
  normCN <- 2
  if (!is.null(ploidyT) & yaxis != "integer"){
    ploidyS <- purity * ploidyT + (1-purity) * normCN
    dataIn[, colName] <- as.numeric(dataIn[, colName]) + log2(ploidyS / 2)
    
    if (!is.null(segs)){
      segs[, colName] <- segs[, colName] + log2(ploidyS / 2)
    }
  }
  
  
  if (!is.null(chr)){
    for (i in chr){
      dataByChr <- dataIn[dataIn[,"Chr"]==as.character(i),]
       ## set y axis labels as either integer or logR copy number
      #avgTumPloidy <- round(ploidyT)
 
      zero <- 0.5  
      cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      #ploidyToUse <- ploidyS
      if (i == "X"){
        normCN <- 1
        zero <- 0.25
        cn <- c(0, 1, 2, 4, `^`(2, 3:(yrange[2]+1)))
      }      
      if (yaxis == "integer"){
        y.ticks <- log2(cn)
        y.ticks[1] <- log2(zero)  
        yrange[1] <- y.ticks[1]    
        ylab <- "Copy Number"
        #dataByChr[, colName] <- log2(logRbasedCN(dataByChr[, colName], purity, ploidyT, cn=normCN))
        dataByChr[, colName] <- log2(dataByChr[, colName])
        if (!is.null(segs)){
          segs[, colName] <- log2(segs[, colName])# + log2(ploidyS / 2)
        }
        centreLine <- log2(normCN)
      }else{      
        #dataByChr[, colName] <- dataByChr[, colName] + log2(ploidyS / 2)
        cnLog <- log2(cn[-which(cn==3)] / normCN)  
        cn <- seq(-2,yrange[2],2)#c(-2, cn)
        y.ticks <- cn
        ylab <- "Copy Number (log2 ratio)"
        centreLine <- 0
      }

      #plot the data
      #if (outfile!=""){ pdf(outfile,width=10,height=6) }
      par(mar=c(spacing,8,4,2))
      #par(xpd=NA)
      coord <- (as.numeric(dataByChr[,"End"]) + as.numeric(dataByChr[,"Start"]))/2
      if (is.null(xlim)){
        xlim <- c(1,as.numeric(dataByChr[dim(dataByChr)[1],"Start"]))
        xaxt <- "n"
      }
      if (is.null(plot.title)){
        plot.title <- paste("Chromosome ",i,sep="")
      }
      ## plot logR for bins ##
      if (plot.new){
        plot(coord,as.numeric(dataByChr[, colName]),col=cnCol[dataByChr[,callColName]],
             pch=pch, ylim=yrange, yaxt="n",
             xlim=xlim, xaxt = xaxt, xlab="",ylab=ylab,
             cex.lab=1.5,cex.axis=1.5, cex=cex,las=1, ...)
      }else{
        points(coord,as.numeric(dataByChr[, colName]),col=cnCol[dataByChr[,callColName]],
             pch=pch, cex=cex, ...)
      }
      axis(2, at=y.ticks, labels=cn, las=2, cex.axis=1.5)
      title(plot.title, line = 1.25, xpd=NA, cex.main=1.5)
      ## plot centre line ##
      lines(c(1,tail(na.omit(dataByChr[,3]), 1)),rep(centreLine,2),type="l",col="grey",lwd=0.75)
      if (!is.null(segs)){
        segsByChr <- segs[segs[,"Chromosome"]==as.character(i),,drop=FALSE]
        #ind <- segsByChr$subclone.status == FALSE
        apply(segsByChr, 1, function(x){
          lines(x[c("Start","End")], rep(x[colName], 2), col = cnCol[x[callColName]], lwd = 3)
          invisible()
        })
        #if (sum(!ind) > 0){
        #  apply(segsByChr[!ind, ], 1, function(x){
        #    lines(x[c("Start","End")], rep(x["Median_logR"], 2), col = subcloneCol, lwd = 3)
        #    invisible()
        #  })
        #}
      }
      
      if (cytoBand==TRUE){
        #require(quantsmooth)
        #par(xpd = NA)
        #paintCytobands(chrom=chr, units="bases", pos=c(0,(yrange[1]-0.1)), width=0.75, legend=F) 
      }
      
      if (!is.null(geneAnnot)){
        #par(xpd=F)
        colnames(geneAnnot) <- c("Gene","Chr","Start","Stop")
        geneAnnot <- geneAnnot[geneAnnot[,"Chr"]==chr,]
        if (nrow(geneAnnot) > 0){
        for (g in 1:dim(geneAnnot)[1]){
          print(geneAnnot[g,"Gene"])
          abline(v=as.numeric(geneAnnot[g,"Start"]),col="black",lty=3,xpd=F)
          abline(v=as.numeric(geneAnnot[g,"Stop"]),col="black",lty=3,xpd=F)     
          atP <- (as.numeric(geneAnnot[g,"Stop"]) - as.numeric(geneAnnot[g,"Start"]))/2 + as.numeric(geneAnnot[g,"Start"])
          if (atP < dataByChr[1,"Start"]){ atP <- dataByChr[1,"Start"] }
          else if (atP > dataByChr[dim(dataByChr)[1],"Start"]){ atP <- dataByChr[dim(dataByChr)[1],"Start"] }
          mtext(geneAnnot[g,"Gene"],side=3,line=0,at=atP,cex=gene.cex)
          }
        }
      }
    }
  }
}


#####################################
########## PLOT CHR RESULTS #########
#####################################	 
for (j in 1:length(chrStr)){
  ###################################
  if (genomeStyle == "NCBI"){
    chrTitle <- paste0("chr", chrStr[j])
  }else{
    chrTitle <- chrStr[j]
  }
  plotTitle <- paste0(id, " (", chrTitle,")")
  if (zoom){
    ylimMax <- ulp[Chr==chrStr[j] & Start >= xlim[1] & Start <= xlim[2], max(logR_Copy_Number, na.rm=T)] + 1
    #outPlot <- paste0(outPlotDir, "/", id, "_CNA-SV-BX_",plotType,"_chr",chrStr[j],"-",startPos,"-",endPos,".pdf")
    plotTitle <- paste0(id, " (",chrTitle,":",startTitle,"-",endTitle, ")")
  }else{
    xlim <- c(1, seqlengths(seqinfo)[chrStr[j]])
    ylimMax <- ulp[, max(LogRatio, na.rm=T)] + 1
  }    

  ylim[2] <- min(max(ylim[2], ceiling(ylimMax)), 10)
  
  if (plotFormat == "png"){
    png(outPlotFile, width = width, height=height)
  }else{
    pdf(outPlotFile, width = width, height=height)
  }
	
  if (plotSegs) { segsToPlot <- segs } else { segsToPlot <- NULL}
  message("Plotting read depth CN")
  plotTitanIchorCNA(as.data.frame(ulp), segs=segsToPlot, chr=chrStr[j], colName=colName, 
      cytoBand=showCytoBand, geneAnnot=genes, purity = purity, ploidyT = NULL, yaxis=yaxis, cnCol = cnCol,
      yrange=ylim, xlim=xlim, spacing=spacing, xaxt=xaxt, cex = cex, gene.cex = 1, 
      plot.title = plotTitle, pch = pch, bg="black",)

  mtext(text = paste0("Tumor Fraction: ", format(purity, digits=2)), 
            side=1, line=-2, at = xlim[2], padj=0, adj=1, cex=cex.txt)
  mtext(text = paste0("Tumor Ploidy: ", format(ploidyT, digits=3)), 
            side=1, line=-1, at = xlim[2], padj=0, adj=1, cex=cex.txt)

  dev.off()
		
}

#save.image(file=outImage)
