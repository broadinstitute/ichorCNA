library(HMMcopy)
library(GenomicRanges)
library(optparse)
library(ggplot2)

options(stringsAsFactors=FALSE, scipen=0)
options(bitmapType='cairo')

option_list <- list(
	make_option(c("--gcWig"), type = "character", help = "GC Wig file for reference genome"),
	make_option(c("--mapWig"), type = "character", default=NULL, help = "Mappabiliy Wig file for reference genome"),
	make_option(c("--repTimeWig"), type = "character", default=NULL, help ="Path to replication timing WIG file. Default: [%default]"),
	make_option(c("-f", "--filelist"), type = "character", help = "List of of wig files."),
	make_option(c("-o", "--outfile"), type = "character", help = "Output file."),
	make_option(c("-c", "--centromere"), type="character", help = "File containing Centromere locations"),
	make_option(c("--rmCentromereFlankLength"), type="numeric", default=1e5, help="Length of region flanking centromere to remove. Default: [%default]"),
	make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze."),
  	make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--genomeBuild"), type="character", default="hg19", help="Geome build. Default: [%default]"),
	make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases"),
	make_option(c("--minMapScore"), type = "numeric", default=0.0, help="Include bins with a minimum mappability score of this value. Default: [%default]."),
	make_option(c("--maleChrXLogRThres"), type="numeric", default=-0.80, help = "ChrX Log ratio threshold to confirm as male gender."),
	make_option(c("--fracReadsInChrYForMale"), type="numeric", default=0.001, help = "Threshold for fraction of reads in chrY to assign as male. Default: [%default]"),
	make_option(c("-e", "--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions."),
	make_option(c("--method"), type = "character", default="median", help="Median or Mean."),
	make_option(c("--libdir"), type = "character", default=NULL, help = "Script library path. Usually exclude this argument unless custom modifications have been made to the ichorCNA R package code and the user would like to source those R files. Default: [%default]"),
	make_option(c("--ylim"), type = "character", default="c(-2,2)", help="Y-limits for plotting of mean/median log ratios"),
	make_option(c("--plotChrPanels"), type = "logical", default = FALSE, help = "Plot PoN values.")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

#id <- opt$id
gcWig <- opt$gcWig
mapWig <- opt$mapWig
repTimeWig <- opt$repTimeWig
filelist <- opt$filelist
exons.bed <- opt$exons.bed  # "0" if none specified
centromere <- opt$centromere
flankLength <- opt$rmCentromereFlankLength
method <- opt$method
outfile <- opt$outfile
genomeStyle <- opt$genomeStyle
genomeBuild <- opt$genomeBuild
minMapScore <- opt$minMapScore
libdir <- opt$libdir
ylim <- eval(parse(text = opt$ylim))
plotChrPanels <- opt$plotChrPanels
maleChrXLogRThres <- opt$maleChrXLogRThres
fracReadsInChrYForMale <- opt$fracReadsInChrYForMale
chrs <- as.character(eval(parse(text = opt$chrs)))
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize))); 
seqlevelsStyle(chrs) <- genomeStyle
seqlevelsStyle(chrNormalize) <- genomeStyle

## load ichorCNA library or source R scripts
if (!is.null(libdir) && libdir != "None"){
	source(paste0(libdir,"/R/utils.R"))
	source(paste0(libdir,"/R/segmentation.R"))
	source(paste0(libdir,"/R/EM.R"))
	source(paste0(libdir,"/R/output.R"))
	source(paste0(libdir,"/R/plotting.R"))
} else {
    library(ichorCNA)
}


if (!is.null(centromere)){
	centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
}
seqlevelsStyle(centromere$Chr) <- genomeStyle

files <- read.delim(filelist, header = FALSE, stringsAsFactors=FALSE, sep ="\t")[, 1]



## LOAD GC/MAP/REPTIME WIG FILES ###
message("Reading GC and mappability files")
gc <- wigToGRanges(gcWig)
if (is.null(gc)){
    stop("GC wig file not provided but is required")
}
map <- wigToGRanges(mapWig)
if (is.null(map)){
  message("No mappability wig file input, excluding from correction")
}
repTime <- wigToGRanges(repTimeWig)
if (is.null(repTime)){
  message("No replication timing wig file input, excluding from correction")
}else{
  if (mean(repTime$value, na.rm = TRUE) > 1){
    repTime$value <- repTime$value / 100 ## values in [0,1] - for LNCaP_repTime_10kb_hg38.txt
  }
}

normalGR <- NULL
info <- data.frame(Sample = character(), sex = character(), chrYcov = numeric(),
	chrXMedian = numeric(), normChrX = numeric())
for (i in 1:length(files)){
	### LOAD NORMAL FILES ###
	sid <- gsub(".wig","",basename(files[i]))
	message("Loading normal file:", files[i])
	normal_reads <- wigToGRanges(files[i])
		
	## FILTER BY EXONS IF PROVIDED ##
	## add gc and map to RangedData object ##
	if (!is.null(exons.bed)){
		targetedSequences <- read.delim(exons.bed, header=F, sep="\t", skip=86)
	}else{
		targetedSequences <- NULL
	}
	# normal_counts <- loadReadCountsFromWig(normal_reads, chrs=chrs, gc=gc, map=map, 
	# 				centromere=centromere, targetedSequences=targetedSequences)
	normal_counts <- loadReadCountsFromWig(normal_reads, chrs = chrs, gc = gc, map = map, 
									   repTime = repTime, centromere = centromere, 
									   flankLength = flankLength, targetedSequences = targetedSequences, 
									   chrXMedianForMale = maleChrXLogRThres, 
									   genomeStyle = genomeStyle, fracReadsInChrYForMale = fracReadsInChrYForMale,
                                       chrNormalize = chrNormalize, mapScoreThres = minMapScore)
	gender <- normal_counts$gender
	### CORRECT TUMOUR DATA FOR GC CONTENT AND MAPPABILITY BIASES ###
	message("Correcting ", sid, " sex: ", gender$gender, 
		" chrYcov: ", gender$chrYCovRatio, 
		" chrXMedian: ", gender$chrXMedian)
	if (is.null(normalGR)){
		normalGR <- normal_counts$counts
		values(normalGR) <- values(normalGR)$copy
		colnames(values(normalGR)) <- sid
	}else{
		values(normalGR)[[sid]] <- normal_counts$counts$copy
	}
	
	## Normalize chrX -- if male sex ##
	chrXStr <- grep("X", chrs, value = TRUE)
	chrXInd <- as.character(seqnames(normalGR)) == chrXStr		
	if (gender$gender == "male"){
		chrXMedian <- gender$chrXMedian		
		values(normalGR)[[sid]][chrXInd] <- values(normalGR)[[sid]][chrXInd] - chrXMedian
	}
	
	info <- rbind(info, data.frame(Sample = sid, sex = gender$gender, 
		chrYcov = gender$chrYCovRatio, chrXMedian = gender$chrXMedian,
		normChrX = median(values(normalGR)[[sid]][chrXInd], na.rm=T)))


}
print (info)

mat <- values(normalGR)
if (method == "median"){
   medianVal <- apply(mat, 1, median, na.rm = TRUE)
}else if (method == "mean"){
   medianVal <- apply(mat, 1, mean, na.rm = TRUE)
}else{
  stop("method is not specified as median or mean.")
}
values(normalGR)[["Median"]] <- medianVal

####################################
## Save the output ##
####################################

## save GR object ##
saveRDS(normalGR, file = paste0(outfile, "_", method, ".rds"))

## output to text file ##
write.table(as.data.frame(normalGR[,"Median"]), file=paste0(outfile, "_", method, ".txt"), col.names=TRUE, row.names=F, quote=F, sep="\t")

####################################
## Plots ##
####################################
## plot the median profile of the PoN
normalGR <- filterEmptyChr(normalGR)
chrLvls <- levels(seqnames(normalGR))
chrsToUse <- as.vector(seqnames(normalGR))
starts <- start(ranges(normalGR))
ends <- end(ranges(normalGR))
midPts <- starts + ((ends - starts) / 2) + 1
# get chr info from USCS using GenomeInfoDB #
seqInfo <- Seqinfo(genome = genomeBuild)
seqlevelsStyle(seqInfo) <- genomeStyle

#print(seqInfo)

seqInfo <- keepSeqlevels(seqInfo, chrs)
if (plotChrPanels){
	coord <- NULL
  	coord$posns <- midPts
}else{
	coord <- getGenomeWidePositions(chrsToUse, midPts, seqInfo)
}
mat <- data.frame(Chrs = chrsToUse, Positions = coord$posns, Median = normalGR$Median)
mat$Chrs <- factor(mat$Chrs, levels = chrLvls)

gp <- ggplot(mat, aes(x = Positions, y = Median)) +
	geom_point(alpha = 0.75, size = 1, shape = 19) +# + geom_line() 
	ylim(ylim) +
	theme(panel.background = element_blank(),
				panel.border = element_rect(colour = "black", fill = NA),
				panel.grid.major = element_blank(),#element_line(colour = "grey75", linetype="dashed"), 
				panel.grid.minor = element_blank(),#element_blank(), axis.title=element_text(size = 14, face = "bold"),
				axis.text = element_text(size = 12))


if (plotChrPanels){ # additional chr panel attributes #
	gp <- gp + facet_wrap( ~ Chrs, ncol = 4, scales = "free_x") + 
	  scale_x_continuous(breaks = NULL) + expand_limits(x = 0)
}else{
	numLines <- length(coord$chrBkpt)
	mid <- (coord$chrBkpt[1:(numLines - 1)] + coord$chrBkpt[2:numLines]) / 2
	gp <- gp + scale_x_continuous(name = "Chromosome", labels = chrs, breaks = mid, limits = c(coord$chrBkpt[1], tail(coord$chrBkpt ,1))) +
	  geom_vline(xintercept = coord$chrBkpt, linetype = "dotted")
}
outplotfile <- paste0(outfile, "_", method, ".png")
ggsave(gp, file = outplotfile, width = 16, height = 4)	





