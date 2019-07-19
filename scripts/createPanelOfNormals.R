library(HMMcopy)
library(GenomicRanges)
library(optparse)
library(ichorCNA)

options(stringsAsFactors=FALSE, scipen=0)
options(bitmapType='cairo')

option_list <- list(
	make_option(c("--gcWig"), type = "character", help = "GC Wig file for reference genome"),
	make_option(c("--mapWig"), type = "character", default=NULL, help = "Mappabiliy Wig file for reference genome"),
	make_option(c("-f", "--filelist"), type = "character", help = "List of of wig files."),
	make_option(c("-o", "--outfile"), type = "character", help = "Output file."),
	make_option(c("-c", "--centromere"), type="character", help = "File containing Centromere locations"),
	make_option(c("--chrs"), type="character", default="c(1:22,\"X\")", help = "Specify chromosomes to analyze."),
  make_option(c("--genomeStyle"), type = "character", default = "NCBI", help = "NCBI or UCSC chromosome naming convention; use UCSC if desired output is to have \"chr\" string. [Default: %default]"),
	make_option(c("--chrNormalize"), type="character", default="c(1:22)", help = "Specify chromosomes to normalize GC/mappability biases"),
	make_option(c("--maleChrXLogRThres"), type="numeric", default=-0.80, help = "ChrX Log ratio threshold to confirm as male gender."),
	make_option(c("-e", "--exons.bed"), type = "character", default=NULL, help = "Path to bed file containing exon regions."),
	make_option(c("--method"), type = "character", default="median", help="Median or Mean.")
	#make_option(c("--ylim"), type = "character", default="c(-2,2)", help="Y-limits for plotting of mean/median log ratios")
)
parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

#id <- opt$id
gcWig <- opt$gcWig
mapWig <- opt$mapWig
filelist <- opt$filelist
exons.bed <- opt$exons.bed  # "0" if none specified
centromere <- opt$centromere
method <- opt$method
outfile <- opt$outfile
genomeStyle <- opt$genomeStyle
#ylim <- eval(parse(text = opt$ylim))
maleChrXLogRThres <- opt$maleChrXLogRThres
chrs <- as.character(eval(parse(text = opt$chrs)))
chrNormalize <- as.character(eval(parse(text=opt$chrNormalize))); 
seqlevelsStyle(chrs) <- genomeStyle
seqlevelsStyle(chrNormalize) <- genomeStyle

if (!is.null(centromere)){
	centromere <- read.delim(centromere,header=T,stringsAsFactors=F,sep="\t")
}
seqlevelsStyle(centromere$Chr) <- genomeStyle

files <- read.delim(filelist, header = FALSE, stringsAsFactors=FALSE, sep ="\t")[, 1]

normalGR <- NULL

for (i in 1:length(files)){
	### LOAD NORMAL FILES ###
	sid <- gsub(".wig","",basename(files[i]))
	message("Loading normal file:", files[i])
	normal_reads <- wigToGRanges(files[i])

	message("Reading GC and mappability files")
	if (is.null(gcWig)) {
	  stop("GC Wig file required")
	} else {
	  gc <- wigToGRanges(gcWig)
	}
	if (is.null(mapWig)) {
	  message("Normalizing without mappability Wig file.")
	  map <- NULL
	} else {
	  map <- wigToGRanges(mapWig)
	}
		
	## FILTER BY EXONS IF PROVIDED ##
	## add gc and map to RangedData object ##
	if (!is.null(exons.bed)){
		targetedSequences <- read.delim(exons.bed, header=F, sep="\t", skip=86)
	}else{
		targetedSequences <- NULL
	}
	normal_counts <- loadReadCountsFromWig(normal_reads, chrs=chrs, gc=gc, map=map, 
					centromere=centromere, targetedSequences=targetedSequences)
	
	gender <- normal_counts$gender

	### CORRECT TUMOUR DATA FOR GC CONTENT AND MAPPABILITY BIASES ###
	message("Correcting ", sid)
	if (is.null(normalGR)){
		normalGR <- normal_counts$counts
		values(normalGR) <- values(normalGR)$copy
		colnames(values(normalGR)) <- sid
	}else{
		values(normalGR)[[sid]] <- normal_counts$counts$copy
	}
	
	chrXMedian <- gender$chrXMedian
	chrXStr <- grep("X", chrs, value = TRUE)
	chrXInd <- as.character(seqnames(normalGR)) == chrXStr
	## Normalize chrX ##
	values(normalGR)[[sid]][chrXInd] <- values(normalGR)[[sid]][chrXInd] - chrXMedian
	
}


mat <- values(normalGR)
if (method == "median"){
   medianVal <- apply(mat, 1, median, na.rm = TRUE)
}else if (method == "mean"){
   medianVal <- apply(mat, 1, mean, na.rm = TRUE)
}else{
  stop("method is not specified as median or mean.")
}
values(normalGR)[["Median"]] <- medianVal

## output to text file ##
write.table(as.data.frame(normalGR[,"Median"]), file=paste0(outfile, "_", method, ".txt"), col.names=TRUE, row.names=F, quote=F, sep="\t")

## save GR object ##
saveRDS(normalGR, file = paste0(outfile, "_", method, ".rds"))
