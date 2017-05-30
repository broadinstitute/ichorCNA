configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule correctDepth:
  input: 
  	expand("results/ichorCNA/{tumor}/{tumor}.cna.seg", tumor=config["pairings"]),
  	expand("results/readDepth/{samples}.bin{binSize}.wig", samples=config["samples"], binSize=str(config["binSize"]))

rule read_counter:
	input:
		lambda wildcards: config["samples"][wildcards.samples]
	output:
		"results/readDepth/{samples}.bin{binSize}.wig"		
	params:
		readCounter=config["readCounterScript"],
		binSize=config["binSize"],
		qual="20",
		chrs=config["chrs"]
	resources:
		mem=4
	log:
		"logs/readDepth/{samples}.log"
	shell:
		"{params.readCounter} {input} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule ichorCNA:
	input:
		tum="results/readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
		norm=lambda wildcards: "results/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
	output:
		#corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
		#param="results/ichorCNA/{tumor}/{tumor}.params.txt",
		#cna="results/ichorCNA/{tumor}/{tumor}.cna.seg",
		#segTxt="results/ichorCNA/{tumor}/{tumor}.seg.txt",
		#seg="results/ichorCNA/{tumor}/{tumor}.seg",
		#rdata="results/ichorCNA/{tumor}/{tumor}.RData",
		outDir="results/ichorCNA/{tumor}/",
	params:
		rscript=config["ichorCNA_rscript"],
		libdir=config["ichorCNA_libdir"],
		datadir=config["ichorCNA_datadir"],
		id="{tumor}",
		ploidy=config["ichorCNA_ploidy"],
		normal=config["ichorCNA_normal"],
		lda="1000",
		maxCN=config["ichorCNA_maxCN"],
		chrs=config["ichorCNA_chrs"],
		chrTrain=config["ichorCNA_chrTrain"],
		centromere=config["ichorCNA_centromere"],
		exons=config["ichorCNA_exons"],
		txnE=config["ichorCNA_txnE"],
		txnStrength=config["ichorCNA_txnStrength"],
		fracReadsChrYMale="0.001",
		plotFileType=config["ichorCNA_plotFileType"],
		plotYlim=config["ichorCNA_plotYlim"]
	resources:
		mem=4
	log:
		"logs/ichorCNA/{tumor}.log"	
	shell:
		"Rscript {params.rscript} --libdir {params.libdir} --datadir {params.datadir} --id {params.id} --WIG {input.tum} --NORMWIG {input.norm} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --lambda {params.lda} --maxCN {params.maxCN} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --fracReadsInChrYForMale {params.fracReadsChrYMale} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {output.outDir} > {log} 2> {log}"

