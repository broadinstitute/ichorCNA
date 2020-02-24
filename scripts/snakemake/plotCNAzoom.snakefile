configfile: "config/configPlotZoom.yaml"
configfile: "config/samples.yaml"

rule all:
  input:
  	expand("results/plotCNAzoom/{plotID}/{tumor}_CNA-SV_{chr}-{start}-{end}.{format}", tumor=config["samples"], plotID=config["plot_id"], chr=config["plot_chr"], start=config["plot_startPos"], end=config["plot_endPos"], format=config["plot_format"])


rule plotCNAzoom:
	input:
		cna="results/ichorCNA/{tumor}/{tumor}.cna.seg",
		params="results/ichorCNA/{tumor}/{tumor}.params.txt",
	output:
		"results/plotCNAzoom/{plotID}/{tumor}_CNA-SV_{chr}-{start}-{end}.{format}"
	params:
		plotCNscript=config["plotCN_script"],
		#plotfuncs=config["plot_funcs"],
		genomeBuild=config["genomeBuild"],
		genomeStyle=config["genomeStyle"],
		zoom=config["plot_zoom"],
		yaxis=config["plot_yaxis"],
		chr=config["plot_chr"],
		start=config["plot_startPos"],
		end=config["plot_endPos"],
		ylim=config["plot_ylim"],
		geneFile=config["plot_geneFile"],
		size=config["plot_size"],
		format=config["plot_format"]
	log:
		"logs/plotCNAzoom/{plotID}/{tumor}_{chr}-{start}-{end}.{format}.log"
	shell:
		"Rscript {params.plotCNscript} --id {wildcards.tumor} --cnFile {input.cna} --paramFile {input.params} --chrs \"{params.chr}\" --genomeBuild {params.genomeBuild} --genomeStyle {params.genomeStyle} --start {params.start} --end {params.end} --zoom {params.zoom} --yaxis {params.yaxis} --plotYlim \"{params.ylim}\" --geneFile {params.geneFile} --plotSize \"{params.size}\" --outPlotFile {output} > {log} 2> {log}"
