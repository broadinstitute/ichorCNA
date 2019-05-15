## Contents

### 1. Snakemake Pipeline ###
[See instructions for snakemake pipeline](snakemake)

### 2. ichorCNA R script details ###

The main [ichorCNA R script](scripts/runIchorCNA.R) can also be run with even more exposed options once WIG files for read counts have been generated.

Here is an example of how to launch the R script from the command line:
```
Rscript ../runIchorCNA.R --libdir ../../R/ --datadir ../../inst/extdata/ \
  --id tumor_sample1 --WIG results/readDepth/tumor_sample1.bin1000000.wig \
  --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --centromere ../../inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --exons.bed None --txnE 0.9999 --txnStrength 10000 --plotFileType png --plotYLim "c(-2,4)" \
  --outDir results/ichorCNA/tumor_sample1/ > logs/ichorCNA/tumor_sample1.log 2> logs/ichorCNA/tumor_sample1.log
```

The list of arguments that can be passed to the script by invoking `--help`:
```
>Rscript runIchorCNA.R --help

Usage: runIchorCNA.R [options]


Options:
        -l LIBDIR, --libdir=LIBDIR
                Script library path

        --datadir=DATADIR
                Reference wig dir path

        -t WIG, --WIG=WIG
                Path to tumor WIG file.

        --NORMWIG=NORMWIG
                Path to normal WIG file.

        --normalPanel=NORMALPANEL
                Median corrected depth from panel of normals

        -e EXONS.BED, --exons.bed=EXONS.BED
                Path to bed file containing exon regions.

        --id=ID
                Patient ID.


        -n NORMAL, --normal=NORMAL
                Initial normal contamination

        --scStates=SCSTATES
                Subclonal states to consider

        -p PLOIDY, --ploidy=PLOIDY
                Initial tumour ploidy

        -m MAXCN, --maxCN=MAXCN
                Total clonal CN states

        --estimateNormal=ESTIMATENORMAL
                Estimate normal.

        --estimateScPrevalence=ESTIMATESCPREVALENCE
                Estimate subclonal prevalence.

        --estimatePloidy=ESTIMATEPLOIDY
                Estimate tumour ploidy.

 (plus more arguments)
```
