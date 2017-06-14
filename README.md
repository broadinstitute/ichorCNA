# *ichorCNA*
ichorCNA is a tool for estimating the fraction of tumor in cell-free DNA from ultra-low-pass whole genome sequencing (ULP-WGS, 0.1x coverage). 

## Table of Contents
* [Description](#description)
* [Installation](#installation)
* [Usage](#usage)
* [Contacts](#contacts)
* [Acknowledgements](#acknowledgements)

## Description
ichorCNA uses a probabilistic model, implemented as a hidden Markov model (HMM), to simultaneously segment the genome, predict large-scale copy number alterations, and estimate the tumor fraction of a ultra-low-pass whole genome sequencing sample (ULP-WGS).  
The analysis workflow consists of 3 major tasks:  
1. Computing read coverage  
  - 
2. Data normalization  
3. CNA prediction and estimation of tumor fraction  

## Installation
1. Using R devtools
	```
	install.packages("devtools")
	library(devtools)
	install_github("broadinstitute/ichorCNA")
	```
2. Manual installation  
    a. Checkout the latest release of ichorCNA from GitHub  
    ```
    git clone git@github.com:broadinstitute/ichorCNA.git  
    ```  
    b. Install R dependencies (in R)  

    ```
    ## install from CRAN
    install.packages("plyr") # version > 1.8.4
    ## install packages from
    source("https://bioconductor.org/biocLite.R")
    biocLite("HMMcopy") # version >= 1.14.0
    biocLite("GenomeInfoDb") # version >= 1.8.7
    ```  

    c. Install the ichorCNA R package  

    ```
    ## from the command line and in the directory where ichorCNA github was cloned.
    R CMD INSTALL ichorCNA  
    ```  
        
3. Other dependencies  
  a. Install the HMMcopy suite from <http://compbio.bccrc.ca/software/hmmcopy/>  
    Please follow instructions on the HMMcopy website.


## Usage
There are 2 main steps in the analysis workflow:
1. Generating read count coverage information using `readCounter` from the HMMcopy suite.
2. Copy number analysis and prediction of tumor fraction using ichorCNA R package.

Both steps have been compiled into a Python Snakemake pipeline that accepts BAM files and performs all steps to generate ichorCNA results. Please see [ichorCNA snakemake pipeline](scripts/snakemake).

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

## Contacts
**Gavin Ha, Ph.D.** <gavinha@broadinstitute.org>  
**Justin Rhoades** <rhoades@broadinstitute.org>  
**Samuel Freeman** <sfreeman@broadinstitute.org>  
*Blood Biopsy Group, Cancer Program*  
Broad Institute of MIT and Harvard  
Dana-Farber Cancer Institute  
Koch Institute for integrative cancer research at MIT  

## Acknowledgements
ichorCNA is developed and maintained by Gavin Ha (<gavinha@broadinstitute.org>), Justin Rhoades (<rhoades@broadinstitute.org>), and Sam Freeman (<sfreeman@broadinstitute.org>).  

This work was done in collaboration between  
- **Blood Biopsy Group** at the Broad Institute of MIT and Harvard
- Laboratory of **Matthew Meyerson**, Medical Oncology, Dana-Farber Cancer Institute
- Laboratory of **J. Christopher Love**, Koch Institute for integrative cancer research at MIT
- Laboratory of **Gad Getz**, Cancer Program, Broad Institute
