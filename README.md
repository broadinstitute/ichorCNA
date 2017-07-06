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

The methodology and probabilistic model are described in:  
Adalsteinsson, Ha, Freeman, et al. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. (2017) Nature Communications, in press.

The analysis workflow consists of 3 major tasks:  
1. Computing read coverage from ULP-WGS  
2. Data normalization  
3. CNA prediction and estimation of tumor fraction of cfDNA

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

Both steps have been compiled into a Python Snakemake pipeline that accepts BAM files and performs all steps to generate ichorCNA results. **Please see [ichorCNA snakemake pipeline](scripts/snakemake).**


## Contacts
**Gavin Ha, Ph.D.** <gavinha@broadinstitute.org>  
**Justin Rhoades** <rhoades@broadinstitute.org>  
**Samuel Freeman** <sfreeman@broadinstitute.org>  
**Blood Biopsy Group, Cancer Program** <bloodbiopsy@broadinstitute.org>  
Broad Institute of MIT and Harvard  
Koch Institute for integrative cancer research at MIT  
Dana-Farber Cancer Institute  

## Acknowledgements
ichorCNA is developed and maintained by Gavin Ha, Justin Rhoades, and Sam Freeman.  

This work was done in collaboration with  
- **Blood Biopsy Group**, Group Leader Viktor Adalsteinsson, Broad Institute of MIT and Harvard
- Laboratory of **Matthew Meyerson**, Medical Oncology, Dana-Farber Cancer Institute
- Laboratory of **J. Christopher Love**, Koch Institute for integrative cancer research at MIT
- Laboratory of **Gad Getz**, Cancer Program, Broad Institute

## Software License
