# *ichorCNA*
ichorCNA is a tool for estimating the fraction of tumor in cell-free DNA from ultra-low-pass whole genome sequencing (ULP-WGS, 0.1x coverage). 

## Contacts
**Gavin Ha, Ph.D.** <gavinha@broadinstitute.org>  
**Justin Rhoades** <rhoades@broadinstitute.org>  
**Samuel Freeman** <sfreeman@broadinstitute.org>  
*Blood Biopsy Group, Cancer Program*  
Broad Institute of MIT and Harvard  
Dana-Farber Cancer Institute  
Koch Institute for integrative cancer research at MIT  

## Table of Contents
* [Description](#description)
* [Installation](#installation)
* [Usage](#usage)
* [Acknowledgements](#acknowledgements)

## Description
ichorCNA uses a probabilistic model, implemented as a hidden Markov model (HMM), to simultaneously segment the genome, predict large-scale copy number alterations, and estimate the tumor fraction of a ultra-low-pass whole genome sequencing sample (ULP-WGS).  
The analysis workflow consists of 3 major tasks:  
1. Computing read coverage  
  - 
2. Data normalization  
3. CNA prediction and estimation of tumor fraction  

## Installation
1. Checkout the latest release of ichorCNA from GitHub  
  ```
  
  git clone git@github.com:broadinstitute/ichorCNA.git  
  
  ```  
2. Install R dependencies (in R)
  ```
  ## install from CRAN
  install.packages("plyr") # version > 1.8.4
  ## install packages from
  source("https://bioconductor.org/biocLite.R")
  biocLite("HMMcopy") # version >= 1.14.0
  biocLite("GenomeInfoDb") # version >= 1.8.7
  biocLite("quantsmooth") # version >= 1.38.0
  ```

3. Install the ichorCNA R package.   
  ```
  ## from the command line where ichorCNA github was cloned.
  R CMD INSTALL ichorCNA  
  
  ```  
4. Install the HMMcopy suite from <http://compbio.bccrc.ca/software/hmmcopy/>  
Please follow instructions on the HMMcopy website.


## Usage
Here is a demo script to run ichorCNA on some example data  
  ```
  sh run_ichorCNA_analysis.sh
  ```

## Acknowledgements
ichorCNA is developed and maintained by Gavin Ha (<gavinha@broadinstitute.org>), Justin Rhoades (<rhoades@broadinstitute.org>), and Sam Freeman (<sfreeman@broadinstitute.org>).  

This work was done in collaboration between  
- **Blood Biopsy Group** at the Broad Institute of MIT and Harvard
- Laboratory of **Matthew Meyerson**, Medical Oncology, Dana-Farber Cancer Institute
- Laboratory of **J. Christopher Love**, Koch Institute for integrative cancer research at MIT
- Laboratory of **Gad Getz**, Cancer Program, Broad Institute
