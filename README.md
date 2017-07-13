# *ichorCNA*
ichorCNA is a tool for estimating the fraction of tumor in cell-free DNA from ultra-low-pass whole genome sequencing (ULP-WGS, 0.1x coverage). 

## Table of Contents
* [Description](#description)
* [Contacts](#contacts)
* [ichorCNA Wiki Page](#ichorcna-wiki-page)
* [Acknowledgements](#acknowledgements)

## Description
ichorCNA uses a probabilistic model, implemented as a hidden Markov model (HMM), to simultaneously segment the genome, predict large-scale copy number alterations, and estimate the tumor fraction of a ultra-low-pass whole genome sequencing sample (ULP-WGS). 

The methodology and probabilistic model are described in:  
Adalsteinsson, Ha, Freeman, et al. Scalable whole-exome sequencing of cell-free DNA reveals high concordance with metastatic tumors. (2017) Nature Communications, in press.

The analysis workflow consists of 3 major tasks:  
1. Computing read coverage from ULP-WGS  
2. Data normalization  
3. CNA prediction and estimation of tumor fraction of cfDNA

## ichorCNA Wiki Page
**For more details on usage/pipelines, outputs, and FAQs, please visit the [GitHub Wiki page for ichorCNA](https://github.com/broadinstitute/ichorCNA/wiki)**

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
