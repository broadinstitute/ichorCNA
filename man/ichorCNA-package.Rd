\name{ichorCNA-package}
\alias{ichorCNA-package}
\alias{ichorCNA}
\docType{package}
\title{
ichorCNA
}
\description{
ichorCNA is a tool for estimating the fraction of tumor in ultra-low-pass whole genome sequencing (ULP-WGS) of cell-free DNA.
}
\details{

\tabular{ll}{
Package: \tab ichorCNA\cr
Type: \tab Package\cr
Version: \tab 1.12.1\cr
Date: \tab 2017-07-03\cr
License: \tab GPL-3\cr
}

\author{
Gavin Ha <gavinha@broadinstitute.org>
Justin Rhoades <rhoades@broadinstitute.org>
Samuel Freeman <sfreeman@broadinstitute.org>
}

\examples{
library(HMMcopy)
tumor_file <- system.file("extdata", "MBC_315.ctDNA.reads.wig",
                      package = "ichorCNA")

tumor_reads <- wigToRangedData(tumor_file)
}

