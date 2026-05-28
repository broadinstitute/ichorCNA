# ichorCNA docker image — for the IchorCNA WDL workflow.
#
# Provides:
#   - R 4.5 with ichorCNA + HMMcopy + GenomicRanges + GenomeInfoDb + plyr + optparse
#   - samtools 1.20 with S3 (HTSlib) support — required to stream CRAMs from S3
#   - Python 3 + scripts/bam_to_wig.py — stdin SAM → WIG, no BAM index needed
#   - Bundled ichorCNA reference files (GC/map WIGs, hg38 PoN) under /opt/ichorCNA/inst/extdata

FROM debian:bookworm-slim

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        wget \
        git \
        gnupg \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libbz2-dev \
        liblzma-dev \
        libdeflate-dev \
        zlib1g-dev \
        libncurses-dev \
        libreadline-dev \
        libpcre2-dev \
        libffi-dev \
        python3 \
        python3-pip \
        gfortran \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# CRAN apt repo for R 4.5
RUN gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' && \
    gpg --armor --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7' > /etc/apt/trusted.gpg.d/cran.asc && \
    echo "deb https://cloud.r-project.org/bin/linux/debian bookworm-cran40/" > /etc/apt/sources.list.d/cran.list && \
    apt-get update && apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# samtools built with HTSlib S3/CRAM support
ARG HTSLIB_VERSION=1.20
ARG SAMTOOLS_VERSION=1.20
RUN wget -q https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --enable-libcurl --enable-s3 --enable-gcs && \
    make -j"$(nproc)" && make install && ldconfig && \
    cd .. && rm -rf htslib-${HTSLIB_VERSION} htslib-${HTSLIB_VERSION}.tar.bz2 && \
    wget -q https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --with-htslib=system && \
    make -j"$(nproc)" && make install && \
    cd .. && rm -rf samtools-${SAMTOOLS_VERSION} samtools-${SAMTOOLS_VERSION}.tar.bz2

# R packages: install pre-compiled CRAN packages where possible, then ichorCNA deps
RUN apt-get update && apt-get install -y --no-install-recommends \
        r-cran-optparse \
        r-cran-plyr \
        r-cran-data.table \
        r-cran-bitops \
        r-cran-rcurl \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# BiocManager + Bioconductor deps required by ichorCNA / HMMcopy
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/', Ncpus=$(nproc))" && \
    R -e "BiocManager::install(c('HMMcopy','GenomicRanges','GenomeInfoDb','IRanges','S4Vectors','BiocGenerics'), ask=FALSE, update=FALSE, Ncpus=$(nproc))"

# Install ichorCNA itself (R package + scripts + reference data) into /opt
COPY . /opt/ichorCNA
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org/'); remotes::install_local('/opt/ichorCNA', upgrade='never')"

# bam_to_wig.py is invoked from /opt/ichorCNA/scripts/bam_to_wig.py — already
# in place via the COPY above. Also expose Rscript and python3 on PATH.
ENV PATH=/usr/local/bin:/usr/bin:/bin

# Sanity check at build time
RUN samtools --version | head -1 && \
    Rscript -e 'stopifnot(requireNamespace("ichorCNA"), requireNamespace("HMMcopy"))' && \
    python3 -c "import sys; assert sys.version_info >= (3,9)" && \
    test -x /opt/ichorCNA/scripts/bam_to_wig.py || chmod +x /opt/ichorCNA/scripts/bam_to_wig.py

WORKDIR /work
