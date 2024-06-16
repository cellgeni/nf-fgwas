FROM ubuntu:22.04

LABEL org.opencontainers.image.authors="Patrick Pett <jan.patick.pett@sanger.ac.uk>; Cellular Genetics Informatics <cellgeni@sanger.ac.uk>"
LABEL org.opencontainers.image.title="Functional GWAS for single cell"
LABEL org.opencontainers.image.description="Tools for integrating GWAS results with scRNA-seq data to identify disease associated cells"
LABEL org.opencontainers.image.source="https://github.com/cellgeni/nf-fgwas"
LABEL org.opencontainers.image.licenses="MIT"


ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
        build-essential perl curl git software-properties-common dirmngr pkg-config \
        python3 python3-venv python3-pip python-is-python3 \
        libblas-dev liblapack-dev gfortran \
        libbz2-dev libcurl4-openssl-dev libgsl0-dev \
        liblzma-dev libncurses5-dev libperl-dev libssl-dev zlib1g-dev

# vcftools
ENV VCFTOOLS_VER="0.1.16"
RUN curl -O -L https://github.com/vcftools/vcftools/releases/download/v${VCFTOOLS_VER}/vcftools-${VCFTOOLS_VER}.tar.gz && \
    tar -zxf vcftools-${VCFTOOLS_VER}.tar.gz -C /tmp && \
    rm -rf vcftools-${VCFTOOLS_VER}.tar.gz && \
    cd /tmp/vcftools-${VCFTOOLS_VER} && \
    ./configure --prefix=/opt/vcftools-${VCFTOOLS_VER} && \
    make && make install && \
    rm -rf /tmp/vcftools-${VCFTOOLS_VER} 
ENV PATH=${PATH}:/opt/vcftools-${VCFTOOLS_VER}/bin
ENV PERL5LIB=/opt/vcftools-${VCFTOOLS_VER}/share/perl:${PERL5LIB}

# htslib
ENV HTSLIB_VER="1.20"
RUN curl -O -L https:/github.com/samtools/htslib/releases/download/${HTSLIB_VER}/htslib-${HTSLIB_VER}.tar.bz2 && \
   tar -xjf htslib-${HTSLIB_VER}.tar.bz2 -C /opt && \
   rm -rf htslib-${HTSLIB_VER}.tar.bz2 && \
   cd /opt/htslib-${HTSLIB_VER} && \
   ./configure --enable-plugins --with-plugin-path='$(libexecdir)/htslib:/usr/libexec/htslib' && \
   make install && ldconfig

# samtools
RUN curl -O -L https://github.com/samtools/samtools/releases/download/${HTSLIB_VER}/samtools-${HTSLIB_VER}.tar.bz2 && \
    tar -xjf samtools-${HTSLIB_VER}.tar.bz2 -C /opt && \
    rm -rf samtools-${HTSLIB_VER}.tar.bz2 && \
    cd /opt/samtools-${HTSLIB_VER} && \
    ./configure --with-htslib=system && \
    make install

# bcftools
RUN curl -O -L https://github.com/samtools/bcftools/releases/download/${HTSLIB_VER}/bcftools-${HTSLIB_VER}.tar.bz2 && \
    tar -xjf bcftools-${HTSLIB_VER}.tar.bz2 -C /opt && \
    rm -rf bcftools-${HTSLIB_VER}.tar.bz2 && \
    cd /opt/bcftools-${HTSLIB_VER} && \
    ./configure --enable-libgsl --enable-perl-filters --with-htslib=system && \
    make install && \
    git clone --depth 1 https://github.com/samtools/htslib-plugins.git /opt/htslib-plugins && \
    cd /opt/htslib-plugins && \
    make PLUGINS='hfile_cip.so hfile_mmap.so' install

# R and R packages
RUN curl -L https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get update && apt-get install -y r-base r-recommended
RUN Rscript -e 'install.packages(c("qvalue","tidyr","plyr","readr","dplyr","tibble","futile.logger","AnnotationHub"))'

# python packages
RUN pip install --no-cache pandas pyarrow

# PHM
RUN git clone https://github.com/natsuhiko/PHM.git /opt/PHM
# RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
#     bash Miniforge3-Linux-x86_64.sh -b -p "${HOME}/conda" && \
#     source "${HOME}/conda/etc/profile.d/conda.sh" && \
#     source "${HOME}/conda/etc/profile.d/mamba.sh" && \
#     mamba activate
# RUN mamba install -y conda-forge::f2c conda-forge::clapack
# RUN apt-get install -y f2c
RUN curl -O -L https://www.netlib.org/clapack/clapack.tgz && \
    tar zxvf clapack.tgz && \
    cd CLAPACK-3.2.1 && \
    mv make.inc.example make.inc && \
    export CFLAGS="$CFLAGS -fcommon" && \
    make && \
    ln -s lapack_LINUX.a liblapack.a && \
    ln -s tmglib_LINUX.a libtmglib.a && \
    ln -s blas_LINUX.a libblas.a
RUN cd /opt/PHM/src && make && make install && ln -s /opt/PHM/bin/hm /opt/PHM/bin/fgwas_hm
ENV PATH=/opt/PHM/bin:$PATH

# getRsq
ADD getRsq.tar.gz /opt
RUN cd /opt/getRsq/src && make
ENV PATH=/opt/getRsq/src:$PATH
