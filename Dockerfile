FROM ubuntu:14.04
MAINTAINER Jonathan Dursi <jonathan@dursi.ca>
LABEL Description="Consensus calling pipeline for PCAWG"

VOLUME /data
WORKDIR /data

# get ubuntu packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libxml2-dev \
        libz-dev \
        python \
        python-dev \
        python-pip \
        tabix \
        wget \
        zlib1g-dev 

# samtools - for indexing reference, etc
RUN mkdir -p /deps && \
    cd /deps && \
    wget -nv https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 && \
    tar -xjvf samtools-1.3.tar.bz2 && \
    rm samtools-1.3.tar.bz2 && \
    cd samtools-1.3 && \
    make prefix=/usr/local/ install && \
    cd .. && \
    rm -rf samtools-1.3

# get pyvcf for annotate_from_readcounts.py
RUN pip install --upgrade pip && \
    pip install pyvcf

# install vcfanno for annotating VCFs with beds/vcfs
RUN cd /tmp \
    && wget -nv https://github.com/brentp/vcfanno/releases/download/v0.1.0/vcfanno_0.1.0_linux_amd64.tar.gz \
    && tar -xzvf vcfanno_0.1.0_linux_amd64.tar.gz \
    && mv vcfanno_0.1.0_linux_amd64/vcfanno /usr/local/bin \
    && rm -rf vcfanno_0.1.0_linux_amd64.tar.gz vcfanno_0.1.0_linux_amd64/vcfanno

# install mergevcf for merging the callsets
RUN cd /tmp \
    && wget -nv https://github.com/ljdursi/mergevcf/archive/0.2.tar.gz \
    && tar -xzf 0.2.tar.gz \
    && rm 0.2.tar.gz \
    && cd mergevcf-0.2 \
    && python setup.py install \
    && cd .. \
    && rm -rf mergevcf-0.2

COPY clean_snv_calls.py /usr/local/bin
COPY clean_indel_calls.py /usr/local/bin
COPY dbsnp_annotate_one.sh /usr/local/bin
COPY merge-one-tumour-snv.sh /usr/local/bin
COPY consensus_snv.sh /usr/local/bin
COPY consensus_indel.sh /usr/local/bin
COPY wrapper.sh /usr/local/bin
COPY wrapper.sh /usr/local/bin
COPY build_dbs.sh /usr/local/bin

# install vt (for normalizing VCFs )
RUN cd /tmp \
    && wget -nv https://github.com/atks/vt/archive/0.5772.tar.gz \
    && tar -xzf 0.5772.tar.gz \
    && rm 0.5772.tar.gz \
    && cd vt-0.5772 \
    && make \
    && mv vt /usr/local/bin \
    && cd .. \
    && rm -rf vt-0.5772

###
### R and various packages are needed for the model filtering
###

RUN add-apt-repository http://cran.utstat.utoronto.ca/bin/linux/ubuntu/ \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base r-base-dev 

COPY Rdeps.R /deps

RUN Rscript /deps/Rdeps.R

COPY models /dbs
COPY analysis /usr/local/bin/
COPY filter /usr/local/bin/

ENTRYPOINT ["/usr/local/bin/wrapper.sh"]
