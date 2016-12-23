FROM ubuntu:14.04
MAINTAINER Jonathan Dursi <jonathan@dursi.ca>
LABEL Description="Consensus calling pipeline for PCAWG"

VOLUME /data
WORKDIR /data

# get ubuntu packages
RUN apt-get update && \
    apt-get install -y \
        automake \
        autotools-dev \
        build-essential \
        cmake \
        git \
        libcurl4-openssl-dev \
        libhts-dev \
        libhts0 \
        libjemalloc-dev \
        libncurses5-dev \
        libsparsehash-dev \
        libxml2-dev \
        libz-dev \
        python \
        python-dev \
        python-pip \
        tabix \
        wget \
        zlib1g-dev 

# build remaining dependencies:
# bamtools - for SGA
RUN mkdir -p /deps && \
    cd /deps && \
    wget -nv https://github.com/pezmaster31/bamtools/archive/v2.4.0.tar.gz && \
    tar -xzvf v2.4.0.tar.gz && \
    rm -rf v2.4.0.tar.gz && \
    cd bamtools-2.4.0 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# bam-readcount: for SNV annotation
RUN mkdir -p /deps && \
    cd /deps && \
    wget -nv https://github.com/genome/bam-readcount/archive/v0.7.4.tar.gz && \
    tar -xzvf v0.7.4.tar.gz && \
    rm v0.7.4.tar.gz && \
    cd bam-readcount-0.7.4 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

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

# vcflib - for tools like vcfbreakmulti
# set a fixed version for reproducibility
RUN mkdir -p /deps && \
    cd /deps && \
    git clone --recursive git://github.com/ekg/vcflib.git && \
    cd vcflib && \
    git checkout d453d91592fe8a74d92b49cd6c7cd73f79a8b70b && \
    make 

# get pancan standard reference
RUN mkdir -p /reference && \
    cd /reference && \
    wget -nv ftp://ftp.sanger.ac.uk/pub/project/PanCancer/genome.fa.gz && \
    gunzip genome.fa.gz ; \
    bgzip genome.fa && \
    samtools faidx genome.fa.gz

# build SGA
RUN mkdir -p /src && \
    cd /src && \
    wget -nv https://github.com/jts/sga/archive/v0.10.14.tar.gz && \
    tar -xzvf v0.10.14.tar.gz && \
    rm v0.10.14.tar.gz && \
    cd sga-0.10.14/src && \
    ./autogen.sh && \
    ./configure --with-bamtools=/deps/bamtools-2.4.0 --with-jemalloc=/usr --prefix=/usr/local && \
    make && \
    make install

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
COPY dbsnp_annotate_one.sh /usr/local/bin
COPY merge-one-tumour-snv.sh /usr/local/bin
COPY consensus_snv.sh /usr/local/bin

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

# copy annotation-database-copying scripts

ENV DBDIR /dbs/annotation_databases/data
RUN mkdir -p $DBDIR

ENV DBSNP $DBDIR/All_20160601.vcf
ENV RMSK $DBDIR/hg19.rmsk.bed

RUN wget -nv ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz -O $DBSNP.gz \
    && gunzip $DBSNP.gz \
    && bgzip $DBSNP \
    && tabix -p vcf $DBSNP.gz

RUN wget -nv http://people.virginia.edu/~arq5x/files/gemini/annotations/hg19.rmsk.bed.gz -O $RMSK.gz \
    && gunzip $RMSK.gz \
    && bgzip $RMSK \
    && tabix -p bed $RMSK.gz 

ENV VCF $DBDIR/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
ENV NEWVCF $DBDIR/1000genomes.phase3.decomposed.normalized.vcf.gz
ENV REF /reference/genome.fa.gz

RUN wget -nv ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -O $VCF \
    && /usr/local/bin/vt decompose -s "$VCF" \
        | /usr/local/bin/vt normalize -r "${REF}" - \
        | grep -v "##INFO=<ID=MEINFO" \
        | sed -e 's/MEINFO=[^;]*//' \
        | bgzip -f > $NEWVCF \
    && tabix -p vcf $NEWVCF \
    && rm -f $VCF

###
### R and various packages are needed for the model filtering
###

RUN echo "deb http://cran.utstat.utoronto.ca/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-base r-base-dev 

COPY Rdeps.R /deps

RUN Rscript /deps/Rdeps.R

COPY models /dbs

ENTRYPOINT ["/usr/local/bin/consensus_snv.sh"]
