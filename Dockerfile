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

# get pancan standard reference
RUN mkdir -p /reference && \
    cd /reference && \
    wget -nv ftp://ftp.sanger.ac.uk/pub/project/PanCancer/genome.fa.gz && \
    gunzip genome.fa.gz ; \
    bgzip genome.fa && \
    samtools faidx genome.fa.gz

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

RUN wget -nv ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz -O - \
    | zcat | bgzip -f > $DBSNP.gz ;\
    tabix -p vcf $DBSNP.gz

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
COPY analysis /usr/local/bin/

COPY filter /usr/local/bin/

ENTRYPOINT ["/usr/local/bin/consensus_snv.sh"]
