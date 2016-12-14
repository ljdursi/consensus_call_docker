#!/bin/bash

# requires wget, vt, bgzip, tabix
readonly DBDIR=/dbs/annotation_databases/data
readonly VCF=${DBDIR}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
readonly REF=/reference/genome.fa

axel -q -n 4 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -o ${VCF}
readonly OUTPUT=1000genomes.phase3.decomposed.normalized.vcf.gz
readonly VT=$( which vt )

${VT} decompose -s "$VCF" \
    | ${VT} normalize -r "${REF}" - \
    | grep -v "##INFO=<ID=MEINFO" \
    | sed -e 's/MEINFO=[^;]*//' \
    | bgzip -f \
    > ${DBDIR}/${OUTPUT}

cd ${DBDIR}
tabix -p vcf ${OUTPUT}
