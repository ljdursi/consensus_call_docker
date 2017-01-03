#!/bin/bash

readonly COMMAND=$1
readonly DIRECTORY=$2

readonly EXECUTABLE_PATH=${USE_EXECUTABLE_PATH:-"/usr/local/bin"}

function usage {
    echo >&2 "$0: downloads and installs necessary databases for consensus calling"
    echo >&2 "    options: "
    echo >&2 "       reference /directory/name: "
    echo >&2 "       annotations /directory/name/: annotation files (dbsnp, 1000 genomes, repeat masker)"
    echo >&2 "       Must download reference before annotations."
    exit 1
}

###
### Make sure options are valid 
###
if [[ -z "${COMMAND}" ]] || [[ -z "${DIRECTORY}" ]] 
then
    echo >&2 "Missing arguments."
    usage
fi

if [[ "${COMMAND}" != "reference" ]] && [[ "${COMMAND}" != "annotations" ]]
then
    echo >&2 "Invalid sub-command ${COMMAND}."
    usage
fi

if [[ ! -d "${DIRECTORY}" ]]
then
    echo >&2 "Invalid directory ${DIRECTORY}."
    usage
fi

###
### download/install pancan standard reference
###

readonly REFERENCE_URL=ftp://ftp.sanger.ac.uk/pub/project/PanCancer/genome.fa.gz
readonly REFERENCE_FILE=genome.fa

if [[ "${COMMAND}" == "reference" ]] 
then
    wget -nv "${REFERENCE_URL}" -O - \
        | zcat \
        | bgzip > "${DIRECTORY}/${REFERENCE_FILE}.gz"
    samtools faidx "${DIRECTORY}/${REFERENCE_FILE}.gz"
fi

###
### download/install annotation files
###
readonly PATH_TO_REFERENCE="${DIRECTORY}/reference/${REFERENCE_FILE}.gz"
if [[ -z "${PATH_TO_REFERENCE}" ]] || [[ ! -f "${PATH_TO_REFERENCE}" ]] 
then
    echo >&2 "No reference found! ${PATH_TO_REFERENCE}"
fi

readonly INTERMEDIATE="${DIRECTORY}/annotation_databases/intermediate.vcf.gz"
rm -f "${INTERMEDIATE}"

readonly DBSNP_URL="ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/All_20160601.vcf.gz"
readonly RMSK_URL="http://people.virginia.edu/~arq5x/files/gemini/annotations/hg19.rmsk.bed.gz"
readonly KGENOMES_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz"

readonly DBSNP="${DIRECTORY}/annotation_databases/All_20160601.vcf"
readonly RMSK="${DIRECTORY}/annotation_databases/hg19.rmsk.bed"
readonly KGENOMES="${DIRECTORY}/annotation_databases/1000genomes.phase3.decomposed.normalized.vcf"

function download_and_normalize {
    local url="$1"
    local intermediate="$2"
    local output="$3"
    local reference="$4"

    wget -nv "$url" -O "$intermediate" 
    vt decompose -s "$intermediate" \
        | vt normalize -r "$reference" - \
        | grep -v "##INFO=<ID=MEINFO" \
        | sed -e 's/MEINFO=[^;]*//' \
        | bgzip -f > "$output" 
    tabix -p vcf "$output" 
    rm -f "$intermediate"
}

if [[ "${COMMAND}" == "annotations" ]] 
then
    # repeat masker
    wget -nv "${RMSK_URL}" -O - \
        | zcat \
        | bgzip \
        > "${RMSK}.gz"
    tabix -p bed "${RMSK}.gz"

    # dbsnp 
    download_and_normalize "${DBSNP_URL}" "${INTERMEDIATE}" "${DBSNP}.gz" "${PATH_TO_REFERENCE}"

    # 1000 genomes
    download_and_normalize "${KGENOMES_URL}" "${INTERMEDIATE}" "${KGENOMES}.gz" "${PATH_TO_REFERENCE}"
fi
