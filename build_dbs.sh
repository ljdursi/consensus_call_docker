#!/bin/bash

readonly COMMAND=$1
readonly DIRECTORY=$2

readonly EXECUTABLE_PATH=${USE_EXECUTABLE_PATH:-"/usr/local/bin"}

function usage {
    echo >&2 "$0: downloads and installs necessary databases for consensus calling"
    echo >&2 "    options: "
    echo >&2 "       reference /directory/name: "
    echo >&2 "       annotations /directory/name/: annotation files (dbsnp, 1000 genomes, repeat masker)"
    echo >&2 "       cosmic /directory/name/: Cosmic VCF files (coding, noncoding: must have password)"
    echo >&2 "       Must download reference before annotations, cosmic."
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

if [[ "${COMMAND}" != "reference" ]] && [[ "${COMMAND}" != "annotations" ]] && [[ "${COMMAND}" != "cosmic" ]]
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
readonly PATH_TO_REFERENCE="${DIRECTORY}/reference/${REFERENCE_FILE}.gz"

mkdir -p "${DIRECTORY}/reference"

if [[ "${COMMAND}" == "reference" ]] 
then
    readonly TMP_REF="${DIRECTORY}/reference/TMP.fa"
    wget -nv "${REFERENCE_URL}" -O "${TMP_REF}.gz"
    gunzip "${TMP_REF}.gz"
    bgzip "${TMP_REF}" > "${PATH_TO_REFERENCE}" \
        && rm -f "${TMP_REF}"
    samtools faidx "${PATH_TO_REFERENCE}"
fi

###
### download/install annotation files of one sort or another
###

if [[ -z "${PATH_TO_REFERENCE}" ]] || [[ ! -f "${PATH_TO_REFERENCE}" ]] 
then
    echo >&2 "No reference found! ${PATH_TO_REFERENCE}"
fi

mkdir -p "${DIRECTORY}/annotation_databases"
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
    local intermediate="$2" local output="$3"
    local reference="$4"
    local method="$5"

    if [[ "$method" == "sftp" ]]
    then
        sftp "$url" "$intermediate"
    else
        wget -nv "$url" -O "$intermediate" 
    fi

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
    download_and_normalize "${DBSNP_URL}" "${INTERMEDIATE}" "${DBSNP}.gz" "${PATH_TO_REFERENCE}" "wget"

    # 1000 genomes
    download_and_normalize "${KGENOMES_URL}" "${INTERMEDIATE}" "${KGENOMES}.gz" "${PATH_TO_REFERENCE}" "wget"
fi

###
### Cosmic annotations (separate because optional for SNVs, and requires registration)
###

if [[ "${COMMAND}" == "cosmic" ]]
then
    # get Cosmic coding and noncoding variants
    echo "Downloading COSMIC database: you will be asked for username (email address) and password."
    echo "You can register for COSMIC access at: https://cancer.sanger.ac.uk/cosmic/register"
    echo "Make sure docker is run with -it (interactive, provide terminal)"
    echo ""
    echo "Enter username (email address):"
    read -r emailaddr

    readonly COSMICADDR=sftp-cancer.sanger.ac.uk
    readonly CODINGPATH=/files/grch37/cosmic/v77/VCF/CosmicCodingMuts.vcf.gz
    readonly NONCODINGPATH=/files/grch37/cosmic/v77/VCF/CosmicNonCodingVariants.vcf.gz

    readonly CODINGPATH=${DIRECTORY}/annotation_databases/CosmicCodingMuts.vcf.gz
    readonly NONCODINGPATH=${DIRECTORY}/annotation_databases/CosmicNonCodingVariants.vcf.gz

    # coding
    download_and_normalize "\"${emailaddr}\"@${COSMICADDR}:${CODING}" "${INTERMEDIATE}" "${CODINGPATH}" "${PATH_TO_REFERENCE}" "sftp"

    # noncoding
    download_and_normalize "\"${emailaddr}\"@${COSMICADDR}:${NONCODING}" "${INTERMEDIATE}" "${NONCODINGPATH}" "${PATH_TO_REFERENCE}" "sftp"
fi
