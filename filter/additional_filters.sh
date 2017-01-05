#!/bin/bash

readonly COMMAND=$1
readonly EXECUTABLE_PATH=${USE_EXECUTABLE_PATH:-"/usr/local/bin"}/filter

function usage {
    echo >&2 "$0: Adds additional filters or annotations post-consensus."
    echo >&2 "       presence_in_maf [--info] -d description /path/to/input.maf sample_ID filter_name -i /path/to/input.vcf -o output.vcf: adds filter/info if variant present in input.maf"
    echo >&2 "       column_of_maf [--info] -c column_name -d description /path/to/input.maf sample_ID filter_name -i /path/to/input.vcf -o output.vcf: adds filter/info to variant from column in input.maf"
    echo >&2 "       sex -s {male|female} /path/to/input.vcf: filters Y-chrom calls if female"
    echo >&2 "       header_from_tsv -n header_name -I sample_ID -c column_number -t /path/to/input.tsv -i /path/to/input.vcf : adds annotation to VCF header entry column in input.tsv"
    exit 1
}

if [[ -z "${COMMAND}" ]] 
then
    echo >&2 "Missing arguments."
    usage
fi

case "${COMMAND}" in
    "presence_in_maf")
        python "${EXECUTABLE_PATH}"/filter_by_presence_in_maf.py "${@:2}"
        ;;
    "column_of_maf")
        python "${EXECUTABLE_PATH}"/filter_from_MAF_column.py "${@:2}"
        ;;
    "sex")
        python "${EXECUTABLE_PATH}"/apply_sex_filter.py "${@:2}"
        ;;
    "header_from_tsv")
        "${EXECUTABLE_PATH}"/vcf_header_from_tsv_column.sh "${@:2}"
        ;;
    *)
        echo >&2 "Invalid command ${COMMAND}."
        usage
        ;;
esac
