#!/bin/bash
#
# Calls an R script to apply a model given a model file,
# an input, an output, and optionally a threshold

readonly DEFTHRESH=0.71
readonly EXECUTABLE_PATH=${USE_EXECUTABLE_PATH:-"/usr/local/bin"}

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ ! -f "$1" ] || [ ! -f "$2" ]
then
    >&2 echo "$0 - Apply an ensemble model (provided) to the VCF"
    >&2 echo "Usage: $0 modelfile input.vcf raw_output.vcf output.vcf [thresh=${DEFTHRESH}]"
    >&2 echo "Invoked as: $0 $1 $2 $3 $4 $5"
    exit 
fi

readonly MODEL="$1"
readonly INPUTVCF="$2"
readonly RAWOUTPUTVCF="$3"
readonly OUTPUTVCF="$4"
readonly THRESH=${5:-${DEFTHRESH}}

Rscript --vanilla "${EXECUTABLE_PATH}/analysis/filter_calls_by_model.R" "$MODEL" "$INPUTVCF" "$RAWOUTPUTVCF" "${EXECUTABLE_PATH}/analysis/" "$THRESH" 

python "${EXECUTABLE_PATH}/clean_indel_calls.py" "${RAWOUTPUTVCF}" \
    | grep -v "^##INFO=<ID=RepeatRefCount" \
    | grep -v "^##INFO=<ID=TumorTotalDepth" \
    | grep -v "^##INFO=<ID=TumorVarDepth" \
    | grep -v "^##INFO=<ID=TumorVAF" \
    | grep -v "^##INFO=<ID=NormalTotalDepth" \
    | grep -v "^##INFO=<ID=NormalVarDepth" \
    | grep -v "^##INFO=<ID=NormalVAF" \
    | grep -v "^##INFO=<ID=dbsnp_VP" \
    | grep -v "^##contig=<ID=" \
    | sed -e 's/t_vaf/VAF/g' \
        > "${OUTPUTVCF}"
bgzip -f "${OUTPUTVCF}"
tabix -p vcf "${OUTPUTVCF}.gz"
