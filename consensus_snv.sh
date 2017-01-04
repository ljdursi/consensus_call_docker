#!/bin/bash
###
### generate consensus SNV callset
###

function usage {
    >&2 echo "usage: $0 -b [/path/to/broad.snv.vcf.gz] -d [dkfz] -m [muse] -s [sanger] -o [outfile: defaults to merged.vcf]"
    >&2 echo "       Generates consensus somatic snv VCFs for sample"
    >&2 echo "       Input VCFs must be bgzip-ed and tabix-ed. Cosmic variant lists are optional"
    exit 1
}

readonly EXECUTABLE_PATH=${USE_EXECUTABLE_PATH:-"/usr/local/bin"}
outfile=consensus.snv.vcf

while getopts "b:d:m:s:o:h" OPTION 
do
    case $OPTION in
        b) readonly broadfile=${OPTARG}
            ;;
        d) readonly dkfzfile=${OPTARG}
            ;;
        m) readonly musefile=${OPTARG}
            ;;
        s) readonly sangerfile=${OPTARG}
            ;;
        o) outfile=${OPTARG}
            ;;
        h) usage
            ;;
    esac
done

##
## make sure required arguments are given
##
if [[ -z "$musefile" ]] || [[ -z "$broadfile" ]] || [[ -z "$dkfzfile" ]] || [[ -z "$sangerfile" ]]
then
    >&2 echo "required argument missing: need muse (-m), broad (-b), dkfz (-d) and sanger (-s) files"
    usage
fi

if [[ -z "${outfile}" ]]
then
    >&2 echo "Invalid empty output filename"
    usage
fi

##
## make sure the files (look to be) bgzipped and have .tbi files
##
for file in "$musefile" "$broadfile" "$dkfzfile" "$sangerfile"
do
    if [[ ! -z "$file" ]] 
    then
        if [[ $file != *.gz ]] || [[ ! -f "${file}.tbi" ]]
        then
            >&2 echo "Input VCF files must be bgziped and tabixed."
            usage
        fi
    fi
done

##
## Merge SNV calls
##
readonly INTERMEDIATE=/tmp/merged.snv.vcf
"${EXECUTABLE_PATH}"/merge-one-tumour-snv.sh \
    -b "${broadfile}" -d "${dkfzfile}" -m "${musefile}" -s "${sangerfile}" -o "$INTERMEDIATE"

##
## annotate with dbsnp, 1kgenomes, repeat_masker, and cosmic if provided
##
dbsnp_args=("${INTERMEDIATE}.gz" "snv" "${outfile}")
"${EXECUTABLE_PATH}"/dbsnp_annotate_one.sh  "${dbsnp_args[@]}"

rm -f "${INTERMEDIATE}"
rm -f "${INTERMEDIATE}.gz"
rm -f "${INTERMEDIATE}.gz.tbi"
