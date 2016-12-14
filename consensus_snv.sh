#!/bin/bash
###
### generate consensus SNV callset
##

function usage {
    >&2 echo "usage: $0 -b [/path/to/broad.snv.vcf.gz] -d [dkfz] -m [muse] -s [sanger] -o [outfile: defaults to merged.vcf] -c [/path/to/cosmicCoding.vcf.gz] -n [/path/to/cosmicNonCoding.vcf.gz]"
    >&2 echo "       Generates consensus somatic snv VCFs for sample"
    >&2 echo "       Input VCFs must be bgzip-ed and tabix-ed. Cosmic variant lists are optional"
    exit 1
}

readonly EXECUTABLE_PATH=${USE_EXECUTABLE_PATH:-"/usr/local/bin"}
outfile=consensus.snv.vcf

while getopts "b:d:m:s:o:c:n:h" OPTION 
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
        n) readonly cosmic_non_coding=${OPTARG}
            ;;
        c) readonly cosmic_coding=${OPTARG}
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
for file in "$musefile" "$broadfile" "$dkfzfile" "$sangerfile" "$cosmic_coding" "$cosmic_noncoding"
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
## make sure optional files are present if provided
##
for cosmicfile in "$cosmic_coding" "$cosmic_non_coding" 
do
    if [[ ! -z "$cosmicfile" ]] && [[ ! -f "$cosmicfile" ]] 
    then
        >&2 echo "Optional cosmic file provided but not present: ${cosmicfile}"
        usage
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
if [[ ! -z "$cosmic_coding" ]] 
then
    dbsnp_args+=("${cosmic_coding}")
fi
if [[ ! -z "$cosmic_non_coding" ]] 
then
    dbsnp_args+=("${cosmic_non_coding}")
fi

"${EXECUTABLE_PATH}"/dbsnp_annotate_one.sh  "${dbsnp_args[@]}"

rm -f ${INTERMEDIATE}
rm -f ${INTERMEDIATE}.gz
rm -f ${INTERMEDIATE}.gz.tbi
