#!/bin/bash
## merges SNV files together
## requires mergevcf, bgzip, tabix

function usage {
    >&2 echo "usage: $0 -b [/path/to/broad.snv.vcf.gz] -d [dkfz] -m [muse] -s [sanger] -o [outfile: defaults to merged.vcf]"
    >&2 echo "       merges snv VCFs for sample"
    >&2 echo "       Input VCFs must be bgzip-ed and tabix-ed"
    exit 1
}

function add_header_and_sort {
    local file=$1
    local header=$2
    local grep_or_zgrep=grep
    if [[ $file == *.gz ]]
    then
        local grep_or_zgrep=zgrep
    fi

    ${grep_or_zgrep} "^##" "$file"
    if [[ ! -z "$header" ]] 
    then
        echo $header
    fi
    ${grep_or_zgrep} -v "^##" "$file" \
        | sort -k1,1 -k2,2n
}

function cleanup {
    local file=$1
    add_header_and_sort "$file" \
        | grep -v '=$' \
        | sed -e 's/Tier[0-9]/PASS/' 
}


outfile=merged.vcf

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
        h) usage
            ;;
        o) outfile=${OPTARG}
            ;;
    esac
done

if [ -z "$musefile" ] || [ -z "$broadfile" ] || [ -z "$dkfzfile" ] || [ -z "$sangerfile" ]
then
    >&2 echo "required argument missing: need muse (-m), broad (-b), dkfz (-d) and sanger (-s) files"
    usage
fi

if [ -z "${outfile}" ]
then
    >&2 echo "Invalid empty output filename"
    usage
fi

if [ ! -f "$musefile" ] || [ ! -f "$broadfile" ] || [ ! -f "$dkfzfile" ] || [ ! -f "$sangerfile" ]
then
    >&2 echo "files missing: one of ${musefile} ${broadfile} ${dkfzfile} ${sangerfile} not found"
    usage
fi

for file in "$musefile" "$broadfile" "$dkfzfile" "$sangerfile"
do
    if [[ $file != *.gz ]] 
    then
        >&2 echo "Input VCF files must be bgziped and tabixed."
        usage
    fi
done

newest=$musefile
if [[ $broadfile -nt $newest ]]; then newest=$broadfile; fi
if [[ $dkfzfile -nt $newest ]]; then newest=$dkfzfile; fi
if [[ $sangerfile -nt $newest ]]; then newest=$sangerfile; fi

if [[ -f "${outfile}.gz" ]] && [[ "${outfile}.gz" -nt "${newest}" ]]
then
    >&2 echo "$0: ${outfile} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

##
## merge input files into one VCF
##

readonly MERGEDFILE=/tmp/merged.vcf
mergevcf -l broad,dkfz,muse,sanger \
        <( cleanup "${broadfile}" ) \
        <( cleanup "${dkfzfile}" ) \
        <( cleanup "${musefile}" ) \
        <( cleanup "${sangerfile}" )\
        --ncallers --mincallers 2 \
    | grep -v "Callers=muse;" \
    > ${MERGEDFILE}

add_header_and_sort ${MERGEDFILE} \
    | bgzip > ${MERGEDFILE}.gz
tabix -p vcf ${MERGEDFILE}.gz

##
## annotate with a filter if there were OXOG info fields in any of the input files
##

# create vcfanno config file to run vcf against
readonly OXOGCONF=/tmp/oxog.conf

rm -f ${OXOGCONF}
touch ${OXOGCONF}

for file in "${broadfile}" "${dkfzfile}" "${musefile}" "${sangerfile}"
do
    cat >> ${OXOGCONF} <<EOF
[[annotation]]
file = "$file"
fields = ["OXOG_Fail"]
names = ["OXOG_Fail"]
ops = ["self"]

EOF
done

vcfanno -p 1 ${OXOGCONF} ${MERGEDFILE}.gz 2> /dev/null \
    > ${outfile} 
bgzip -f ${outfile}
tabix -p vcf ${outfile}.gz

rm ${OXOGCONF}
rm ${MERGEDFILE}
rm ${MERGEDFILE}.gz
rm ${MERGEDFILE}.gz.tbi
