#!/bin/bash -x
module load vcfanno
module load tabix/0.2.6
module load python-packages/2

function usage {
    >&2 echo "usage: $0 input-file"
    >&2 echo "       annotates one merged snv tumor vcf with OxoG status"
    exit 1
}

readonly VARIANT=snv_mnv
readonly INDIR=dbsnp_annotated/${VARIANT}

readonly INPUT_FILE=$1
if [[ -z "$INPUT_FILE" ]] || [[ ! -f "$INPUT_FILE" ]]
then
    >&2 echo "input file missing $INPUT_FILE"
    usage
fi

readonly ID=$( basename $INPUT_FILE | cut -f 1 -d . )
if [[ -z "$ID" ]]
then
    >&2 echo "no ID"
    usage
fi

readonly OUTDIR=annotated/${VARIANT}
mkdir -p $OUTDIR

readonly DATE=$( date +%Y%m%d )
readonly OUTPUT_FILE=${OUTDIR}/${ID}.consensus.${DATE}.somatic.${VARIANT}.vcf
sed -e "s/@@SAMPLE@@/${ID}/" -e "s/@@PATH@@/processed/" annotation/vaf_oxog.annotations.conf.template > annotation/vaf.${ID}.conf

if [[ -f $OUTPUT_FILE ]] && [[ $outfile -nt $INPUT_FILE ]]
then
    >&2 echo "$0: ${OUTPUT_FILE} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

vcfanno -p 1 annotation/vaf.${ID}.conf ${INPUT_FILE} 2> /dev/null \
    | python ./scripts/clean_snv_calls.py -o ${OUTPUT_FILE} 
rm annotation/vaf.${ID}.conf

bgzip -f ${OUTPUT_FILE}
tabix -p vcf ${OUTPUT_FILE}.gz
