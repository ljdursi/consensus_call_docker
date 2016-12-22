#!/bin/bash -l
#
# Calls an R script to apply a model given a model file,
# an input, an output, and optionally a threshold

module purge 
module use /.mounts/labs/simpsonlab/modules/
module load R/3.3.0
module load python-packages/2
module load tabix/0.2.6

readonly DEFTHRESH=0.725

if [ $# -eq 0 ] || [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ] || [ ! -f "$1" ] || [ ! -f "$2" ]
then
    >&2 echo "$0 - Apply an ensemble model (provided) to the VCF"
    >&2 echo "Usage: $0 modelfile input.vcf raw_output.vcf output.vcf [thresh=${DEFTHRESH}]"
    >&2 echo "Invoked as: $0 $1 $2 $3 $4 $5"
    exit 
fi

readonly MODEL=$1
readonly INPUTVCF=$2
readonly RAWOUTPUTVCF=$3
readonly OUTPUTVCF=$4
readonly THRESH=${5:-${DEFTHRESH}}

readonly VARIANT="indel"

Rscript --vanilla scripts/filter_calls_by_model.R $MODEL $INPUTVCF $RAWOUTPUTVCF $THRESH

ID=$( basename $INPUTVCF | cut -f 1 -d . )
PCAWG1ID=$( ./scripts/pancanid_to_pcawg1id.sh ${ID} )
readonly STARS=/oicr/data/pancanxfer/consensus/annotations/star/PAWG_QC_Summary_of_Measures.tsv
readonly classification_maf=/oicr/data/pancanxfer/consensus/annotations/final_classifications/${VARIANT}/${ID}.${VARIANT}.maf
readonly normalpanel=/oicr/data/pancanxfer/consensus/filters/panel_of_normals/${VARIANT}/${ID}.${VARIANT}.maf
readonly TIN=/oicr/data/pancanxfer/consensus/annotations/TiN/release_may2016.v1.1.TiN__donor.TiNsorted.20Jul2016.tsv
readonly VALIDATION_FILE=/oicr/data/pancanxfer/validation/vcfs/quality-filtered/${PCAWG1ID}.${VARIANT}.vcf
readonly RELEASE=/oicr/data/pancanxfer/consensus/annotations/release/by_tumour_id.tsv
SEX=$( ./scripts/sex-from-id.sh ${ID} )

python ./scripts/clean_indel_calls.py ${RAWOUTPUTVCF} \
    | python ./scripts/apply_validation_calls.py ${VALIDATION_FILE} \
    | python ./scripts/filter_by_presence_in_maf.py -d "Presence in Panel of Normals" ${normalpanel} ${ID} NORMALPANEL \
    | python ./scripts/info_or_filter_from_MAF.py -a info -c Variant_Classification -d "Variant Classification" ${classification_maf} Variant_Classification \
    | ./scripts/annotate_vcf_from_tsv_column.sh -c 11 -i ${ID} -n TumourInNormalEstimate -t ${TIN} \
    | ./scripts/annotate_vcf_from_tsv_column.sh -c 3 -i ${ID} -n BAMQCStars -t ${RELEASE} \
    | ./scripts/annotate_vcf_from_tsv_column.sh -c 2 -i ${ID} -n ContEST -t ${RELEASE} \
    | python ./scripts/apply_sex_filter.py -s "${SEX}" \
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
        > ${OUTPUTVCF}
bgzip -f ${OUTPUTVCF}
tabix -p vcf ${OUTPUTVCF}.gz
