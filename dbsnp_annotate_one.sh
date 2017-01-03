#!/bin/bash 

readonly DB_PATH=${USE_DB_PATH:-"/dbs/annotation_databases"}
readonly EXECUTABLE_PATH=${USE_EXECUTABLE_PATH:-"/usr/local/bin"}

function usage {
    >&2 echo "usage: $0 input.vcf.gz snv|indel output_file [/path/to/cosmic_coding] [/path/to/cosmic_noncoding]"
    >&2 echo "       annotates one merged tumor"
    exit 1
}

readonly input_file="$1"
readonly variant_type="$2"
readonly output_file="$3"
readonly cosmic_coding="$4"
readonly cosmic_non_coding="$5"

if [[ -z "$input_file" ]] || [[ -z "$output_file" ]] 
then
    >&2 echo "argument missing"
    >&2 echo "invocation: $0 \"$1\" \"$2\""
    usage
fi

if [[ ! -f "$input_file" ]]
then
    >&2 echo "file missing: ${input_file} not found"
    usage
fi

if [[ -z "$variant_type" ]]
then
    >&2 echo "Variant type missing"
    usage
fi

if [[ "$variant_type" != "snv" ]] && [[ "$variant_type" != "indel" ]]
then
    >&2 echo "Invalid variant type $variant_type"
    >&2 echo "Valid variant types: snv indel"
    usage
fi

if [[ ! -d "$DB_PATH" ]]
then
    >&2 echo "database directory missing: ${DB_PATH} not found"
    usage
fi

if [[ -f "${output_file}.gz" ]] && [[ "${output_file}.gz" -nt "$input_file" ]]
then
    >&2 echo "$0: ${output_file} exists and is newer than inputs; cowardly refusing to overwrite."
    exit 1
fi

function fix_vcfanno_header {
    sed \
        -e 's/^##INFO=<ID=1000genomes_AF.*$/##INFO=<ID=1000genomes_AF,Number=1,Type=Float,Description="Thousand Genomes phase 3 occurance fraction if found: ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz">/' \
        -e 's/^##INFO=<ID=1000genomes_ID.*$/##INFO=<ID=1000genomes_ID,Number=1,Type=String,Description="Thousand Genomes phase 3 ID if found: ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz">/' \
        -e 's/^##INFO=<ID=cosmic,.*$/##INFO=<ID=cosmic,Number=1,Type=String,Description="(first) cosmic ID if found, COSMICv76">/' \
        -e 's/^##INFO=<ID=dbsnp,.*$/##INFO=<ID=dbsnp,Number=1,Type=String,Description="(first) dbSNP ID if found, build 147, All_20160408.vcf.gz">/' \
        -e 's/^##INFO=<ID=repeat_masker,.*$/##INFO=<ID=repeat_masker,Number=1,Type=String,Description="Repeat masker region if in one">/' 
}

readonly DBSNP_ANNOTATIONS=/tmp/dbsnp_annotations.conf
rm -f $DBSNP_ANNOTATIONS
cat > $DBSNP_ANNOTATIONS <<EOF
[[annotation]]
file = "${DB_PATH}/All_20160601.vcf.gz"
fields = ["ID", "VP"]
names = ["dbsnp", "dbsnp_VP"]
ops = ["first", "first"]

[[annotation]]
file = "${DB_PATH}/hg19.rmsk.bed.gz"
columns = [4]
names = ["repeat_masker"]
ops = ["first"]

[[annotation]]
file = "${DB_PATH}/1000genomes.phase3.decomposed.normalized.vcf.gz"
fields = ["ID", "AF"]
names = ["1000genomes_ID", "1000genomes_AF"]
ops = ["self", "self"]

EOF

for cosmicfile in "$cosmic_coding" "$cosmic_non_coding"
do
    if [[ -f "$cosmicfile" ]]
    then
        cat >> $DBSNP_ANNOTATIONS <<EOF
[[annotation]]
file = "${cosmicfile}"
fields = ["ID"]
names = ["cosmic"]
ops = ["first"]

EOF
    fi
done

function cleancalls {
    if [[ $variant_type == "snv" ]]
    then
        "${EXECUTABLE_PATH}/clean_snv_calls.py" 
    else
        cat
    fi
}

vcfanno $DBSNP_ANNOTATIONS "${input_file}" \
    | fix_vcfanno_header \
    | cleancalls \
    > "${output_file}"

bgzip "${output_file}"
tabix -p vcf "${output_file}.gz"

rm $DBSNP_ANNOTATIONS
