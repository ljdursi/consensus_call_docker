#!/bin/bash

readonly UNDEF=None
COLUMN=$UNDEF
ID=$UNDEF
TSV=$UNDEF
NAME=$UNDEF

function usage {
    >&2 echo "Usage: $0 -t tsvfile -n annotation_name -c column -i id_to_search_for"
    >&2 echo "       Adds an annotation to a VCF header (from stdin) based on presence of id in file tsvfile"
    >&2 echo "       Annotation is of the form ##annotation_name=[value in column c of row containing id_to_search_for"
    exit 1
}


while getopts ":c:i:n:t:" o
do
    case "${o}" in
        c) COLUMN=${OPTARG} ;;
        i) ID=${OPTARG} ;;
        n) NAME=${OPTARG} ;;
        t) TSV=${OPTARG} ;;
    esac
done

if [[ "$COLUMN" == "$UNDEF" ]] || [[ $ID == $UNDEF ]] || [[ $NAME == $UNDEF ]] || [[ $TSV == $UNDEF ]]
then
    usage
fi

key=$( grep $ID $TSV | cut -f$COLUMN -d$'\t')
echo "##fileformat=VCFv4.1"
if [[ "$key" != "" ]]
then
    echo "##${NAME}=$key"
fi
tail -n +2
