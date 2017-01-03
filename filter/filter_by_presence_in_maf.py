#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import vcf.parser
import csv
import sys

def variant_tuple(sample, record):
    return (sample, record.CHROM, record.POS)

def get_entries_from_MAF(maf):
    entries = set()
    try:
        with open(maf, 'r') as maffile:
            mafreader = csv.DictReader(maffile, delimiter='\t')
            for record in mafreader:
                items = ['tumor_aliquot_id', 'Chromosome', 'Start_position' ]
                if not all([item in record for item in items]):
                    continue
                variant=(record['tumor_aliquot_id'], record['Chromosome'], int(record['Start_position']))
                entries.add(variant)
                if 'End_position' in record:
                    variant=(record['tumor_aliquot_id'], record['Chromosome'], int(record['End_position']))
                    entries.add(variant)
    except:
        pass
    return entries

def main():
    parser = argparse.ArgumentParser(description='Fix dbsnp VP calls and add OXOG filter')
    parser.add_argument('MAF', type=str, help="MAF file for filtering")
    parser.add_argument('sample', type=str, help="tumour aliquot id")
    parser.add_argument('filtername', type=str, help="Filter name to apply")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin, help="Merged and annotated VCF file (default: stdin)")
    parser.add_argument('-d', '--desc', type=str, default="", help="Description of filter")
    parser.add_argument('-n', '--info', action='store_true', help="Add info flag rather than filter")
    args = parser.parse_args()

    reader = vcf.Reader(args.input)
    if args.info:
        reader.infos[args.filtername] = vcf.parser._Info(id=args.filtername, num=0, type='Flag', desc=args.desc, source=None, version=None)
    else:
        reader.filters[args.filtername] = vcf.parser._Filter(id=args.filtername, desc=args.desc)
    writer = vcf.Writer(args.output, reader)

    entries = get_entries_from_MAF(args.MAF)
    for record in reader:
        variants = [variant_tuple(args.sample, record)]
        if len(record.ALT[0]) != len(record.REF):
            variants += [(args.sample, record.CHROM, record.POS+abs(len(record.ALT[0])-len(record.REF)))]
        for variant in variants:
            if variant in entries:
                if not args.info:
                    if not record.FILTER:
                        record.FILTER = [args.filtername]
                    elif args.filtername not in record.FILTER:
                        record.FILTER = record.FILTER + [args.filtername]
                else:
                    record.INFO[args.filtername]=True
        writer.write_record(record)

    return 0

if __name__ == "__main__":
    sys.exit(main())
