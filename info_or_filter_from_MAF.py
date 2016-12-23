#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import vcf.parser
import csv
import sys

# MAF has different convention than VCFs, removes common prefixes
def remove_common_prefix(ref, alt):
    prefix_len = 0
    for r,a in zip(ref, alt):
        if r != a:
            break
        prefix_len += 1
    newref = ref[prefix_len:]
    newalt = alt[prefix_len:]
    if newref == "":
        newref = "-"
    if newalt == "":
        newalt = "-"
    return newref, newalt

def variant_tuple(record, alt):
    newref, newalt = remove_common_prefix(record.REF, str(alt))
    pos = record.POS
    if newalt == "-":
        pos += 1
    return (record.CHROM, pos, newref, newalt)

def get_classification_dict_from_MAF(maf, column):
    classifications = {}
    try:
        with open(maf, 'r') as maffile:
            mafreader = csv.DictReader(maffile, delimiter='\t')
            for record in mafreader:
                items = ['Chromosome', 'Start_position', 'Reference_Allele', 'Tumor_Seq_Allele2']
                if column is not None:
                    items += [column]
                if not all([item in record for item in items]):
                    continue
                variant=(record['Chromosome'], int(record['Start_position']), record['Reference_Allele'], record['Tumor_Seq_Allele2'])
                if column is not None:
                    classifications[variant] = record[column]
                else:
                    classifications[variant] = 1
    except:
        pass
    return classifications

def main():
    parser = argparse.ArgumentParser(description='Add info or filter tag based on presence in MAF')
    parser.add_argument('MAF', type=str, help="MAF file for variant classification annotation")
    parser.add_argument('name', help='Name of info field or filter')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin, help="VCF file to be processed (default: stdin)")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Output file (default:stdout)")
    parser.add_argument('-a', '--action', choices=['info','filter'], default='info', help='add tag (info) or filter based on presence in MAF (default:info)')
    parser.add_argument('-c', '--column', type=str, help='column in MAF to use for info, if present')
    parser.add_argument('-d', '--description', default="", type=str, help='description of new info/filter field')
    args = parser.parse_args()

    reader = vcf.Reader(args.input)
    if args.action == "info":
        reader.infos[args.name] = vcf.parser._Info(id=args.name, num='1', type='String', desc=args.description, source=None, version=None)
    else:
        reader.filters[args.name] = vcf.parser._Filter(id=args.name, desc=args.description)
    writer = vcf.Writer(args.output, reader)
    classification_dict = get_classification_dict_from_MAF(args.MAF, args.column)

    for record in reader:
        assert len(record.ALT) == 1
        variant = variant_tuple(record, record.ALT[0])
        if variant in classification_dict:
            if args.action == "info":
                record.INFO[args.name] = classification_dict[variant]
            else:
                if not record.FILTER:
                    record.FILTER = [args.name]
                else:
                    record.FILTER += [args.name]

        writer.write_record(record)

    return 0

if __name__ == "__main__":
    sys.exit(main())
