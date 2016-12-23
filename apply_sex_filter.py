#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import vcf.parser
import sys

def main():
    parser = argparse.ArgumentParser(description='Filter Y-chromosome calls if sex is female')
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=sys.stdin, help="input VCF file")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    parser.add_argument('-s', '--sex', type=str, default="male", help="Apply filter if this flag is exactly 'female' (default:male)")
    args = parser.parse_args()

    apply_filter = args.sex == "female"
    filtername = 'SEXF'

    reader = vcf.Reader(args.input)
    if apply_filter:
        reader.filters[filtername] = vcf.parser._Filter(id=filtername, desc='Likely artifact or call in PAR region: Y-chromosome variant in female donor')
    writer = vcf.Writer(args.output, reader)

    for record in reader:
        if (apply_filter and record.CHROM in ['Y', 'chrY']):
            if not record.FILTER:
                record.FILTER = [filtername]
            else:
                record.FILTER += [filtername]
        writer.write_record(record)
    return 0

if __name__ == "__main__":
    sys.exit(main())
