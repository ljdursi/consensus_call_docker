#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import vcf.parser
import sys

def main():
    parser = argparse.ArgumentParser(description='Fix dbsnp VP calls and add OXOG filter')
    parser.add_argument('inputvcf', type=argparse.FileType('r'), default=sys.stdin, nargs='?', help="Merged and annotated VCF file")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    args = parser.parse_args()

    reader = vcf.Reader(args.inputvcf)
    reader.infos['dbsnp_somatic'] = vcf.parser._Info(id='dbsnp_somatic', num=None, type='Flag', desc='Known-somatic dbSNP variant', source=None, version=None)
    del reader.infos['OXOG_Fail']
    reader.filters['OXOGFAIL'] = vcf.parser._Filter(id='OXOGFAIL', desc="Failed OXOG oxidative artifact filter")
    reader.metadata['reference']='ftp://ftp.sanger.ac.uk/pub/project/PanCancer/genome.fa.gz'
    writer = vcf.Writer(args.output, reader)

    for record in reader:
        new_info = {}

        # copy some records over directly
        for item in ['VAF', 't_alt_count', 't_ref_count']:
            if item in record.INFO and record.INFO[item] > 0:
                new_info[item] = record.INFO[item]

        for item in ['dbsnp', 'cosmic', 'Callers', 'NumCallers', 'repeat_masker', '1000genomes_AF', '1000genomes_ID']:
            if item in record.INFO:
                new_info[item] = record.INFO[item]

        if 'dbsnp_VP' in record.INFO:
            qualbyte = int(record.INFO['dbsnp_VP'],16) & 255
            somatic = (qualbyte & 2**5) > 0
            if somatic:
                new_info['dbsnp_somatic'] = True

        if 'OXOG_Fail' in record.INFO:
            if record.INFO['OXOG_Fail'] == 'True':
                if record.FILTER is None:
                    record.FILTER = ['OXOGFAIL']
                else:
                    record.FILTER.append('OXOGFAIL')

        record.INFO = new_info
        writer.write_record(record)

    return 0

if __name__ == "__main__":
    sys.exit(main())
