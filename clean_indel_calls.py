#!/usr/bin/env python
from __future__ import print_function
import argparse
import vcf  
import vcf.parser
import sys

def round_dig(value, ndigits=2):
    scale=10**ndigits
    return int(value*scale)*1./scale

def main():
    parser = argparse.ArgumentParser(description='Annotate merged vcf with VAF information where available')
    parser.add_argument('inputvcf', type=argparse.FileType('r'), default=sys.stdin, nargs='?', help="Merged and annotated VCF file")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help="Specify output file (default:stdout)")
    args = parser.parse_args()

    reader = vcf.Reader(args.inputvcf)
    reader.metadata['reference']='ftp://ftp.sanger.ac.uk/pub/project/PanCancer/genome.fa.gz'
    reader.infos['dbsnp_somatic'] = vcf.parser._Info(id='dbsnp_somatic', num=None, type='Flag', desc='Known-somatic dbSNP variant', source=None, version=None)
    reader.infos['t_vaf'] = vcf.parser._Info(id='t_vaf', num=1, type='Float', desc='VAF in tumor from sga', source=None, version=None)
    reader.infos['n_vaf'] = vcf.parser._Info(id='n_vaf', num=1, type='Float', desc='VAF in normal from sga', source=None, version=None)
    reader.infos['t_alt_count'] = vcf.parser._Info(id='t_alt_count', num=1, type='Integer', desc='Tumor alt count from sga', source=None, version=None)
    reader.infos['t_ref_count'] = vcf.parser._Info(id='t_ref_count', num=1, type='Integer', desc='Tumor ref count from sga if non-zero', source=None, version=None)
    reader.infos['n_alt_count'] = vcf.parser._Info(id='n_alt_count', num=1, type='Integer', desc='Normal alt count from sga if non-zero', source=None, version=None)
    reader.infos['n_ref_count'] = vcf.parser._Info(id='n_ref_count', num=1, type='Integer', desc='Normal ref count from sga if non-zero', source=None, version=None)
    reader.infos['model_score'] = vcf.parser._Info(id='model_score', num=1, type='Float', desc='consensus model score, 0-1', source=None, version=None)
    reader.filters['LOWSUPPORT'] = vcf.parser._Filter(id='LOWSUPPORT', desc='Insufficient support in consensus model')
    writer = vcf.Writer(args.output, reader)

    for record in reader:
        new_info = {}
        # skip broad private calls
        if record.INFO['Callers'] == ['broad']:
            continue

        # skip calls that are defined to be SVs
        varlen = abs(len(record.REF) - len(record.ALT[0]))
        if varlen >= 100:
            continue

        # copy some records over directly
        for item in ['dbsnp', 'cosmic', 'Callers', 'NumCallers', 'repeat_masker', '1000genomes_AF', '1000genomes_ID']:
            if item in record.INFO:
                new_info[item] = record.INFO[item]

        if 'model_score' in record.INFO:
            new_info['model_score'] = round_dig(float(record.INFO['model_score']),3)

        if 'dbsnp_VP' in record.INFO:
            qualbyte = int(record.INFO['dbsnp_VP'],16) & 255
            somatic = (qualbyte & 2**5) > 0
            if somatic:
                new_info['dbsnp_somatic'] = True

        if ('TumorVarDepth' in record.INFO) and (record.INFO['TumorVarDepth'] > 0):
            new_info['t_vaf'] = record.INFO['TumorVAF']
            new_info['t_alt_count'] = record.INFO['TumorVarDepth']
            new_info['t_ref_count'] = record.INFO['TumorTotalDepth']-record.INFO['TumorVarDepth']

        if ('NormalVarDepth' in record.INFO) and (record.INFO['NormalVarDepth'] > 0):
            new_info['n_vaf'] = record.INFO['NormalVAF']
            new_info['n_alt_count'] = record.INFO['NormalVarDepth']
            new_info['n_ref_count'] = record.INFO['NormalTotalDepth']-record.INFO['NormalVarDepth']

        record.INFO = new_info
        writer.write_record(record)

    return 0

if __name__ == "__main__":
    sys.exit(main())
