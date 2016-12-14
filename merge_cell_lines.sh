#!/bin/bash

mkdir -p merged
mkdir -p dbsnp_annotated
mkdir -p consensus_processed

module load python/2.7.10
module load python-packages/2
module load gcc/4.8.1 
module load bam-readcount/0.7.4
module load samtools/1.2
module load tabix
module load vcfanno

readonly DATESTRING=$( date +"%Y%m%d" )
readonly ROOT="/scratch2/users/jdursi/cell-lines/"

#donor_unique_id submitter_donor_id icgc_donor_id normal_wgs_alignment_gnos_repo normal_wgs_alignment_bam_file_name tumor_wgs_aliquot_id tumor_wgs_alignment_bam_file_name
#TCGA_MUT_BENCHMARK_4::G15511 G15511  https://gtrepo-bsc.annailabs.com/ b50ceff5cf365343616d724324dcb445.bam 519b8381-95d5-4fce-a90c-7576cce2110c 257235f3926b2be84e8a9e80acdfb345.bam
#TCGA_MUT_BENCHMARK_4::G15512 G15512  https://gtrepo-bsc.annailabs.com/ a144319a81bbf0fc74942ef93273516b.bam 3d2edf87-6ec5-4c9f-9212-e8a751cc33e8 2fa459a3f50e3dfc2103ee9c233db716.bam

# SNVs
# sample 3d2edf87-6ec5-4c9f-9212-e8a751cc33e8 - G15512
${ROOT}/scripts/merge-one-tumour-snv.sh \
    -b processed/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.broad.snv_mnv.vcf.gz \
    -d processed/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.dkfz.snv_mnv.vcf.gz \
    -m processed/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.muse.snv_mnv.vcf.gz \
    -s processed/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.sanger.snv_mnv.vcf.gz \
    -o merged/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.merged.${DATESTRING}.somatic.snv_mnv.vcf

${ROOT}/scripts/dbsnp_annotate_one.sh \
    merged/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.merged.${DATESTRING}.somatic.snv_mnv.vcf.gz \
    dbsnp_annotated/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.annotated.${DATESTRING}.somatic.snv_mnv.vcf \
    snv

${ROOT}/scripts/vaf_oxog_annotate_one_snv.sh \
    dbsnp_annotated/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.annotated.${DATESTRING}.somatic.snv_mnv.vcf.gz \

# sample 519b8381-95d5-4fce-a90c-7576cce2110c - G15511
${ROOT}/scripts/merge-one-tumour-snv.sh \
    -b processed/519b8381-95d5-4fce-a90c-7576cce2110c.broad.snv_mnv.vcf.gz \
    -d processed/519b8381-95d5-4fce-a90c-7576cce2110c.dkfz.snv_mnv.vcf.gz \
    -m processed/519b8381-95d5-4fce-a90c-7576cce2110c.muse.snv_mnv.vcf.gz \
    -s processed/519b8381-95d5-4fce-a90c-7576cce2110c.sanger.snv_mnv.vcf.gz \
    -o merged/519b8381-95d5-4fce-a90c-7576cce2110c.merged.${DATESTRING}.somatic.snv_mnv.vcf

${ROOT}/scripts/dbsnp_annotate_one.sh \
    merged/519b8381-95d5-4fce-a90c-7576cce2110c.merged.${DATESTRING}.somatic.snv_mnv.vcf.gz \
    dbsnp_annotated/519b8381-95d5-4fce-a90c-7576cce2110c.annotated.${DATESTRING}.somatic.snv_mnv.vcf \
    snv

${ROOT}/scripts/vaf_oxog_annotate_one_snv.sh \
    dbsnp_annotated/519b8381-95d5-4fce-a90c-7576cce2110c.annotated.${DATESTRING}.somatic.snv_mnv.vcf.gz

exit 1

# INDELs
# sample 3d2edf87-6ec5-4c9f-9212-e8a751cc33e8 - G15512
${ROOT}/scripts/merge-one-tumour-indel.sh \
    -b processed/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.broad.indel.vcf \
    -d processed/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.dkfz.indel.vcf \
    -s processed/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.sanger.indel.vcf \
    -m processed/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.smufin.indel.vcf \
    -o merged/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.merged.${DATESTRING}.somatic.indel.vcf

${ROOT}/scripts/dbsnp_annotate_one_indel.sh \
    merged/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.merged.${DATESTRING}.somatic.indel.vcf.gz \
    dbsnp_annotated/3d2edf87-6ec5-4c9f-9212-e8a751cc33e8.processed.${DATESTRING}.somatic.indel.vcf

# sample 519b8381-95d5-4fce-a90c-7576cce2110c - G15512
${ROOT}/scripts/merge-one-tumour-indel.sh \
    -b processed/519b8381-95d5-4fce-a90c-7576cce2110c.broad.indel.vcf \
    -d processed/519b8381-95d5-4fce-a90c-7576cce2110c.dkfz.indel.vcf \
    -s processed/519b8381-95d5-4fce-a90c-7576cce2110c.sanger.indel.vcf \
    -m processed/519b8381-95d5-4fce-a90c-7576cce2110c.smufin.indel.vcf \
    -o merged/519b8381-95d5-4fce-a90c-7576cce2110c.merged.${DATESTRING}.somatic.indel.vcf

${ROOT}/scripts/dbsnp_annotate_one_indel.sh \
    merged/519b8381-95d5-4fce-a90c-7576cce2110c.merged.${DATESTRING}.somatic.indel.vcf.gz \
    dbsnp_annotated/519b8381-95d5-4fce-a90c-7576cce2110c.processed.${DATESTRING}.somatic.indel.vcf
