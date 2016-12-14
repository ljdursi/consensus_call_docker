#!/usr/bin/env cwl-runner

class: CommandLineTool
id: pcawg-merge-annotate
label: PCAWG pre-merge annotation
cwlVersion: v1.0
dct:creator:
  '@id': http://orcid.org/0000-0002-4697-798X
  foaf:name: Jonathan Dursi
  foaf:mbox: mailto:jonathan@dursi.ca
requirements:
- class: DockerRequirement
  dockerPull: quay.io/ljdursi/pcawg-merge-annotate:1.0.0
#hints:
#  - class: ResourceRequirement
#    coresMin: 1
#    ramMin: 4092
#    outdirMin: 512000
#    description: "the process requires at least 4G of RAM"

inputs:
  output: string
  variant_type:
    type: string
    default: SNV
    inputBinding:
      position: 1

    doc: Annotate SNV (SNV) or indel (indel) variants
  input_vcf:
    type: File
    format: http://edamontology.org/format_3016
    inputBinding:
      position: 2

    doc: The VCF file to be annotated
  normal_bam:
    type: File
    secondaryFiles: .bai
    format: http://edamontology.org/format_2572
    inputBinding:
      position: 3

    doc: The normal BAM file
  tumour_bam:
    type: File
    secondaryFiles: .bai
    format: http://edamontology.org/format_2572
    inputBinding:
      position: 4

    doc: The tumour BAM file
stdout: $(inputs.output)

outputs:
  annotated_vcf:
    type: File
    format: http://edamontology.org/format_3016
    outputBinding:
      glob: $(inputs.output)
    doc: A zip file that contains the HTML report and various graphics.
baseCommand: []
doc: |
  A Docker container for the uniform annotation as part of the merging process.
