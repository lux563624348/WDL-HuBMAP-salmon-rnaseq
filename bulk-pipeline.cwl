#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.2
label: bulk scRNA-seq pipeline using Salmon
requirements:
  ScatterFeatureRequirement: {}
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  fastq_dir:
    label: "Directory containing FASTQ files"
    type: Directory
  threads:
    label: "Number of threads for Salmon"
    type: int
    default: 1
  organism:
    type: string?
outputs:
  fastqc_dir:
    outputSource: fastqc/fastqc_dir
    type: Directory
    label: "Directory of FastQC output files, mirroring input directory structure"
  salmon_output:
      outputSource: salmon-bulk/output_dir
      type: Directory
      label: "Full output of `salmon quant`"
  expression_matrix:
    outputSource: make_expression_matrix/expression_matrix
    type: File
    label: "A hd5 file containing transcript by sample matrices of TPM and number of reads"

steps:

  - id: fastqc
    in:
      - id: fastq_dir
        source: fastq_dir
      - id: threads
        source: threads
    out:
    - fastqc_dir

    run: steps/fastqc.cwl
    label: "Run fastqc on all fastq files in fastq directory"

  - id: salmon-bulk
    in:
      - id: fastq_dir
        source: fastq_dir
      - id: threads
        source: threads
      - id: organism
        source: organism
    out:
      - output_dir
    run: steps/bulk-salmon.cwl
    when: $(inputs.organism == 'human')
    label: "Salmon quant 1.0.0, with index from GRCh38 transcriptome"

  - id: salmon-bulk-mouse
    in:
      - id: fastq_dir
        source: fastq_dir
      - id: threads
        source: threads
      - id: organism
        source: organism
    out:
      - output_dir
    run: steps/bulk-salmon-mouse.cwl
    when: $(inputs.organism == 'mouse')
    label: "Salmon quant 1.0.0, with index from GRCm39 transcriptome"

  - id: make_expression_matrix
    in:
      - id: quant_dir
        source: [salmon-bulk/output_dir, salmon-bulk-mouse/output_dir]
        pickValue: first_non_null

    out:
      - expression_matrix
    run: steps/expression-matrix.cwl
    label: "Make expression matrix from quant vectors"
