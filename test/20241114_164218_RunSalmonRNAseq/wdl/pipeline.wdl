version 1.0
## Use double '#' for workflow-level comments
## This workflow implements a one-task workflow

import "steps/fastqc.wdl" as FastQC
import "steps/salmon-quantification.wdl" as SalmonQuantification

workflow RunSalmonRNAseq {
    meta {
        description: "The SalmonRNAseq pipeline processes"
        author: "Xiang Li"
    }

	input {
        File fastq1
        File fastq2
        Int threads = 1
        String assay = "scrnaseq"
        String? organism
        Int? expected_cell_count
        File? limits
        Boolean? keep_all_barcodes
	}

    parameter_meta {
        fastq1: "fastq1"
        fastq2: "fastq2"
        threads: "# cpu for compute"
        expected_cell_count: "expected_cell_count"
        keep_all_barcodes: "keep_all_barcodes"
    }

    scatter (fastq in [fastq1, fastq2]) {
        call FastQC.fastqc as fastqc {
            input:
                fastqs = [fastq],
                threads = threads,
                limits = limits
        }
    }

    scatter (fastq in [fastq1, fastq2]) {
        call SalmonQuantification.salmon_quantification as SalmonQuantificationCall {
            input:
                fastqs = [fastq],
                threads = threads,
                assay = assay,
                organism = organism,
                expected_cell_count = expected_cell_count,
                keep_all_barcodes = keep_all_barcodes
        }
    }

	output {
		Array[File] reports = flatten(fastqc.reports)
        Array[File] salmon_output = SalmonQuantificationCall.salmon_output
        Array[File] count_matrix_h5ad = SalmonQuantificationCall.count_matrix_h5ad
        Array[File] genome_build_json = SalmonQuantificationCall.genome_build_json
        Array[File?] raw_count_matrix = SalmonQuantificationCall.raw_count_matrix
	}

}
