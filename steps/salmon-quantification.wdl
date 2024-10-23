version 1.0

workflow SalmonQuantification {
    input {
        Array[Directory] fastq_dir
        String assay
        Int threads = 1
        Int? expected_cell_count
        Boolean? keep_all_barcodes
        Directory? img_dir
        Directory? metadata_dir
        String? organism
    }

    call AdjustBarcodes {
        input:
            fastq_dir = fastq_dir,
            assay = assay
    }

    call TrimReads {
        input:
            orig_fastq_dirs = fastq_dir,
            adj_fastq_dir = AdjustBarcodes.adj_fastq_dir,
            assay = assay,
            threads = threads
    }

    call Salmon {
        input:
            orig_fastq_dirs = fastq_dir,
            trimmed_fastq_dir = TrimReads.trimmed_fastq_dir,
            assay = assay,
            threads = threads,
            expected_cell_count = expected_cell_count,
            keep_all_barcodes = keep_all_barcodes,
            organism = organism
    }

    call SalmonMouse {
        input:
            orig_fastq_dirs = fastq_dir,
            trimmed_fastq_dir = TrimReads.trimmed_fastq_dir,
            assay = assay,
            threads = threads,
            expected_cell_count = expected_cell_count,
            keep_all_barcodes = keep_all_barcodes,
            organism = organism
    }

    call AlevinToAnnData {
        input:
            assay = assay,
            alevin_dir = [Salmon.output_dir, SalmonMouse.output_dir],
            organism = organism
    }

    call AnnotateCells {
        input:
            orig_fastq_dirs = fastq_dir,
            assay = assay,
            h5ad_file = AlevinToAnnData.expr_h5ad,
            img_dir = img_dir,
            metadata_dir = metadata_dir,
            metadata_json = AdjustBarcodes.metadata_json
    }

    output {
        Directory salmon_output = Salmon.output_dir
        File count_matrix_h5ad = AnnotateCells.annotated_h5ad_file
        File? raw_count_matrix = AlevinToAnnData.raw_expr_h5ad
        File genome_build_json = AlevinToAnnData.genome_build_json
    }
}

task AdjustBarcodes {
    input {
        Array[Directory] fastq_dir
        String assay
    }

    output {
        Directory adj_fastq_dir = "adj_fastq"
        File? metadata_json = "metadata.json"
    }

    command {
        # Command to adjust barcodes (replace with actual command)
        /opt/adjust_barcodes.py ~{assay} directory ~{join(fastq_dir, " ")}
    }

    runtime {
        container: "your_container_image"
    }
}

task TrimReads {
    input {
        Array[Directory] orig_fastq_dirs
        Directory adj_fastq_dir
        String assay
        Int threads
    }

    output {
        Directory trimmed_fastq_dir = "trimmed"
    }

    command {
        # Command to trim reads (replace with actual command)
        /opt/trim_reads.py ~{assay} ~{adj_fastq_dir} ~{join(orig_fastq_dirs, " ")} --threads ~{threads}
    }

    runtime {
        container: "your_container_image"
    }
}

task Salmon {
    input {
        Array[Directory] orig_fastq_dirs
        Directory trimmed_fastq_dir
        String assay
        Int threads
        Int? expected_cell_count
        Boolean? keep_all_barcodes
        String? organism
    }

    output {
        Directory output_dir = "salmon_out"
    }

    command {
        # Command for Salmon (replace with actual command)
        /opt/salmon_wrapper.py ~{assay} ~{trimmed_fastq_dir} ~{join(orig_fastq_dirs, " ")} --threads ~{threads} ~{if defined(expected_cell_count) then "--expected-cell-count " + expected_cell_count else ""} ~{if defined(keep_all_barcodes) && keep_all_barcodes then "--keep-all-barcodes" else ""}
    }

    runtime {
        container: "your_container_image"
    }
}

task SalmonMouse {
    input {
        Array[Directory] orig_fastq_dirs
        Directory trimmed_fastq_dir
        String assay
        Int threads
        Int? expected_cell_count
        Boolean? keep_all_barcodes
        String? organism
    }

    output {
        Directory output_dir = "salmon_mouse_out"
    }

    command {
        # Command for Salmon Mouse (replace with actual command)
        /opt/salmon_mouse_wrapper.py ~{assay} ~{trimmed_fastq_dir} ~{join(orig_fastq_dirs, " ")} --threads ~{threads} ~{if defined(expected_cell_count) then "--expected-cell-count " + expected_cell_count else ""} ~{if defined(keep_all_barcodes) && keep_all_barcodes then "--keep-all-barcodes" else ""}
    }

    runtime {
        container: "your_container_image"
    }
}

task AlevinToAnnData {
    input {
        String assay
        Array[Directory] alevin_dir
        String? organism
    }

    output {
        File expr_h5ad = "expr.h5ad"
        File? raw_expr_h5ad = "raw_expr.h5ad"
        File genome_build_json = "genome_build.json"
    }

    command {
        # Command to convert Alevin output to AnnData
        /opt/alevin_to_anndata.py ~{assay} ~{join(alevin_dir, " ")} --organism ~{organism}
    }

    runtime {
        container: "your_container_image"
    }
}

task AnnotateCells {
    input {
        Array[Directory] orig_fastq_dirs
        String assay
        File h5ad_file
        Directory? img_dir
        Directory? metadata_dir
        File? metadata_json
    }

    output {
        File annotated_h5ad_file = "annotated_expr.h5ad"
    }

    command {
        # Command to annotate cells
        /opt/annotate_cells.py ~{assay} ~{h5ad_file} ~{join(orig_fastq_dirs, " ")} ~{if defined(img_dir) then "--img_dir " + img_dir else ""} ~{if defined(metadata_dir) then "--metadata_dir " + metadata_dir else ""} ~{if defined(metadata_json) then "--metadata_json " + metadata_json else ""}
    }

    runtime {
        container: "your_container_image"
    }
}
