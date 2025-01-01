version 1.0

workflow salmon_quantification {
    input {
        #Array[File] fastqs
        Int threads
        Int mem_gb
        File fastq1
        File fastq2
        String assay
        String? run_id
        String? protocol
        String? species
        String? img_dir
        String? metadata_dir
        String? organism
        Boolean? keep_all_barcodes
        Int? expected_cell_count
    }

    # Step 1: Adjust Barcodes
    call adjust_barcodes {
        input:
            #fastqs = fastqs,
            fastqs = [fastq1, fastq2],
            assay = assay
    }

    # Step 2: Trim Reads
    call trim_reads {
        input:
            #orig_fastqss = fastqs,
            orig_fastqss = [fastq1, fastq2],
            adj_fastqs = adjust_barcodes.adj_fastqs,
            assay = assay
    }

    # Step 3: Salmon Quantification
    call salmon {
        input:
            #orig_fastqss = fastqs,
            orig_fastqss = [fastq1, fastq2],
            trimmed_fastqs = trim_reads.trimmed_fastqs,
            assay = assay,
            threads = threads,
            mem_gb = mem_gb,
            expected_cell_count = expected_cell_count,
            keep_all_barcodes = keep_all_barcodes,
            organism = organism
    }

    # Step 4: Convert Alevin to AnnData
    call alevin_to_anndata {
        input:
            assay = assay,
            alevin_dir = select_first([salmon.output_dir]),
            species = species
    }

    # Step 5: Annotate Cells
    call annotate_cells {
        input:
            #orig_fastqss = fastqs,
            orig_fastqss = [fastq1, fastq2],
            assay = assay,
            h5ad_file = alevin_to_anndata.expr_h5ad,
            img_dir = select_first([img_dir, "default_value"]),
            metadata_dir = select_first([metadata_dir, "default_value"]),
            metadata_json = select_first([adjust_barcodes.metadata_json, "default_value"])
    }

    output {
        File salmon_output = salmon.output_dir
        File count_matrix_h5ad = alevin_to_anndata.expr_h5ad
        File? raw_count_matrix = alevin_to_anndata.raw_expr_h5ad
        File genome_build_json = alevin_to_anndata.genome_build_json
        File? annotated_h5ad = annotate_cells.annotated_h5ad_file
    }
}

# Task Definitions

task adjust_barcodes {
    input {
        Array[File] fastqs
        String assay
    }

    output {
        String adj_fastqs = "adj_fastq"
        File? metadata_json = "metadata.json"
    }

    command {
        # Command for barcode adjustment
        /opt/adjust_barcodes.py ~{assay} ~{sep=" " fastqs}
    }

    runtime {
        container: "hubmap/scrna-barcode-adj:latest"
    }
}

task trim_reads {
    input {
        Array[File] orig_fastqss
        String adj_fastqs
        String assay
    }

    output {
        String trimmed_fastqs = "trimmed"
    }

    command {
        # Command for trimming reads
        /opt/trim_reads.py ~{assay} ~{adj_fastqs} ~{sep=" " orig_fastqss}
    }

    runtime {
        container: "hubmap/scrna-trim-reads:latest"
    }
}

task salmon {
    input {
        Array[File] orig_fastqss
        String trimmed_fastqs
        String assay
        Int threads
        Int mem_gb
        String? expected_cell_count
        String? keep_all_barcodes
        String? organism
    }

    output {
        String output_dir = "salmon_out"
    }

    command {
        # Command for Salmon quantification (human)
        /opt/salmon_wrapper.py ~{assay} ~{trimmed_fastqs} ~{sep=" " orig_fastqss} --threads ~{threads} ~{if defined(expected_cell_count) then "--expected-cell-count " + expected_cell_count else ""} ~{if defined(keep_all_barcodes) then "--keep-all-barcodes " + keep_all_barcodes else ""}
    }

    runtime {
        container: "hubmap/salmon-grch38:latest"
    }
}

task alevin_to_anndata {
    input {
        String assay
        String alevin_dir
        String? species
    }

    output {
        File? raw_expr_h5ad = "raw_expr.h5ad"
        File expr_h5ad = "expr.h5ad"
        File genome_build_json = "genome_build.json"
    }

    command {
        # Command to convert Alevin output to AnnData object
        /opt/alevin_to_anndata.py ~{assay} ~{alevin_dir}
    }

    runtime {
        container: "hubmap/scrna-analysis:latest"
    }
}

task annotate_cells {
    input {
        Array[File] orig_fastqss
        String assay
        File h5ad_file
        String img_dir
        String metadata_dir
        String metadata_json
    }

    output {
        File annotated_h5ad_file = "expr.h5ad"
    }

    command {
        # Command to annotate cells in the AnnData file
        /opt/annotate_cells.py ~{assay} ~{h5ad_file} ~{sep=" " orig_fastqss} ~{if defined(img_dir) then "--img_dir " + img_dir else ""} ~{if defined(metadata_dir) then "--metadata_dir " + metadata_dir else ""} ~{if defined(metadata_json) then "--metadata_json " + metadata_json else ""}
    }

    runtime {
        container: "hubmap/scrna-analysis:latest"
    }
}
