version 1.0

workflow SalmonQuantification {
  
  # Define inputs
  input {
    Array[File] fastq_dir
    String assay
    Int threads = 1
    Optional[Int] expected_cell_count
    Optional[Boolean] keep_all_barcodes
    Optional[Array[File]] img_dir
    Optional[Array[File]] metadata_dir
    Optional[String] organism
  }

  # Define outputs
  output {
    Directory salmon_output = {Salmon.output_dir, SalmonMouse.output_dir}[0]
    File count_matrix_h5ad = AnnotateCells.annotated_h5ad_file
    Optional[File] raw_count_matrix = AlevinToAnnData.raw_expr_h5ad
    File genome_build_json = AlevinToAnnData.genome_build_json
  }

  # Call AdjustBarcodes task
  scatter (fastq in fastq_dir) {
    output {
      Directory adj_fastq = "adj_fastq"
      Optional[File] metadata_json = "metadata.json"
    }
    command {
      """
      /opt/adjust_barcodes.py {assay} {fastq} {write_lines(adj_fastq)}
      """
    }
    runtime {
      docker "your_container_image"
    }
  }

  # Call TrimReads task
  scatter (fastq in fastq_dir) {
    input {
      Directory adj_fastq = AdjustBarcodes.adj_fastq[fastq]
    }
    output {
      Directory trimmed_fastq = "trimmed"
    }
    command {
      """
      /opt/trim_reads.py {assay} {AdjustBarcodes.adj_fastq[fastq]} {fastq} --threads {threads}
      """
    }
    runtime {
      docker "your_container_image"
    }
  }

  # Call Salmon and SalmonMouse tasks in parallel
  work {
    call Salmon {
      input {
        Array[File] orig_fastq_dirs = fastq_dir
        Directory trimmed_fastq = TrimReads.trimmed_fastq[fastq]
      }
      output {
        Directory output_dir = "salmon_out"
      }
      command {
        """
        /opt/salmon_wrapper.py {assay} {TrimReads.trimmed_fastq[fastq]} {fastq_dir} --threads {threads} 
          ${if defined expected_cell_count then ' --expected-cell-count ' + expected_cell_count else ''}
          ${if defined keep_all_barcodes && keep_all_barcodes then ' --keep-all-barcodes' else ''}
        """
      }
      runtime {
        docker "your_container_image"
      }
    }
    call SalmonMouse {
      input {
        Array[File] orig_fastq_dirs = fastq_dir
        Directory trimmed_fastq = TrimReads.trimmed_fastq[fastq]
      }
      output {
        Directory output_dir = "salmon_mouse_out"
      }
      command {
        """
        /opt/salmon_mouse_wrapper.py {assay} {TrimReads.trimmed_fastq[fastq]} {fastq_dir} --threads {threads} 
          ${if defined expected_cell_count then ' --expected-cell-count ' + expected_cell_count else ''}
          ${if defined keep_all_barcodes && keep_all_barcodes then ' --keep-all-barcodes' else ''}
        """
      }
      runtime {
        docker "your_container_image"
      }
    }
  }

  # Call AlevinToAnnData task with combined outputs from Salmon and SalmonMouse
  task AlevinToAnnData {
    input {
      String assay
      Array[File] alevin_dir = {Salmon.output_dir, SalmonMouse.output_dir}
      Optional[String] organism
    }
    output {
      File expr_h5ad = "expr.h5ad"
      Optional[File] raw_expr_h5ad = "raw_expr.h5ad"
      File genome_build_json = "genome_build.json"
    }
    command {
      """
      /opt/alevin_to_anndata.py {assay} {write_lines(alefin_dir)} --organism {organism}
      """
    }
    runtime {
      docker "your_container_image"
    }
  }

  # Call AnnotateCells task
  task AnnotateCells {
    input {
      Array[File] orig_fastq_dirs
      String assay
      File h5ad_file
      Optional[Array[File]] img_dir
      Optional[Array[File]] metadata_dir
      Optional[File] metadata_json
    }
    output {
      File annotated_h5ad_file = "annotated_expr.h5ad"
    }
    command {
      """
      /opt/annotate_cells.py {assay} {h5ad_file} {write_tsv(select_all(orig_fastq_dirs))} 
        ${if defined img_dir then ' --img_dir ' + write_tsv(select_all(img_dir)) else ''}
        ${if defined metadata_dir then ' --metadata_dir ' + write_tsv(select_all(metadata_dir)) else ''}
        ${if defined metadata_json then ' --metadata_json ' + metadata_json else ''}
      """
    }
    runtime {
      docker "your_container_image"
    }
  }
}