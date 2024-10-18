## Use double '#' for workflow-level comments
## This workflow implements a one-task workflow

# write the WDL version number 'version 1.0' -- 1
# possible to write 'WDL developent' as a version number as well

version 1.0

workflow scRNA_seq_pipeline {
  input {
    Array[Directory] fastq_dir
    Directory? img_dir
    Directory? metadata_dir
    String assay
    Int threads = 1
    Int? expected_cell_count
    Boolean? keep_all_barcodes
    String? organism
  }

  call salmon_quantification {
    input:
      fastq_dir = fastq_dir,
      img_dir = img_dir,
      metadata_dir = metadata_dir,
      assay = assay,
      threads = threads,
      expected_cell_count = expected_cell_count,
      keep_all_barcodes = keep_all_barcodes,
      organism = organism
  }

  call fastqc scatter=fq in {
    fastq_dir = fastq_dir,
    threads = threads
  }

  call scanpy_analysis {
    input:
      assay = assay,
      h5ad_file = salmon_quantification.count_matrix_h5ad
  }

  call scvelo_analysis {
    input:
      spliced_h5ad_file = salmon_quantification.count_matrix_h5ad,
      assay_name = assay
  }

  call squidpy_analysis {
    input:
      assay = assay,
      h5ad_file = scanpy_analysis.filtered_data_h5ad,
      img_dir = img_dir
  }

  call compute_qc_results {
    input:
      assay = assay,
      primary_matrix_path = salmon_quantification.count_matrix_h5ad,
      secondary_matrix_path = scanpy_analysis.filtered_data_h5ad,
      salmon_dir = salmon_quantification.salmon_output
  }

  output {
    Directory salmon_output = salmon_quantification.salmon_output
    File count_matrix_h5ad = salmon_quantification.count_matrix_h5ad
    File? raw_count_matrix = salmon_quantification.raw_count_matrix
    File genome_build_json = salmon_quantification.genome_build_json
    Array[Directory] fastqc_dir = fastqc.fastqc_dir
    File scanpy_qc_results = compute_qc_results.scanpy_qc_results
    File qc_report = compute_qc_results.qc_metrics
    File dispersion_plot = scanpy_analysis.dispersion_plot
    File umap_plot = scanpy_analysis.umap_plot
    File umap_density_plot = scanpy_analysis.umap_density_plot
    File? spatial_plot = scanpy_analysis.spatial_plot
    File filtered_data_h5ad = scanpy_analysis.filtered_data_h5ad
    File marker_gene_plot_t_test = scanpy_analysis.marker_gene_plot_t_test
    File marker_gene_plot_logreg = scanpy_analysis.marker_gene_plot_logreg
    File? scvelo_annotated_h5ad = scvelo_analysis.annotated_h5ad_file
    File? scvelo_embedding_grid_plot = scvelo_analysis.embedding_grid_plot
    File? squidpy_annotated_h5ad = squidpy_analysis.squidpy_annotated_h5ad
    File? neighborhood_enrichment_plot = squidpy_analysis.neighborhood_enrichment_plot
    File? co_occurrence_plot = squidpy_analysis.co_occurrence_plot
    File? interaction_matrix_plot = squidpy_analysis.interaction_matrix_plot
    File? centrality_scores_plot = squidpy_analysis.centrality_scores_plot
    File? ripley_plot = squidpy_analysis.ripley_plot
    File? squidpy_spatial_plot = squidpy_analysis.spatial_plot
  }
}

task salmon_quantification {
  input {
    Array[Directory] fastq_dir
    Directory? img_dir
    Directory? metadata_dir
    String assay
    Int threads
    Int? expected_cell_count
    Boolean? keep_all_barcodes
    String? organism
  }

  command {
    # Command to run Salmon quantification
  }

  output {
    Directory salmon_output = "path/to/salmon_output"
    File count_matrix_h5ad = "path/to/count_matrix.h5ad"
    File? raw_count_matrix = "path/to/raw_count_matrix"
    File genome_build_json = "path/to/genome_build.json"
  }
}

task fastqc {
  input {
    Array[Directory] fastq_dir
    Int threads
  }

  command {
    # Command to run FastQC
  }

  output {
    Array[Directory] fastqc_dir = "path/to/fastqc_output"
  }
}

task scanpy_analysis {
  input {
    String assay
    File h5ad_file
  }

  command {
    # Command for ScanPy analysis
  }

  output {
    File filtered_data_h5ad = "path/to/filtered_data.h5ad"
    File umap_plot = "path/to/umap_plot"
    File marker_gene_plot_t_test = "path/to/marker_gene_plot_t_test"
    File marker_gene_plot_logreg = "path/to/marker_gene_plot_logreg"
    File dispersion_plot = "path/to/dispersion_plot"
    File umap_density_plot = "path/to/umap_density_plot"
    File? spatial_plot = "path/to/spatial_plot"
  }
}

task scvelo_analysis {
  input {
    File spliced_h5ad_file
    String assay_name
  }

  command {
    # Command for scVelo analysis
  }

  output {
    File? annotated_h5ad_file = "path/to/annotated_h5ad_file"
    File? embedding_grid_plot = "path/to/embedding_grid_plot"
  }
}

task squidpy_analysis {
  input {
    String assay
    File h5ad_file
    Directory? img_dir
  }

  command {
    # Command for SquidPy analysis
  }

  output {
    File? squidpy_annotated_h5ad = "path/to/squidpy_annotated_h5ad"
    File? neighborhood_enrichment_plot = "path/to/neighborhood_enrichment_plot"
    File? co_occurrence_plot = "path/to/co_occurrence_plot"
    File? interaction_matrix_plot = "path/to/interaction_matrix_plot"
    File? centrality_scores_plot = "path/to/centrality_scores_plot"
    File? ripley_plot = "path/to/ripley_plot"
    File? spatial_plot = "path/to/spatial_plot"
  }
}

task compute_qc_results {
  input {
    String assay
    File primary_matrix_path
    File secondary_matrix_path
    Directory salmon_dir
  }

  command {
    # Command to compute QC results
  }

  output {
    File scanpy_qc_results = "path/to/scanpy_qc_results"
    File qc_metrics = "path/to/qc_metrics"
  }
}
