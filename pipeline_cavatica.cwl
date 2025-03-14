{
    "class": "Workflow",
    "cwlVersion": "v1.2",
    "id": "lux563624348/hubmapconsortium-salmon-rnaseq/hubmap-salmon-rnaseq/3",
    "label": "scRNA-seq pipeline using Salmon and Alevin",
    "$namespaces": {
        "sbg": "https://www.psc.edu/"
    },
    "inputs": [
        {
            "id": "fastq_dir",
            "type": "Directory[]",
            "label": "Directory containing FASTQ files",
            "sbg:fileId": [],
            "loadListing": "deep_listing",
            "sbg:x": 0,
            "sbg:y": 1076.0001220703125
        },
        {
            "id": "img_dir",
            "type": "Directory?",
            "label": "Directory containing TIFF image data (for Visium assay)",
            "sbg:fileId": null,
            "loadListing": "deep_listing",
            "sbg:x": 0,
            "sbg:y": 969.3333740234375
        },
        {
            "id": "metadata_dir",
            "type": "Directory?",
            "label": "Directory containing metadata, including gpr slide file (for Visium assay)",
            "sbg:fileId": null,
            "loadListing": "deep_listing",
            "sbg:x": 0,
            "sbg:y": 755.3333129882812
        },
        {
            "id": "assay",
            "type": "string",
            "label": "scRNA-seq assay",
            "sbg:x": 0,
            "sbg:y": 1290.0001220703125
        },
        {
            "id": "threads",
            "type": "int",
            "label": "Number of threads for Salmon",
            "default": 1,
            "sbg:x": 0,
            "sbg:y": 421.3333435058594
        },
        {
            "id": "expected_cell_count",
            "type": "int?",
            "sbg:x": 0,
            "sbg:y": 1183.0001220703125
        },
        {
            "id": "keep_all_barcodes",
            "type": "boolean?",
            "sbg:x": 0,
            "sbg:y": 862.3333129882812
        },
        {
            "id": "organism",
            "type": "string?",
            "sbg:x": 0,
            "sbg:y": 648.6666259765625
        }
    ],
    "outputs": [
        {
            "id": "salmon_output",
            "outputSource": [
                "salmon/output_dir"
            ],
            "type": "Directory",
            "label": "Full output of `salmon alevin`",
            "sbg:fileId": true,
            "sbg:x": 606.128173828125,
            "sbg:y": 1193.5364990234375
        },
        {
            "id": "count_matrix_h5ad",
            "outputSource": [
                "alevin_to_anndata/expr_h5ad"
            ],
            "type": "File",
            "label": "Unfiltered count matrix from Alevin, converted to H5AD, spliced and unspliced counts",
            "sbg:fileId": true,
            "sbg:x": 540,
            "sbg:y": 1899.184814453125
        },
        {
            "id": "raw_count_matrix",
            "outputSource": [
                "alevin_to_anndata/raw_expr_h5ad"
            ],
            "type": "File?",
            "label": "Unfiltered count matrix from Alevin, converted to H5AD, with intronic counts as separate columns",
            "sbg:fileId": true,
            "sbg:x": 180,
            "sbg:y": 1854.1728515625
        },
        {
            "id": "genome_build_json",
            "outputSource": [
                "alevin_to_anndata/genome_build_json"
            ],
            "type": "File",
            "label": "Genome build information in JSON format",
            "sbg:fileId": true,
            "sbg:x": 360,
            "sbg:y": 1876.6788330078125
        },
        {
            "id": "fastqc_dir",
            "outputSource": [
                "fastqc/fastqc_dir"
            ],
            "type": "Directory[]",
            "label": "Directory of FastQC output files, mirroring input directory structure",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 1390.0001220703125
        },
        {
            "id": "scanpy_qc_results",
            "outputSource": [
                "compute_qc_results/scanpy_qc_results"
            ],
            "type": "File",
            "label": "Quality control metrics from Scanpy",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 534.6666870117188
        },
        {
            "id": "qc_report",
            "outputSource": [
                "compute_qc_results/qc_metrics"
            ],
            "type": "File",
            "label": "Quality control report in JSON format",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 748.6666259765625
        },
        {
            "id": "dispersion_plot",
            "outputSource": [
                "scanpy_analysis/dispersion_plot"
            ],
            "type": "File",
            "label": "Gene expression dispersion plot",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 1496.6668701171875
        },
        {
            "id": "umap_plot",
            "outputSource": [
                "scanpy_analysis/umap_plot"
            ],
            "type": "File",
            "label": "UMAP dimensionality reduction plot",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 0
        },
        {
            "id": "umap_density_plot",
            "outputSource": [
                "scanpy_analysis/umap_density_plot"
            ],
            "type": "File",
            "label": "UMAP dimensionality reduction plot, colored by cell density",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 106.66667175292969
        },
        {
            "id": "spatial_plot",
            "outputSource": [
                "scanpy_analysis/spatial_plot"
            ],
            "type": "File?",
            "label": "Slide-seq bead plot, colored by Leiden cluster",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 428
        },
        {
            "id": "filtered_data_h5ad",
            "outputSource": [
                "scanpy_analysis/filtered_data_h5ad"
            ],
            "type": "File",
            "label": "Full data set of filtered results",
            "doc": "Full data set of filtered results: expression matrix, coordinates in dimensionality-reduced space (PCA and UMAP), cluster assignments via the Leiden algorithm, and marker genes for one cluster vs. rest",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 1283.33349609375
        },
        {
            "id": "marker_gene_plot_t_test",
            "outputSource": [
                "scanpy_analysis/marker_gene_plot_t_test"
            ],
            "type": "File",
            "label": "Cluster marker genes, t-test",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 962.6666870117188
        },
        {
            "id": "marker_gene_plot_logreg",
            "outputSource": [
                "scanpy_analysis/marker_gene_plot_logreg"
            ],
            "type": "File",
            "label": "Cluster marker genes, logreg method",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 1069.33349609375
        },
        {
            "id": "scvelo_annotated_h5ad",
            "outputSource": [
                "scvelo_analysis/annotated_h5ad_file"
            ],
            "type": "File?",
            "label": "scVelo-annotated h5ad file, including cell RNA velocity",
            "sbg:fileId": true,
            "sbg:x": 411.95135498046875,
            "sbg:y": 530.3333740234375
        },
        {
            "id": "scvelo_embedding_grid_plot",
            "outputSource": [
                "scvelo_analysis/embedding_grid_plot"
            ],
            "type": "File?",
            "label": "scVelo velocity embedding grid plot",
            "sbg:fileId": true,
            "sbg:x": 411.95135498046875,
            "sbg:y": 423.6667175292969
        },
        {
            "id": "squidpy_annotated_h5ad",
            "outputSource": [
                "squidpy_analysis/squidpy_annotated_h5ad"
            ],
            "type": "File?",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 321
        },
        {
            "id": "neighborhood_enrichment_plot",
            "outputSource": [
                "squidpy_analysis/neighborhood_enrichment_plot"
            ],
            "type": "File?",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 855.6666259765625
        },
        {
            "id": "co_occurrence_plot",
            "outputSource": [
                "squidpy_analysis/co_occurrence_plot"
            ],
            "type": "File?",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 1603.6668701171875
        },
        {
            "id": "interaction_matrix_plot",
            "outputSource": [
                "squidpy_analysis/interaction_matrix_plot"
            ],
            "type": "File?",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 1176.33349609375
        },
        {
            "id": "centrality_scores_plot",
            "outputSource": [
                "squidpy_analysis/centrality_scores_plot"
            ],
            "type": "File?",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 1711.0001220703125
        },
        {
            "id": "ripley_plot",
            "outputSource": [
                "squidpy_analysis/ripley_plot"
            ],
            "type": "File?",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 641.6666259765625
        },
        {
            "id": "squidpy_spatial_plot",
            "outputSource": [
                "squidpy_analysis/spatial_plot"
            ],
            "type": "File?",
            "sbg:fileId": true,
            "sbg:x": 820.28466796875,
            "sbg:y": 213.6666717529297
        }
    ],
    "steps": [
        {
            "id": "fastqc",
            "in": [
                {
                    "id": "fastq_dir",
                    "source": "fastq_dir",
                    "loadListing": "deep_listing"
                },
                {
                    "id": "threads",
                    "source": "threads"
                }
            ],
            "out": [
                {
                    "id": "fastqc_dir"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "$namespaces": {
                    "sbg": "https://github.com/hubmapconsortium"
                },
                "baseCommand": [
                    "/opt/fastqc_wrapper.py"
                ],
                "inputs": [
                    {
                        "sbg:fileId": [],
                        "loadListing": "deep_listing",
                        "id": "fastq_dir",
                        "type": "Directory",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        },
                        "doc": "Directory containing fastq files to be evaluated"
                    },
                    {
                        "id": "threads",
                        "type": "int",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 2
                        },
                        "doc": "The number of threads to use for fastqc"
                    }
                ],
                "outputs": [
                    {
                        "id": "fastqc_dir",
                        "doc": "Individual graph files and additional data files containing the raw data from which plots were drawn.",
                        "type": "Directory",
                        "outputBinding": {
                            "glob": "fastqc_output",
                            "loadListing": "deep_listing"
                        },
                        "sbg:fileId": true
                    }
                ],
                "label": "Runs fastQC on each fastq file in fastq directory",
                "requirements": [
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/scrna-analysis:latest"
                    }
                ],
                "sbg:appVersion": [
                    "v1.2"
                ],
                "sbg:revision": 0,
                "sbg:revisionNotes": null,
                "sbg:modifiedOn": 1741739028,
                "sbg:modifiedBy": "lux563624348",
                "sbg:createdOn": 1741739028,
                "sbg:createdBy": "lux563624348",
                "sbg:sbgMaintained": false,
                "sbg:contributors": [
                    "lux563624348"
                ],
                "sbg:latestRevision": 0,
                "sbg:revisionsInfo": [
                    {
                        "sbg:revision": 0,
                        "sbg:modifiedBy": "lux563624348",
                        "sbg:modifiedOn": 1741739028,
                        "sbg:revisionNotes": null
                    }
                ],
                "sbg:content_hash": "some_hash_value_replace_with_real_hash",
                "sbg:workflowLanguage": "CWL"
            },
            "label": "Run fastqc on all fastq files in fastq directory",
            "scatter": [
                "fastq_dir"
            ],
            "scatterMethod": "dotproduct",
            "sbg:x": 411.95135498046875,
            "sbg:y": 1011.3333740234375
        },
        {
            "id": "adjust_barcodes",
            "in": [
                {
                    "id": "assay",
                    "source": "assay"
                },
                {
                    "id": "fastq_dir",
                    "source": [
                        "fastqc/fastqc_dir"
                    ],
                    "loadListing": "deep_listing"
                }
            ],
            "out": [
                {
                    "id": "adj_fastq_dir"
                },
                {
                    "id": "metadata_json"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "baseCommand": [
                    "/opt/adjust_barcodes.py"
                ],
                "inputs": [
                    {
                        "id": "assay",
                        "type": "string",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 0
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "fastq_dir",
                        "type": "Directory[]",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "adj_fastq_dir",
                        "type": "Directory",
                        "outputBinding": {
                            "glob": "adj_fastq",
                            "loadListing": "deep_listing"
                        }
                    },
                    {
                        "id": "metadata_json",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "metadata.json"
                        }
                    }
                ],
                "label": "Assay-specific adjustment of cell barcodes",
                "requirements": [
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/scrna-barcode-adj:latest"
                    }
                ]
            },
            "label": "Assay-specific adjustment of cell barcodes",
            "sbg:x": 411.95135498046875,
            "sbg:y": 1578.0001220703125
        },
        {
            "id": "trim_reads",
            "in": [
                {
                    "id": "assay",
                    "source": "assay"
                },
                {
                    "id": "adj_fastqc_dir",
                    "source": "adjust_barcodes/adj_fastq_dir",
                    "loadListing": "deep_listing"
                },
                {
                    "id": "orig_fastq_dirs",
                    "source": [
                        "fastq_dir"
                    ],
                    "loadListing": "deep_listing"
                },
                {
                    "id": "threads",
                    "source": "threads"
                }
            ],
            "out": [
                {
                    "id": "trimmed_fastq_dir"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "baseCommand": [
                    "/opt/trim_reads.py"
                ],
                "inputs": [
                    {
                        "id": "assay",
                        "type": "string",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 0
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "adj_fastqc_dir",
                        "type": "Directory",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "orig_fastq_dirs",
                        "type": "Directory[]",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 2
                        }
                    },
                    {
                        "id": "threads",
                        "type": "int",
                        "inputBinding": {
                            "prefix": "--threads",
                            "shellQuote": true,
                            "position": 3
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "trimmed_fastq_dir",
                        "type": "Directory",
                        "outputBinding": {
                            "glob": "trimmed",
                            "loadListing": "deep_listing"
                        }
                    }
                ],
                "label": "Trim FASTQ files",
                "requirements": [
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/scrna-trim-reads:latest"
                    }
                ]
            },
            "label": "Trim FASTQ files",
            "sbg:x": 411.95135498046875,
            "sbg:y": 105.33338928222656
        },
        {
            "id": "salmon",
            "in": [
                {
                    "id": "assay",
                    "source": "assay"
                },
                {
                    "id": "trimmed_fastq_dir",
                    "source": "trim_reads/trimmed_fastq_dir",
                    "loadListing": "deep_listing"
                },
                {
                    "id": "orig_fastq_dirs",
                    "source": [
                        "fastq_dir"
                    ],
                    "loadListing": "deep_listing"
                },
                {
                    "id": "threads",
                    "source": "threads"
                },
                {
                    "id": "expected_cell_count",
                    "source": "expected_cell_count"
                },
                {
                    "id": "keep_all_barcodes",
                    "source": "keep_all_barcodes"
                }
            ],
            "out": [
                {
                    "id": "output_dir"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "baseCommand": [
                    "/opt/salmon_wrapper.py"
                ],
                "inputs": [
                    {
                        "id": "assay",
                        "type": "string",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 0
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "trimmed_fastq_dir",
                        "type": "Directory",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "orig_fastq_dirs",
                        "type": "Directory[]",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 2
                        }
                    },
                    {
                        "id": "threads",
                        "type": "int",
                        "inputBinding": {
                            "prefix": "--threads",
                            "shellQuote": true,
                            "position": 3
                        }
                    },
                    {
                        "id": "expected_cell_count",
                        "type": "int?",
                        "inputBinding": {
                            "prefix": "--expected-cell-count",
                            "shellQuote": true,
                            "position": 4
                        }
                    },
                    {
                        "id": "keep_all_barcodes",
                        "type": "boolean?",
                        "inputBinding": {
                            "prefix": "--keep-all-barcodes",
                            "shellQuote": true,
                            "position": 5
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "output_dir",
                        "type": "Directory",
                        "outputBinding": {
                            "glob": "salmon_out",
                            "loadListing": "deep_listing"
                        }
                    }
                ],
                "label": "Run Salmon Alevin tool on FASTQ input",
                "requirements": [
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "ResourceRequirement",
                        "ramMin": 28672
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/salmon-grch38:latest"
                    }
                ]
            },
            "label": "Run Salmon Alevin tool on FASTQ input",
            "sbg:x": 411.95135498046875,
            "sbg:y": 862.6666870117188
        },
        {
            "id": "alevin_to_anndata",
            "in": [
                {
                    "id": "assay",
                    "source": "assay"
                },
                {
                    "id": "alevin_dir",
                    "source": "salmon/output_dir",
                    "loadListing": "deep_listing"
                },
                {
                    "id": "organism",
                    "source": "organism"
                }
            ],
            "out": [
                {
                    "id": "raw_expr_h5ad"
                },
                {
                    "id": "expr_h5ad"
                },
                {
                    "id": "genome_build_json"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "baseCommand": [
                    "/opt/alevin_to_anndata.py"
                ],
                "inputs": [
                    {
                        "id": "assay",
                        "type": "string",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 0
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "alevin_dir",
                        "type": "Directory",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        }
                    },
                    {
                        "id": "organism",
                        "type": "string?",
                        "inputBinding": {
                            "prefix": "--organism",
                            "shellQuote": true,
                            "position": 2
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "raw_expr_h5ad",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "raw_expr.h5ad"
                        }
                    },
                    {
                        "id": "expr_h5ad",
                        "type": "File",
                        "outputBinding": {
                            "glob": "expr.h5ad"
                        }
                    },
                    {
                        "id": "genome_build_json",
                        "type": "File",
                        "outputBinding": {
                            "glob": "genome_build.json"
                        }
                    }
                ],
                "label": "Convert Alevin sparse output to anndata.AnnData object, save as h5ad",
                "requirements": [
                    {
                        "class": "MultipleInputFeatureRequirement"
                    },
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/scrna-analysis:latest"
                    }
                ]
            },
            "label": "Convert Alevin sparse output to anndata.AnnData object, save as h5ad",
            "sbg:x": 411.95135498046875,
            "sbg:y": 1450.33349609375
        },
        {
            "id": "annotate_cells",
            "in": [
                {
                    "id": "assay",
                    "source": "assay"
                },
                {
                    "id": "h5ad_file",
                    "source": "alevin_to_anndata/expr_h5ad"
                },
                {
                    "id": "orig_fastq_dirs",
                    "source": [
                        "fastq_dir"
                    ],
                    "loadListing": "deep_listing"
                },
                {
                    "id": "img_dir",
                    "source": "img_dir",
                    "loadListing": "deep_listing"
                },
                {
                    "id": "metadata_dir",
                    "source": "metadata_dir",
                    "loadListing": "deep_listing"
                },
                {
                    "id": "metadata_json",
                    "source": "alevin_to_anndata/genome_build_json"
                }
            ],
            "out": [
                {
                    "id": "annotated_h5ad_file"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "baseCommand": [
                    "/opt/annotate_cells.py"
                ],
                "inputs": [
                    {
                        "id": "assay",
                        "type": "string",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 0
                        }
                    },
                    {
                        "id": "h5ad_file",
                        "type": "File",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "orig_fastq_dirs",
                        "type": "Directory[]",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 2
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "img_dir",
                        "type": "Directory?",
                        "inputBinding": {
                            "prefix": "--img_dir",
                            "shellQuote": true,
                            "position": 3
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "metadata_dir",
                        "type": "Directory?",
                        "inputBinding": {
                            "prefix": "--metadata_dir",
                            "shellQuote": true,
                            "position": 4
                        }
                    },
                    {
                        "id": "metadata_json",
                        "type": "File?",
                        "inputBinding": {
                            "prefix": "--metadata_json",
                            "shellQuote": true,
                            "position": 5
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "annotated_h5ad_file",
                        "type": "File",
                        "outputBinding": {
                            "glob": "expr.h5ad"
                        }
                    }
                ],
                "label": "Assay-specific annotation of cell barcodes after quantification",
                "requirements": [
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/scrna-analysis:latest"
                    }
                ]
            },
            "label": "Assay-specific annotation of cell barcodes after quantification",
            "sbg:x": 411.95135498046875,
            "sbg:y": 1301.666748046875
        },
        {
            "id": "scanpy_analysis",
            "in": [
                {
                    "id": "assay",
                    "source": "assay"
                },
                {
                    "id": "h5ad_file",
                    "source": "alevin_to_anndata/expr_h5ad"
                }
            ],
            "out": [
                {
                    "id": "filtered_data_h5ad"
                },
                {
                    "id": "dispersion_plot"
                },
                {
                    "id": "umap_plot"
                },
                {
                    "id": "spatial_plot"
                },
                {
                    "id": "umap_density_plot"
                },
                {
                    "id": "marker_gene_plot_t_test"
                },
                {
                    "id": "marker_gene_plot_logreg"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "baseCommand": [
                    "/opt/scanpy_entry_point.py"
                ],
                "inputs": [
                    {
                        "id": "assay",
                        "type": "string",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 0
                        }
                    },
                    {
                        "id": "h5ad_file",
                        "type": "File",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "filtered_data_h5ad",
                        "type": "File",
                        "outputBinding": {
                            "glob": "secondary_analysis.h5ad"
                        }
                    },
                    {
                        "id": "dispersion_plot",
                        "type": "File",
                        "outputBinding": {
                            "glob": "dispersion_plot.pdf"
                        }
                    },
                    {
                        "id": "umap_plot",
                        "type": "File",
                        "outputBinding": {
                            "glob": "umap_by_leiden_cluster.pdf"
                        }
                    },
                    {
                        "id": "spatial_plot",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "spatial_pos_by_leiden_cluster.pdf"
                        }
                    },
                    {
                        "id": "umap_density_plot",
                        "type": "File",
                        "outputBinding": {
                            "glob": "umap_embedding_density.pdf"
                        }
                    },
                    {
                        "id": "marker_gene_plot_t_test",
                        "type": "File",
                        "outputBinding": {
                            "glob": "marker_genes_by_cluster_t_test.pdf"
                        }
                    },
                    {
                        "id": "marker_gene_plot_logreg",
                        "type": "File",
                        "outputBinding": {
                            "glob": "marker_genes_by_cluster_logreg.pdf"
                        }
                    }
                ],
                "label": "Dimensionality reduction and clustering",
                "requirements": [
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/scrna-analysis:latest"
                    }
                ]
            },
            "label": "Dimensionality reduction and clustering",
            "sbg:x": 411.95135498046875,
            "sbg:y": 679
        },
        {
            "id": "scvelo_analysis",
            "in": [
                {
                    "id": "spliced_h5ad_file",
                    "source": "scanpy_analysis/filtered_data_h5ad"
                },
                {
                    "id": "assay_name",
                    "source": "assay"
                }
            ],
            "out": [
                {
                    "id": "annotated_h5ad_file"
                },
                {
                    "id": "embedding_grid_plot"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "baseCommand": [
                    "/opt/scvelo_analysis.py"
                ],
                "inputs": [
                    {
                        "id": "spliced_h5ad_file",
                        "type": "File",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        }
                    },
                    {
                        "id": "assay_name",
                        "type": "string",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 2
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "annotated_h5ad_file",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "scvelo_annotated.h5ad"
                        }
                    },
                    {
                        "id": "embedding_grid_plot",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "scvelo_embedding_grid.pdf"
                        }
                    }
                ],
                "label": "RNA velocity analysis via scVelo",
                "requirements": [
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/scrna-analysis:latest"
                    },
                    {
                        "class": "EnvVarRequirement",
                        "envDef": [
                            {
                                "envName": "TMPDIR",
                                "envValue": "/tmp"
                            }
                        ]
                    }
                ]
            },
            "label": "RNA velocity analysis via scVelo",
            "sbg:x": 0,
            "sbg:y": 535
        },
        {
            "id": "squidpy_analysis",
            "in": [
                {
                    "id": "assay",
                    "source": "assay"
                },
                {
                    "id": "h5ad_file",
                    "source": "alevin_to_anndata/expr_h5ad"
                },
                {
                    "id": "img_dir",
                    "source": "img_dir",
                    "loadListing": "deep_listing"
                }
            ],
            "out": [
                {
                    "id": "squidpy_annotated_h5ad"
                },
                {
                    "id": "neighborhood_enrichment_plot"
                },
                {
                    "id": "co_occurrence_plot"
                },
                {
                    "id": "spatial_plot"
                },
                {
                    "id": "interaction_matrix_plot"
                },
                {
                    "id": "centrality_scores_plot"
                },
                {
                    "id": "ripley_plot"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "baseCommand": [
                    "/opt/squidpy_entry_point.py"
                ],
                "inputs": [
                    {
                        "id": "assay",
                        "type": "string",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 0
                        }
                    },
                    {
                        "id": "h5ad_file",
                        "type": "File",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "img_dir",
                        "type": "Directory?",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 2
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "squidpy_annotated_h5ad",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "squidpy_annotated.h5ad"
                        }
                    },
                    {
                        "id": "neighborhood_enrichment_plot",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "neighborhood_enrichment.pdf"
                        }
                    },
                    {
                        "id": "co_occurrence_plot",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "co_occurrence.pdf"
                        }
                    },
                    {
                        "id": "spatial_plot",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "spatial_scatter.pdf"
                        }
                    },
                    {
                        "id": "interaction_matrix_plot",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "interaction_matrix.pdf"
                        }
                    },
                    {
                        "id": "centrality_scores_plot",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "centrality_scores.pdf"
                        }
                    },
                    {
                        "id": "ripley_plot",
                        "type": "File?",
                        "outputBinding": {
                            "glob": "ripley.pdf"
                        }
                    }
                ],
                "label": "Dimensionality reduction and clustering",
                "requirements": [
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/squidpy-analysis:latest"
                    }
                ]
            },
            "label": "Dimensionality reduction and clustering",
            "sbg:x": 411.95135498046875,
            "sbg:y": 275.00006103515625
        },
        {
            "id": "compute_qc_results",
            "in": [
                {
                    "id": "assay",
                    "source": "assay"
                },
                {
                    "id": "primary_matrix_path",
                    "source": "alevin_to_anndata/expr_h5ad"
                },
                {
                    "id": "secondary_matrix_path",
                    "source": "scanpy_analysis/filtered_data_h5ad"
                },
                {
                    "id": "salmon_dir",
                    "source": "salmon/output_dir",
                    "loadListing": "deep_listing"
                }
            ],
            "out": [
                {
                    "id": "scanpy_qc_results"
                },
                {
                    "id": "qc_metrics"
                }
            ],
            "run": {
                "class": "CommandLineTool",
                "cwlVersion": "v1.2",
                "baseCommand": [
                    "/opt/compute_qc_metrics.py"
                ],
                "inputs": [
                    {
                        "id": "assay",
                        "type": "string",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 0
                        }
                    },
                    {
                        "id": "primary_matrix_path",
                        "type": "File",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 1
                        }
                    },
                    {
                        "id": "secondary_matrix_path",
                        "type": "File",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 2
                        }
                    },
                    {
                        "loadListing": "deep_listing",
                        "id": "salmon_dir",
                        "type": "Directory",
                        "inputBinding": {
                            "shellQuote": true,
                            "position": 3
                        }
                    }
                ],
                "outputs": [
                    {
                        "id": "scanpy_qc_results",
                        "type": "File",
                        "outputBinding": {
                            "glob": "qc_results.hdf5"
                        }
                    },
                    {
                        "id": "qc_metrics",
                        "type": "File",
                        "outputBinding": {
                            "glob": "qc_results.json"
                        }
                    }
                ],
                "label": "Compute QC metrics",
                "requirements": [
                    {
                        "class": "LoadListingRequirement"
                    },
                    {
                        "class": "DockerRequirement",
                        "dockerPull": "hubmap/scrna-analysis:latest"
                    }
                ]
            },
            "label": "Compute QC metrics",
            "sbg:x": 411.95135498046875,
            "sbg:y": 1146.0001220703125
        }
    ],
    "requirements": [
        {
            "class": "LoadListingRequirement"
        },
        {
            "class": "ScatterFeatureRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "StepInputExpressionRequirement"
        }
    ],
    "sbg:projectName": "hubmapconsortium_salmon_rnaseq",
    "sbg:revisionsInfo": [
        {
            "sbg:revision": 0,
            "sbg:modifiedBy": "lux563624348",
            "sbg:modifiedOn": 1741739028,
            "sbg:revisionNotes": null
        },
        {
            "sbg:revision": 1,
            "sbg:modifiedBy": "lux563624348",
            "sbg:modifiedOn": 1741740527,
            "sbg:revisionNotes": "pipeline.cwl"
        },
        {
            "sbg:revision": 2,
            "sbg:modifiedBy": "lux563624348",
            "sbg:modifiedOn": 1741832127,
            "sbg:revisionNotes": "2"
        },
        {
            "sbg:revision": 3,
            "sbg:modifiedBy": "lux563624348",
            "sbg:modifiedOn": 1741832143,
            "sbg:revisionNotes": ""
        }
    ],
    "sbg:image_url": "https://cavatica.sbgenomics.com/ns/brood/images/lux563624348/hubmapconsortium-salmon-rnaseq/hubmap-salmon-rnaseq/3.png",
    "sbg:appVersion": [
        "v1.2"
    ],
    "sbg:id": "lux563624348/hubmapconsortium-salmon-rnaseq/hubmap-salmon-rnaseq/3",
    "sbg:revision": 3,
    "sbg:revisionNotes": "",
    "sbg:modifiedOn": 1741832143,
    "sbg:modifiedBy": "lux563624348",
    "sbg:createdOn": 1741739028,
    "sbg:createdBy": "lux563624348",
    "sbg:project": "lux563624348/hubmapconsortium-salmon-rnaseq",
    "sbg:sbgMaintained": false,
    "sbg:validationErrors": [
        "Required input is not set: #fastqc_dir",
        "Required input is not set: #fastqc_dir",
        "Required input is not set: #fastqc_dir",
        "Required input is not set: #fastqc_dir",
        "Required input is not set: #alevin_dir",
        "Required input is not set: #h5ad_file",
        "Required input is not set: #fastqc_dir",
        "Required input is not set: #h5ad_file",
        "Required input is not set: #spliced_h5ad_file",
        "Required input is not set: #assay_name",
        "Required input is not set: #h5ad_file",
        "Required input is not set: #primary_matrix_path",
        "Required input is not set: #secondary_matrix_path",
        "Required input is not set: #salmon_dir"
    ],
    "sbg:contributors": [
        "lux563624348"
    ],
    "sbg:latestRevision": 3,
    "sbg:publisher": "sbg",
    "sbg:content_hash": "a101f72239e5dea893b33c9fac4fa78218e7072673d2fe4afafccec73b255c674",
    "sbg:workflowLanguage": "CWL"
}