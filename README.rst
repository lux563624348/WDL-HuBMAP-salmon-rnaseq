.. image:: https://travis-ci.com/hubmapconsortium/salmon-rnaseq.svg?branch=master
    :target: https://travis-ci.com/hubmapconsortium/salmon-rnaseq
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

WDL version of HuBMAP scRNA-seq pipeline: Salmon, Scanpy, scVelo
=================================================

Overview
--------

The HuBMAP scRNA-seq pipeline is built on Salmon, Scanpy, and scVelo, and this is 
implemented as a WDL workflow wrapping command-line tools encapsulated in
Docker containers.



Requirements
------------
If you run it in your machine, download code and demo data:
    git clone https://github.com/lux563624348/WDL-HuBMAP-salmon-rnaseq.git

Running the pipeline requires a WDL workflow execution engine and container runtime;
we recommend Docker and the ``dockstore`` reference implementation.

``dockstore`` can be installed by offical tutorial https://dockstore.org/quick-start 

Afterward, clone this repository and enter, and invoke the pipeline as::

    dockstore workflow launch --local-entry ./pipeline.wdl --json test.wdl.json

To avoid docker images error, highly recommend you to download in advance by:: 
    docker pull hubmap/salmon-grch38:latest
    
    >= 28GB RAM is required.   
Salmon quantification step contributes highest memory usage due to inclusion of the entire GRCh38 reference genome as
decoy sequences in the quantification index. See
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8 for more details.

(The ``master`` branch and ``latest`` published Docker images may not always
be in sync; checking out a version like ``v2.0.6`` is *highly* recommended
before running the pipeline, unless building Docker images locally.)

Supported assays:

* ``10x_v2`` (single-cell)
* ``10x_v2_sn`` (single-nucleus)
* ``10x_v3`` (single-cell)
* ``10x_v3_sn`` (single-nucleus)
* ``snareseq``
* ``sciseq``
* ``slideseq``
* ``multiome_10x``
