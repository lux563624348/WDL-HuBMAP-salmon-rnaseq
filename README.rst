.. image:: https://travis-ci.com/hubmapconsortium/salmon-rnaseq.svg?branch=master
    :target: https://travis-ci.com/hubmapconsortium/salmon-rnaseq
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black

HuBMAP scRNA-seq pipeline: Salmon, Scanpy, scVelo
=================================================

Overview
--------

The HuBMAP scRNA-seq pipeline is built on Salmon, Scanpy, and scVelo, and this is 
implemented as a WDL workflow wrapping command-line tools encapsulated in
Docker containers.

Requirements
------------

Running the pipeline requires a WDL workflow execution engine and container
runtime; we recommend Docker and the ``miniwdl`` reference implementation.
``miniwdl`` is written in Python and can be installed into a sufficiently
recent Python environment with ``pip3 install miniwdl``. Afterward, clone this
repository, check out a tag, and invoke the pipeline as::

  miniwdl pipeline.wdl -i test.wdl.json

At least 28GB memory is required for the Salmon quantification step; this
memory usage is due to inclusion of the entire GRCh38 reference genome as
decoy sequences in the quantification index. See
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02151-8
for more details.

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
