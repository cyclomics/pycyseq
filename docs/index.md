# Introduction

`pycyseq` is a Python library for reading, parsing, and analyzing CySeq sequencing data.

CySeq data is organized according to the [SAM format specification](https://samtools.github.io/hts-specs/SAMv1.pdf) and stored in the standard `.sam` or `.bam` file format.

CySeq files contain read alignments derived from both single-strand DNA (ssDNA) and double-strand DNA (dsDNA). In addition to the alignment information, CySeq reads derived from dsDNA encode strand-specific mutations and overhang type and length information.

`pycyseq` is built on top of `pysam`. It provides a high-level interface for working with CySeq data, making it easy to read and parse alignments, and extract a variety of molecular features from every CySeq read. When lower-level control is required, the underlying `pysam` API remains fully accessible.

!!! warning

    `pycyseq` is still under active delopment, expect code breaking changes.