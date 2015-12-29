barseq_counter
==============

barseq_counter is collection of python (and currently some MATLAB) scripts for processing barcode sequencing data from multiplexed experiments. Originally designed for chemical genomics experiments performed in the Myers/Boone Labs, it is applicable to any experiment in which pools of genetically barcoded cells are grown under different conditions, with the resulting barcode DNA isolated from those cells combined into one 2nd-gen sequencing run via the use of indexed PCR primers.

Features
--------
- Written in python and MATLAB; will ultimately be pure python
- Provides clear formats for data and information files
- Easy to use for processing raw fastq data to matrices of z-scores
- Modular post-processing steps are highly customizable and can be used in any order
- Efficient numpy implementations
