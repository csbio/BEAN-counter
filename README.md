# **BEAN-counter**

BEAN-counter is a pipeline written in Python for processing barcode sequencing data from multiplexed experiments. Originally designed for chemical genomics experiments performed in the Myers/Boone Labs, it is applicable to any experiment in which pools of genetically barcoded cells are grown under different conditions, with the resulting barcode DNA isolated from those cells combined into one 2nd-gen sequencing run via the use of indexed PCR primers.

## Features
- Written in Python
- Provides clear formats for data and information files
- Easy to use for processing raw fastq data to matrices of z-scores
- Modular post-processing steps are highly customizable and can be used in any order
- Efficient numpy implementations

## Installation

### Requirements

BEAN-counter is written in Python 2. It is recommended to download the latest version of Python 2.7. The following Python libraries are also required:

    numpy (>=1.12.1)
    scipy (>=0.19.0)
    pandas (>=0.20.1)
    matplotlib (>=2.0.2)
    scikit-learn (>=0.18.1)
    biopython (>=1.68)
    fastcluster (>=1.1.20)
    pyyaml (>= 3.11)
    networkx (>=1.11)
    jellyfish (>=0.2.0)
    - Two python libraries exist with the name "jellyfish." BEAN-counter
      requires the installation of the jellyfish library for string distance
      computation, not the library for kmer counting.

Conda (https://conda.io/docs/user-guide/install/index.html) provides a way to
quickly install python and the required libraries. If you install conda, then
the required packages above can be obtained from the Anaconda Cloud
(https://anaconda.org/).

### Downloading BEAN-counter

#### Basic

Download the latest release from https://github.com/csbio/BEAN-counter/releases/.

#### Advanced

If you know what you are doing and want to keep up-to-date with the latest version, clone the repository (git clone https://github.com/csbio/BEAN-counter.git or windows equivalent).

### Setting up environment variables

Required: **BARSEQ_PATH**
Set the value of this environment variable to the path of the BEAN-counter folder you downloaded and extracted. The scripts from BEAN-counter will look for this variable in your environment, so it must be set!

Optional, but strongly recommended: adding `$BARSEQ_PATH/master_scripts/` to your PATH
Adding the master_scripts folder inside of the BEAN-counter folder to your **PATH** environment variable allows you to execute the top-level scripts in the pipeline by calling them only by their names.

## How to Run
```
process_screen.py [configuration file]
visualize_zscore_matrices.py [configuration file] [barcode_table_columns] [sample_table_columns]
```

## Documentation
Details about the software and how to install and run it are here:
https://sites.google.com/site/barseqcounterdocs/
