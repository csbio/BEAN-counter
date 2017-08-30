# **BEAN-counter**

BEAN-counter is collection of python (and currently some MATLAB) scripts for processing barcode sequencing data from multiplexed experiments. Originally designed for chemical genomics experiments performed in the Myers/Boone Labs, it is applicable to any experiment in which pools of genetically barcoded cells are grown under different conditions, with the resulting barcode DNA isolated from those cells combined into one 2nd-gen sequencing run via the use of indexed PCR primers.

## Features
- Written in python and MATLAB; will ultimately be pure python
- Provides clear formats for data and information files
- Easy to use for processing raw fastq data to matrices of z-scores
- Modular post-processing steps are highly customizable and can be used in any order
- Efficient numpy implementations

## Installation

### Requirements

BEAN-counter is written in python 2 (and currently MATLAB). It is recommended to download one of the latest versions of python 2.7. The following python libraries are also required (install them using: pip install PACKAGE, or kindly ask your sysadmin):

    numpy
    scipy
    pandas
    matplotlib
    Bio (biopython)
    fastcluster
    networkx (for replicate collapsing only)


### Downloading BEAN-counter

#### Basic

Head on over to https://github.com/csbio/BEAN-counter/releases/ and download the latest release. Extract the compressed folder to a good location from which to run the software (i.e. get it out of your downloads folder!)

#### Advanced

If you know what you are doing and want to keep up-to-date with the latest version, clone the repository (git clone https://github.com/csbio/BEAN-counter.git or windows equivalent).


### Setting up environment variables

Required: **BARSEQ_PATH**
Set the value of this environment variable to the path of the BEAN-counter folder you downloaded and extracted. The scripts from BEAN-counter will look for this variable in your environment, so it must be set!

Optional, but strongly recommended: adding `$BARSEQ_PATH/master_scripts/` to your PATH
Adding the master_scripts folder inside of the BEAN-counter folder to your **PATH** environment variable allows you to execute the top-level scripts in the pipeline by calling them only by their names.

Optional: adding `$BARSEQ_PATH/scripts/` to your **PATH**
You may never need to run the scripts inside the "scripts" folder, but many of them can be run individually. This will make that process easier.

#### How do I set my environment variables?

##### Linux/Mac
The best way to do this is by adding code to the scripts that run every time you open a new shell. If you use the bash shell, then add the following line to either your ~/.bashrc or ~/.bash_profile files:

```
export BARSEQ_PATH=/your/path/to/BEAN-counter/
```

If you use the c shell (csh), then add the following line to your ~/.cshrc file:

```
setenv BARSEQ_PATH /your/path/to/BEAN-counter/
```

To append a directory to your PATH variable, add this line to your ~/.bashrc or ~/.bash_profile (or equivalent for ~/.cshrc):

```
export PATH=$PATH:/your/path/to/BEAN-counter/master_scripts
```

##### Windows

[Tutorial for changing path variable](http://www.computerhope.com/issues/ch000549.htm)

## How to Run
```
process_screen.py [configuration file]
visualize_zscore_matrices.py [configuration file] [formatting file]
```

## Documentation
Details about the software and how to install and run it are here:
https://sites.google.com/site/barseqcounterdocs/
