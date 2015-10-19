#!/usr/bin/env python
# This script takes barseq sequencing data through the following steps:
# 1) Counting barcodes and index tags (individual lanes)
# 2) Calculating chemical-genetic interaction z-scores (individual lanes)
# 3) Removing index tags with strong profile correlations
# 4) Merging the count matrices
# 5) Filtering the dataset based on specified thresholds
# 6) Calculating chemical-genetic interaction z-scores on the entire large dataset

import pandas as pd
import numpy as np
import scipy
import sys, os, gzip
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import itertools as it
import cPickle
import time

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append(os.path.join(barseq_path, 'scripts'))
sys.path.append(os.path.join(barseq_path, 'lib'))

import config_file_parser as cfp
import compressed_file_opener as cfo
import cg_file_tools as cg_file
import cluster_dataset_wrappers as clus_wrap

# Import all of the processing scripts as libraries
import raw_fastq_to_count_matrix
import counts_to_zscores
import mtag_correlations
import merge_count_matrices
import filter_final_count_matrix
import mtag_normalization
import check_mtag_effect

# Function definitions
def get_all_lane_ids(sample_table):
    return np.unique(np.array(sample_table['lane']))

def get_sample_table(config_params):
 
    filename = config_params['sample_table_file']

    # Read everything in as a string, to prevent vexing
    # number interpretation problems! Methods further down
    # can coerce to different types.
    tab = pd.read_table(filename, dtype = 'S')
    return tab

###########################################
#######  Here is the main script  #########
###########################################

# Get the config file, which is the only argument needed for the pipeline
config_file = sys.argv[1]
config_params = cfp.parse(config_file)

# Read in the sample table
sample_table = get_sample_table(config_params)

# Grab all of the ids of the lanes to process
lane_ids = get_all_lane_ids(sample_table)

## Or, if you were silly and ran all lanes but the newest 4, add this in
#lane_ids = ['lane51', 'lane52', 'lane53', 'lane54']

# First, get one strain X condition count matrix per lane
# This only needs to be run once, unless the barcodes
# or index tags change for some reason.
for lane_id in lane_ids:
    print 'generating read count matrix for lane {}...'.format(lane_id)
    start = time.time()
    raw_fastq_to_count_matrix.main(config_file, lane_id)
    end = time.time()
    print 'time to process = {}'.format(time.strftime('%H:%M:%S', time.gmtime(end - start)))
#    
# An initial round of cg interaction scoring is performed
# at the lane level, 
for lane_id in lane_ids:    
    print 'generating z-score matrix for lane {}...'.format(lane_id)
    start = time.time()
    counts_to_zscores.main(config_file, lane_id)
    end = time.time()
    print 'time to process = {}'.format(time.strftime('%H:%M:%S', time.gmtime(end - start)))

# Calculate index tag (condition) correlations on
# the DMSO profiles, for removal in the matrix
# filtering step
print 'computing index tag correlations...'
mtag_correlations.main(config_file)

# Merge all of the count matrices into one big count matrix
print 'merging count matrices...'
merge_count_matrices.main(config_file)

# Filter out all of the strains and conditions with the following issues:
# 1) They were flagged to be excluded a priori
# 2) For conditions: if their associated index tag was too self-correlated
#    in the control conditions
# 3) If the strains or conditions did not meet the count degree thresholds
#    specified in the config file (advanced options)
print 'filtering the final count matrix...'
filter_final_count_matrix.main(config_file)

# Calculate chemical-genetic interaction z-scores on the entire dataset
print 'generating z-score matrix for all lanes'
start = time.time()
counts_to_zscores.main(config_file, 'all_lanes_filtered')
end = time.time()
print 'time to process = {}'.format(time.strftime('%H:%M:%S', time.gmtime(end - start)))

