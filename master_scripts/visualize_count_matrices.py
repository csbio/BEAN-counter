#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.5.0'

# This script converts all dumped count matrices into CDTs,
# which are visualized using Java Treeview. I anticipate
# future versions will be smart and avoid re-clustering
# datasets if the data themselves have not changed.

import pandas as pd
import numpy as np
import scipy
import sys, os, gzip
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import itertools as it
import cPickle

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'scripts'))
sys.path.append(os.path.join(barseq_path, 'lib'))

import compressed_file_opener as cfo
import cg_file_tools as cg_file
import cluster_dataset_wrappers as clus_wrap
from cg_common_functions import get_sample_table, parse_yaml
from version_printing import update_version_file

import argparse

# Parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-str_cols', '--strain_columns', help = 'the columns from the barcode table to be included in the visualization')
parser.add_argument('-cond_cols', '--condition_columns', help = 'the columns from the sample table to be included in the visualization')
parser.add_argument('config_file', help = 'the config file used to analyze the dataset')

args = parser.parse_args()

# Function definitions
def get_all_lane_ids(sample_table):
    return np.unique(np.array(sample_table['lane']))

#def get_sample_table(config_params):
#
#    filename = config_params['sample_table_file']
#
#    # Read everything in as a string, to prevent vexing
#    # number interpretation problems! Methods further down
#    # can coerce to different types.
#    tab = pd.read_table(filename, dtype = 'S')
#    return tab

###########################################
#######  Here is the main script  #########
###########################################

# Get the config file, which is the only argument needed for the pipeline
config_file = args.config_file
config_params = parse_yaml(config_file)

# Read in the sample table
sample_table = get_sample_table(config_params)

# Grab all of the ids of the lanes to process
lane_ids = get_all_lane_ids(sample_table)
lane_ids = np.append(lane_ids, ['all_lanes', 'all_lanes_filtered'])

# First, get one strain X condition count matrix per lane
# This only needs to be run once, unless the barcodes
# or index tags change for some reason.
for lane_id in lane_ids:
    print 'clustering {} count matrix...'.format(lane_id)
    clus_wrap.cluster_count_matrix(config_file, lane_id, args.strain_columns, args.condition_columns)

