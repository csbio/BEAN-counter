#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.2.2'

# This script converts all dumped z-score matrices into CDTs,
# which are visualized using Java Treeview. I anticipate
# future versions will be smart and avoid re-clustering
# datasets if the data themselves have not changed.

import pandas as pd
import numpy as np
import scipy
import sys, os, gzip
import tarfile
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import itertools as it
import cPickle
import shutil

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'scripts'))
sys.path.append(os.path.join(barseq_path, 'lib'))

import compressed_file_opener as cfo
import cg_file_tools as cg_file
import cluster_dataset_wrappers as clus_wrap
from cg_common_functions import get_verbosity, get_sample_table, get_temp_clustergram_name, parse_yaml
from version_printing import update_version_file

import argparse

# Parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('config_file', help = 'The config file used to analyze the dataset')
parser.add_argument('strain_columns', help = 'The columns from the barcode table to be included in the visualization')
parser.add_argument('condition_columns', help = 'The columns from the sample table to be included in the visualization')

args = parser.parse_args()

# Function definitions
def get_all_lane_ids(sample_table):
    return np.unique(np.array(sample_table['lane']))

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
lane_ids = np.append(lane_ids, ['all_lanes_filtered'])

# First, get one strain X condition count matrix per lane
# This only needs to be run once, unless the barcodes
# or index tags change for some reason.
filenames_noext = []
for lane_id in lane_ids:
    print 'clustering {} z-score matrix...'.format(lane_id)
    filename_noext = clus_wrap.cluster_zscore_matrix(config_file, lane_id, args.strain_columns, args.condition_columns)
    filenames_noext.append(filename_noext)

###### Ultimate goal: combine all z-score CDTs into one big tarred/gzipped folder, for distribution!
# First, get all of the filenames (cdt, atr, gtr)
filenames = []
for f_noext in filenames_noext:
    for ext in ['cdt', 'atr', 'gtr']:
        filenames.append('{}.{}'.format(f_noext, ext))
if get_verbosity(config_params) >= 2:
    print filenames
# Now, create a temporary folder to hold all of the clustergrams for tarring/gzipping
# This temp folder name is also the name of the tarred archive. Cool, right?!
tmp_dir = get_temp_clustergram_name(config_params['output_folder'], 'all-zscore-clustergrams')
if get_verbosity(config_params) >= 2:
    print tmp_dir
# Since it's timestamped, directory should never exist with this exact name
os.makedirs(tmp_dir)
# And copy all of the files to that temp directory!
for f in filenames:
    shutil.copy(f, tmp_dir)

# Add version to that directory
update_version_file(tmp_dir, VERSION)

# Now tar the directory!
tar_fname = '{}.tar.gz'.format(tmp_dir)
tar = tarfile.open(tar_fname, 'w:gz')
tar.add(tmp_dir, arcname = os.path.basename(tmp_dir))
tar.close()
# And, remove the temporary directory!
shutil.rmtree(tmp_dir, ignore_errors = True)


