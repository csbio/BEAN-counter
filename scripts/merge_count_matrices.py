#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script reads in all of the count matrices from a screen and
# merges them together. Filtering steps are performed in the next script

import pandas as pd
import numpy as np
import scipy
import sys, os, gzip
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import gzip
import cPickle

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append(os.path.join(barseq_path, 'lib'))

import config_file_parser as cfp
import compressed_file_opener as cfo
import cg_file_tools as cg_file
from cg_common_functions import get_verbosity, get_sample_table

#def get_sample_table(config_params):
#
#    filename = config_params['sample_table_file']
#
#    # Read everything in as a string, to prevent vexing
#    # number interpretation problems! Methods further down
#    # can coerce to different types.
#    tab = pd.read_table(filename, dtype = 'S')
#    return tab

def get_all_lane_ids(sample_table):

    return np.unique(np.array(sample_table['lane']))

def get_lane_data_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'intermediate', lane_id)

def get_dumped_count_matrix_filename(config_params, lane_id):

    lane_folder = get_lane_data_path(config_params, lane_id)
    return os.path.join(lane_folder, '{}_barseq_matrix.dump.gz'.format(lane_id))

def combine_count_matrices(config_params):

    sample_table = get_sample_table(config_params)
    all_lane_ids = get_all_lane_ids(sample_table)

    count_matrix_list = []
    condition_id_list = []

    for lane_id in all_lane_ids:
        count_matrix_filename = get_dumped_count_matrix_filename(config_params, lane_id)
        f = gzip.open(count_matrix_filename)

        gene_barcode_ids, condition_ids, count_matrix = cPickle.load(f)
        f.close()
        count_matrix_list.append(count_matrix)
        condition_id_list.append(condition_ids)

    all_count_matrix = np.hstack(count_matrix_list)
    all_condition_ids = np.vstack(condition_id_list)

    # Note: for all matrices, the gene_barcode_ids, should be the same. When generating the count matrix,
    # I do not let any gene_barcode_ids disappear, even if they have no counts whatsoever. This should
    # simplify initial processing steps, and the filtering does not need to happen until later.
    return gene_barcode_ids, all_condition_ids, all_count_matrix

def dump_dataset(dataset, filename):

    f = gzip.open(filename, 'wb')
    cPickle.dump(dataset, f)
    f.close()

def main(config_file):

    # Read in the config params
    config_params = cfp.parse(config_file)
    sample_table = get_sample_table(config_params)

    # Read in all of the z-score matrices and combine into one matrix
    dataset = combine_count_matrices(config_params)

    # Get a new folder to house the combined count matrix
    combined_count_folder = get_lane_data_path(config_params, 'all_lanes')
    if not os.path.isdir(combined_count_folder):
        os.makedirs(combined_count_folder)

    # Dump out the combined count matrix!
    combined_count_filename = get_dumped_count_matrix_filename(config_params, 'all_lanes')
    dump_dataset(dataset, combined_count_filename)

# call: python merge_count_matrices.py <config_file>
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: python merge_count_matrices.py <config_file>'
    else:
        config_file = sys.argv[1]
        main(config_file)
