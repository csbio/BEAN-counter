#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.2.1'

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
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'lib'))

import compressed_file_opener as cfo
import cg_file_tools as cg_file
from cg_common_functions import get_verbosity, get_sample_table, parse_yaml
from version_printing import update_version_file

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
    gene_barcode_id_list = []
    all_gene_barcode_ids = set()

    for lane_id in all_lane_ids:
        count_matrix_filename = get_dumped_count_matrix_filename(config_params, lane_id)
        f = gzip.open(count_matrix_filename)

        gene_barcode_ids, condition_ids, count_matrix = cPickle.load(f)
	print "number of nonzero barcodes: {}".format(len(gene_barcode_ids))

        f.close()
        count_matrix_list.append(count_matrix)
        condition_id_list.append(condition_ids)
        gene_barcode_id_list.append(gene_barcode_ids)
        all_gene_barcode_ids.update(gene_barcode_ids)

    # Because these matrices may not share all of the same rows
    # in the same order, here I must align the rows of the
    # matrices and add missing rows if needed. Using zeros to
    # fill the rows in some matrices that didn't have any counts.
    all_gene_barcode_ids = list(all_gene_barcode_ids)
    all_gene_barcode_indices = {x:i for i,x in enumerate(all_gene_barcode_ids)}
    aligned_count_matrix_list = []
    for i in range(len(count_matrix_list)):
        aligned_count_matrix_list.append(np.zeros((len(all_gene_barcode_ids), len(condition_id_list[i]))))
        orig_barcodes = gene_barcode_id_list[i]
        for j, barcode in enumerate(orig_barcodes):
            aligned_count_matrix_list[i][all_gene_barcode_indices[barcode], :] = count_matrix_list[i][j, :]

    all_count_matrix = np.hstack(aligned_count_matrix_list)
    all_condition_ids = np.vstack(condition_id_list)

    print all_count_matrix.shape
    assert len(all_gene_barcode_ids) == all_count_matrix.shape[0], "Number of barcodes does not match the number of rows in the matrix"
    assert len(all_condition_ids) == all_count_matrix.shape[1], "Number of conditions does not match the number of columns in the matrix"

    return all_gene_barcode_ids, all_condition_ids, all_count_matrix

def dump_dataset(dataset, filename):

    f = gzip.open(filename, 'wb')
    cPickle.dump(dataset, f)
    f.close()

def main(config_file):

    # Read in the config params
    config_params = parse_yaml(config_file)
    sample_table = get_sample_table(config_params)

    # Read in all of the z-score matrices and combine into one matrix
    dataset = combine_count_matrices(config_params)

    # Get a new folder to house the combined count matrix
    combined_count_folder = get_lane_data_path(config_params, 'all_lanes')
    if not os.path.isdir(combined_count_folder):
        os.makedirs(combined_count_folder)
    update_version_file(combined_count_folder, VERSION)

    # Dump out the combined count matrix!
    combined_count_filename = get_dumped_count_matrix_filename(config_params, 'all_lanes')
    dump_dataset(dataset, combined_count_filename)
    
    update_version_file(config_params['output_folder'], VERSION)

# call: python merge_count_matrices.py <config_file>
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: python merge_count_matrices.py <config_file>'
    else:
        config_file = sys.argv[1]
        main(config_file)
