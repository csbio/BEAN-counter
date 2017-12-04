#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

import os, sys
import gzip, cPickle
import numpy as np

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."

sys.path.append(os.path.join(barseq_path, 'master_scripts'))
sys.path.append(os.path.join(barseq_path, 'lib'))

from testing_lib import get_data_dir, get_results_dir, copy_data_to_results_dir
from cg_common_functions import a_is_row_in_b
import process_screen

# How do I make sure I'm running from the correct directory?!
# Or do I find a way to manipulate the config file...
#def test_filename():
#    print "running test_filename"
#    testing_lib.get_testing_dirname()
#    print __file__

data_dir = get_data_dir()
results_dir = get_results_dir()

copy_data_to_results_dir(data_dir, results_dir)

os.chdir(results_dir)

config_file = os.path.join(results_dir, 'config_files', 'config_file.yaml')

process_screen.main(config_file, 1, 1)

with gzip.open(os.path.join(data_dir, 'output', 'intermediate', 'lane01', 'lane01_barseq_matrix.dump.gz')) as orig_count_mat_f:
    orig_count_dataset = cPickle.load(orig_count_mat_f)

with gzip.open(os.path.join(results_dir, 'output', 'intermediate', 'lane01', 'lane01_barseq_matrix.dump.gz')) as res_count_mat_f:
    res_count_dataset = cPickle.load(res_count_mat_f)

def test_count_mat_rows():

    '''
    Ensures that all rows in the new "results" count matrix are found in
    the original, reference count matrix.
    '''
    assert all(x in orig_count_dataset[0] for x in res_count_dataset[0])

def test_count_mat_cols():

    '''
    Ensures that all columns in the new "results" count matrix are found in
    the original, reference count matrix.
    '''
    assert all(a_is_row_in_b(x, orig_count_dataset[1]) for x in res_count_dataset[1])

def test_count_mat_vals():

    '''
    Ensures that all values in the new "results" count matrix are found in
    the original, reference count matrix.
    '''
    # First, line up rows and columns of the matrix
    # If there's an error stating that only length one arrays can be converted to
    # Python scalars, then that's because there's a row ID in the new result that
    # wasn't there before.
    row_inds = np.array([int(np.where(np.array(res_count_dataset[0]) == x)[0]) for x in orig_count_dataset[0]])
    col_inds = np.array([i for i, x in enumerate(res_count_dataset[1]) for y in orig_count_dataset[1] if all(x == y) ])

    print "row inds: ", row_inds
    print "col inds: ", col_inds

    print "row inds shape: ", row_inds.shape
    print "col inds shape: ", col_inds.shape

    res_mat_reordered = res_count_dataset[2][np.ix_(row_inds, col_inds)]

    assert np.all(np.isclose(res_mat_reordered, orig_count_dataset[2]))
    #assert all(a_is_row_in_b(x, orig_count_dataset[1]) for x in res_count_dataset[1])

#def test_folder():
#    print data_dir
#    print results_dir



