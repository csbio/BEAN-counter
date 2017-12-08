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

# Specific imports
from scipy.stats import spearmanr

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


def get_aligned_new_and_old_datasets():

    '''
    Returns the old and new matrices only where they overlap
    (both rows and columns), so that the data comparison tests don't error out.
    This should only error out if there is no overlap whatsoever.
    '''
    
    # Align the rows
    new_row_ids = [x for x in res_count_dataset[0]]
    old_row_ids = [x for x in orig_count_dataset[0]]
    common_rows = list(set(new_row_ids).intersection(set(old_row_ids)))

    new_matrix_row_inds = [new_row_ids.index(x) for x in common_rows]
    old_matrix_row_inds = [old_row_ids.index(x) for x in common_rows]

    # Align the columns
    new_col_tuples = [tuple(x) for x in res_count_dataset[1]]
    old_col_tuples = [tuple(x) for x in orig_count_dataset[1]]
    common_cols = list(set(new_col_tuples).intersection(set(old_col_tuples)))

    new_matrix_col_inds = [new_col_tuples.index(x) for x in common_cols]
    old_matrix_col_inds = [old_col_tuples.index(x) for x in common_cols]

    # Get the new aligned matrices
    aligned_old_mat = orig_count_dataset[2][np.ix_(old_matrix_row_inds, old_matrix_col_inds)]
    aligned_new_mat = res_count_dataset[2][np.ix_(new_matrix_row_inds, new_matrix_col_inds)]

    print "Original matrix shape: ", orig_count_dataset[2].shape
    print "New matrix shape: ", res_count_dataset[2].shape
    print "Aligned matrix shape: ", (len(common_rows), len(common_cols))

    return [common_rows, common_cols, aligned_new_mat], [common_rows, common_cols, aligned_old_mat]

# Get intersected versions of old and new matrices
aligned_new_dataset, aligned_old_dataset = get_aligned_new_and_old_datasets()


def test_count_mat_vals_abs_error():

    '''
    Compares the values in the new "results" count matrix to those in
    the original, reference count matrix, using the "aligned" aka
    "intersected" of the old/new matrices.
    '''
    abs_error = np.abs(aligned_new_dataset[2] - aligned_old_dataset[2])

    # Test for < +/- 20 count error of All values
    assert np.all(abs_error <= 20), "{} count values differed by more than 20 counts".format(np.sum(abs_error > 20))

    # Test for < +/- 10 count error of All values
    assert np.all(abs_error <= 10), "{} count values differed by more than 10 counts".format(np.sum(abs_error > 10))

    # Test for < +/- 5 count error of All values
    assert np.all(abs_error <= 5), "{} count values differed by more than 5 counts".format(np.sum(abs_error > 5))

    # Test for < +/- 2 count error of All values
    assert np.all(abs_error <= 2), "{} count values differed by more than 2 counts".format(np.sum(abs_error > 2))
    
    # Test for < +/- 1 count error of All values
    assert np.all(abs_error <= 1), "{} count values differed by more than 1 counts".format(np.sum(abs_error > 1))
    
    # Test for equality of ALL of the values
    assert np.all(np.isclose(aligned_new_dataset[2], aligned_old_dataset[2]))


def test_count_mat_profiles_spearman():
    
    # Get Spearman rank correlation coefficient of each old/new pair
    spearman_coefs = np.array([spearmanr(aligned_new_dataset[2][:, i], aligned_old_dataset[2][:, i])[0] for i in range(aligned_old_dataset[2].shape[1])])

    # Test that each new count profile has a spearman rank correlation > 0.80 to the old count profile
    assert np.all(spearman_coefs >= 0.80), "{} count profiles possessed spearman correlations < 0.8".format(np.sum(spearman_coefs > 0.8))
    
    # Test that each new count profile has a spearman rank correlation > 0.90 to the old count profile
    assert np.all(spearman_coefs >= 0.90), "{} count profiles possessed spearman correlations < 0.9".format(np.sum(spearman_coefs > 0.9))
    
    # Test that each new count profile has a spearman rank correlation > 0.95 to the old count profile
    assert np.all(spearman_coefs >= 0.95), "{} count profiles possessed spearman correlations < 0.95".format(np.sum(spearman_coefs > 0.95))

    # Test that each new count profile has a spearman rank correlation > 0.99 to the old count profile
    assert np.all(spearman_coefs >= 0.99), "{} count profiles possessed spearman correlations < 0.99".format(np.sum(spearman_coefs > 0.99))



