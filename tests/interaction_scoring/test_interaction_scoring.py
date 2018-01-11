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
sys.path.append(os.path.join(barseq_path, 'scripts'))
sys.path.append(os.path.join(barseq_path, 'lib'))

from testing_lib import get_data_dir, get_results_dir, copy_data_to_results_dir
from cg_common_functions import a_is_row_in_b
import counts_to_zscores
from compare import compare_2d_datasets

# Specific imports
from scipy.stats import pearsonr

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

config_file = os.path.join(results_dir, 'config_files', 'config_file_filtered.yaml')

lane = 'lane01'
counts_to_zscores.main(config_file, lane)

with gzip.open(os.path.join(data_dir, 'output', 'interactions', lane, '{}_scaled_dev.dump.gz'.format(lane))) as orig_int_score_mat_f:
    orig_int_score_dataset = cPickle.load(orig_int_score_mat_f)

with gzip.open(os.path.join(results_dir, 'output', 'interactions', lane, '{}_scaled_dev.dump.gz'.format(lane))) as res_int_score_mat_f:
    res_int_score_dataset = cPickle.load(res_int_score_mat_f)

reports_dir = os.path.join(results_dir, 'reports')
os.makedirs(reports_dir)
aligned_old_dataset, aligned_new_dataset, comparison_stats, align_stats = compare_2d_datasets(orig_int_score_dataset,
        res_int_score_dataset, reports_dir, 'interaction-scoring')

def test_rows_align():
    '''
    Ensures the rows align (should not be
    an issue, but just making sure)
    '''
    n_mismatching_rows = np.array([len(align_stats.get(x)) for x in ['rows_1_not_2', 'rows_2_not_1']])
    assert all(n_mismatching_rows == 0),'\n{} rows in the reference matrix '\
            'and {} rows in the new matrix are not found in the other'.format(*n_mismatching_rows)

def test_cols_align():
    '''
    Ensures the cols align (should not be
    an issue, but just making sure)
    '''
    n_mismatching_cols = np.array([len(align_stats.get(x)) for x in ['cols_1_not_2', 'cols_2_not_1']])
    assert all(n_mismatching_cols == 0),'\n{} cols in the reference matrix '\
            'and {} cols in the new matrix are not found in the other'.format(*n_mismatching_cols)

def test_abs_error():

    '''
    Compares the values in the new "results" z-score matrix to those in
    the original, reference z-score matrix, using the "aligned" aka
    "intersected" of the old/new matrices.
    '''
    assert comparison_stats['all_close']


#def test_int_score_mat_profiles_pearson():
#    
#    # Get Pearson correlation coefficient of each old/new pair
#    pearson_coefs = np.array([pearsonr(res_int_score_dataset[2][:, i], orig_int_score_dataset[2][:, i])[0] for i in range(orig_int_score_dataset[2].shape[1])])
#    print pearson_coefs
#
#    # Test that each new z-score profile has a pearson correlation > 0.80 to the old z-score profile
#    assert np.all(pearson_coefs >= 0.80), "{} z-score profiles possessed pearson correlations < 0.8".format(np.sum(pearson_coefs > 0.8))
#    
#    # Test that each new z-score profile has a pearson correlation > 0.90 to the old z-score profile
#    assert np.all(pearson_coefs >= 0.90), "{} z-score profiles possessed pearson correlations < 0.9".format(np.sum(pearson_coefs > 0.9))
#    
#    # Test that each new z-score profile has a pearson correlation > 0.95 to the old z-score profile
#    assert np.all(pearson_coefs >= 0.95), "{} z-score profiles possessed pearson correlations < 0.95".format(np.sum(pearson_coefs > 0.95))
#
#    # Test that each new z-score profile has a pearson correlation > 0.99 to the old z-score profile
#    assert np.all(pearson_coefs >= 0.99), "{} z-score profiles possessed pearson correlations < 0.99".format(np.sum(pearson_coefs > 0.99))



