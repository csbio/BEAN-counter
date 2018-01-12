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
import batch_correction
from compare import compare_2d_datasets

# Specific imports
from scipy.stats import spearmanr

# Ensures that we're working with the right folders
data_dir = get_data_dir()
results_dir = get_results_dir()
copy_data_to_results_dir(data_dir, results_dir)
os.chdir(results_dir)

# Gets correct arguments to batch_correction.py
dataset_file = os.path.join(data_dir, 'batch_effect_dataset.dump.gz')
sample_table_file = os.path.join(data_dir, 'sample_table.txt')
batch_col = "lane"
nondup_cols = "name".split(',')
max_comps = 1
output_folder = os.path.abspath(results_dir)
dataset_file = os.path.abspath(dataset_file)
verbosity = 1

# Runs batch_correction.py
dataset = batch_correction.load_dataset(dataset_file)
sample_table = batch_correction.read_sample_table(sample_table_file)
batch_correction.main(dataset, sample_table, batch_col, nondup_cols, max_comps, output_folder, dataset_file, verbosity)

with gzip.open(os.path.join(data_dir, 'batch_corrected_datasets_matlab.dump.gz')) as orig_batch_mat_f:
    orig_batch_data = cPickle.load(orig_batch_mat_f)
orig_batch_data = orig_batch_data[1:4]
orig_batch_data_1 = [orig_batch_data[0], orig_batch_data[1], orig_batch_data[2][0,:,:]]
orig_batch_data_2 = [orig_batch_data[0], orig_batch_data[1], orig_batch_data[2][1,:,:]]

with gzip.open(os.path.join(results_dir, 'batch_corrected_datasets.dump.gz')) as res_batch_mat_f:
    res_batch_data = cPickle.load(res_batch_mat_f)
res_batch_data = res_batch_data[1:4]
res_batch_data_1 = [res_batch_data[0], res_batch_data[1], res_batch_data[2][0,:,:]]
res_batch_data_2 = [res_batch_data[0], res_batch_data[1], res_batch_data[2][1,:,:]]
	
reports_dir = os.path.join(results_dir, 'reports')
os.makedirs(reports_dir)
aligned_old_dataset_1, aligned_new_dataset_1, comparison_stats_1, align_stats_1 = compare_2d_datasets(orig_batch_data_1,
        res_batch_data_1, reports_dir, 'orig_data')
aligned_old_dataset_2, aligned_new_dataset_2, comparison_stats_2, align_stats_2 = compare_2d_datasets(orig_batch_data_2,
        res_batch_data_2, reports_dir, '1_comp_removed', new_log = False)

def test_rows_align():
    '''
    Ensures the rows align (should not be
    an issue, but just making sure)
    '''
    n_mismatching_rows = np.array([len(align_stats_1.get(x)) for x in ['rows_1_not_2', 'rows_2_not_1']])
    assert all(n_mismatching_rows == 0),'\n{} rows in the reference matrix '\
            'and {} rows in the new matrix are not found in the other'.format(*n_mismatching_rows)
    n_mismatching_rows = np.array([len(align_stats_2.get(x)) for x in ['rows_1_not_2', 'rows_2_not_1']])
    assert all(n_mismatching_rows == 0),'\n{} rows in the reference matrix '\
            'and {} rows in the new matrix are not found in the other'.format(*n_mismatching_rows)

def test_cols_align():
    '''
    Ensures the cols align (should not be
    an issue, but just making sure)
    '''
    n_mismatching_cols = np.array([len(align_stats_1.get(x)) for x in ['cols_1_not_2', 'cols_2_not_1']])
    assert all(n_mismatching_cols == 0),'\n{} cols in the reference matrix '\
            'and {} cols in the new matrix are not found in the other'.format(*n_mismatching_cols)
    n_mismatching_cols = np.array([len(align_stats_2.get(x)) for x in ['cols_1_not_2', 'cols_2_not_1']])
    assert all(n_mismatching_cols == 0),'\n{} cols in the reference matrix '\
            'and {} cols in the new matrix are not found in the other'.format(*n_mismatching_cols)

def test_abs_error():

    '''
    Compares the values in the new "results" z-score matrix to those in
    the original, reference z-score matrix, using the "aligned" aka
    "intersected" of the old/new matrices.
    '''
    assert comparison_stats_1['all_close']
    assert comparison_stats_2['all_close']
