#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.2.1'

import pandas as pd
import numpy as np
import matplotlib
import os
import sys
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import gzip


def align_2d_datasets(dataset_1, dataset_2):
    '''
    Returns the intersection of both datasets, aligned row-
    and columnwise.
    '''
    
    # Align the rows
    row_ids_1 = [x for x in dataset_1[0]]
    row_ids_2 = [x for x in dataset_2[0]]
    common_rows = list(set(row_ids_1).intersection(set(row_ids_2)))

    row_inds_1 = [row_ids_1.index(x) for x in common_rows]
    row_inds_2 = [row_ids_2.index(x) for x in common_rows]

    # Align the columns
    col_tuples_1 = [tuple(x) for x in dataset_1[1]]
    col_tuples_2 = [tuple(x) for x in dataset_2[1]]
    common_cols = list(set(col_tuples_1).intersection(set(col_tuples_2)))

    col_inds_1 = [col_tuples_1.index(x) for x in common_cols]
    col_inds_2 = [col_tuples_2.index(x) for x in common_cols]

    # Get the new aligned matrices
    aligned_mat_1 = dataset_1[2][np.ix_(row_inds_1, col_inds_1)]
    aligned_mat_2 = dataset_2[2][np.ix_(row_inds_2, col_inds_2)]

    # Take care of no overlap case - I think it's taken care of by default...

    # Compute/compile alignment statistics
    align_stats = {
            'n_rows_1': dataset_1[2].shape[0],
            'n_rows_2': dataset_2[2].shape[0],
            'rows_1_not_2': set(row_ids_1).difference(set(row_ids_2)),
            'rows_2_not_1': set(row_ids_2).difference(set(row_ids_1)),
            'n_cols_1': dataset_1[2].shape[1],
            'n_cols_2': dataset_2[2].shape[1],
            'cols_1_not_2': set(row_ids_1).difference(set(row_ids_2)),
            'cols_2_not_1': set(row_ids_2).difference(set(row_ids_1))
            }

    return [common_rows, common_cols, aligned_mat_1], [common_rows, common_cols, aligned_mat_2], align_stats

def compare_2d_datasets(dataset_1, dataset_2, output, prefix, new_log = True):
    '''
    Takes in two datasets, prints out both text and plot-
    based reports to the folder "output", with filenames
    prepended with "prefix."
    '''

    aligned_dataset_1, aligned_dataset_2, align_stats = align_2d_datasets(dataset_1, dataset_2)

    if not os.path.isdir(output):
        os.makedirs(output)
    
    log = None
    log_file = os.path.join(output, 'report.txt')
    if os.path.exists(log_file) and new_log is False:
        log = open(log_file, 'a')
    else:
        log = open(log_file, 'w+')

    # Get stats on nans, then replace nans with zeros
    aligned_dims = aligned_dataset_1[2].shape
    nan_inds_1 = np.where(np.isnan(aligned_dataset_1[2]))
    nan_inds_2 = np.where(np.isnan(aligned_dataset_2[2]))
    nan_inds_1_flat = np.ravel_multi_index(nan_inds_1, aligned_dims)
    nan_inds_2_flat = np.ravel_multi_index(nan_inds_2, aligned_dims)
    nan_inds_1_not_2 = np.unravel_index(np.setdiff1d(nan_inds_1_flat, nan_inds_2_flat), dims = aligned_dims)
    nan_inds_2_not_1 = np.unravel_index(np.setdiff1d(nan_inds_2_flat, nan_inds_1_flat), dims = aligned_dims)

    matrix_1 = aligned_dataset_1[2].copy()
    matrix_2 = aligned_dataset_2[2].copy()
   
    matrix_1[np.isnan(matrix_1)] = 0
    matrix_2[np.isnan(matrix_2)] = 0
    
    # Get some stats on individual value differences
    # Error is relative to matrix 1
    abs_error_mat = matrix_2 - matrix_1
    rel_error_mat = abs_error_mat / matrix_1
    log10_abs_rel_error_mat = np.log10(np.abs(rel_error_mat))
    #print np.sum(np.isinf(rel_error_mat))
    #print np.sum(np.isnan(rel_error_mat))
    #print np.unique(np.array(rel_error_mat[~np.isfinite(rel_error_mat)], dtype = np.str))
    
    mat_corr_pearson = pearsonr(matrix_1.flatten(), matrix_2.flatten())[0]
    mat_corr_spearman = spearmanr(matrix_1.flatten(), matrix_2.flatten())[0]
    col_corrs_pearson = np.array([pearsonr(matrix_1[:, i], matrix_2[:, i])[0] for i in range(matrix_1.shape[1])])
    col_corrs_spearman = np.array([spearmanr(matrix_1[:, i], matrix_2[:, i])[0] for i in range(matrix_1.shape[1])])

    export_histogram(abs_error_mat, 'Absolute error', 'Frequency', output, prefix, 'abs-error')
    export_histogram(rel_error_mat, 'Relative error', 'Frequency', output, prefix, 'rel-error')
    export_histogram(log10_abs_rel_error_mat, 'log10 ( abs ( relative error ) )', 'Frequency', output, prefix, 'log10-abs-rel-error')

    # Summarize individual value differences per column
    # Must work in absolute value space so differences don't cancel
    abs_error_columnwise_abs_mean = np.nanmean(np.abs(abs_error_mat), axis = 0)
    rel_error_columnwise_abs_mean = np.array([np.mean(np.abs(x[np.isfinite(x)])) for x in rel_error_mat.T])
    #print rel_error_columnwise_mean

    export_histogram(abs_error_columnwise_abs_mean, 'Columnwise mean of abs( absolute error )', 'Frequency', output, prefix, 'abs-error-columnwise-abs-mean')
    export_histogram(rel_error_columnwise_abs_mean, 'Columnwise mean of abs( relative error )', 'Frequency', output, prefix, 'rel-error-columnwise-abs-mean')

    # Plot relative error vs. absolute error
    fig = plt.figure()
    plt.scatter(abs_error_mat, log10_abs_rel_error_mat, s = 1, marker = '.')
    plt.xlabel('Absolute error')
    plt.ylabel('log10 ( abs ( relative error ) )')
    plt.title(prefix)
    fig.savefig(os.path.join(output, prefix + '_rel-error_vs_abs-error.png'))
    plt.close()
    
    # Get stats on Pearson/Spearman corrs between the datasets
    export_histogram(col_corrs_pearson, 'Columnwise Pearson correlation coefficient', 'Frequency', output, prefix, 'columnwise_pearson', x_range = (-1, 1))
    
    export_histogram(col_corrs_spearman, 'Columnwise Spearman correlation coefficient', 'Frequency', output, prefix, 'columnwise_spearman', x_range = (-1, 1))
    
    # Compute/compile comparison stats
    comparison_stats = {
            'n_nans_1' : len(nan_inds_1[0]),
            'n_nans_2' : len(nan_inds_2[0]),
            'n_nans_in_1_not_2' : len(nan_inds_1_not_2[0]),
            'n_nans_in_2_not_1' : len(nan_inds_2_not_1[0]),
            'n_vals' : matrix_1.size,
            'n_equiv_vals' : np.sum(np.isclose(matrix_1, matrix_2)),
            'all_close' : np.allclose(matrix_1, matrix_2),
            'matrix_pearson': mat_corr_pearson,
            'matrix_spearman': mat_corr_spearman
            }

    print >>log, '=============================================='
    print >>log, prefix
    print >>log, '----------------------------------------------'
    print >>log, '\n'

    print >>log, 'matrix alignment statistics'
    print >>log, '----------------------------------------------'
    print >>log, 'num rows in matrix 1:', align_stats['n_rows_1']
    print >>log, 'num rows in matrix 2:', align_stats['n_rows_2']
    print >>log, 'num rows in matrix 1 but not 2:', len(align_stats['rows_1_not_2'])
    print >>log, 'num rows in matrix 2 but not 1:', len(align_stats['rows_2_not_1'])
    print >>log, 'num cols in matrix 1:', align_stats['n_cols_1']
    print >>log, 'num cols in matrix 2:', align_stats['n_cols_2']
    print >>log, 'num cols in matrix 1 but not 2:', len(align_stats['cols_1_not_2'])
    print >>log, 'num cols in matrix 2 but not 1:', len(align_stats['cols_2_not_1'])
    print >>log, '\n'

    print >>log, 'matrix comparison statistics'
    print >>log, '(performed on aligned matrices)'
    print >>log, '----------------------------------------------'
    print >>log, 'num values in each matrix:', comparison_stats['n_vals']
    print >>log, 'num equivalent values between matrices 1 and 2:', comparison_stats['n_equiv_vals']
    print >>log, 'matrix correlation, Pearson:', comparison_stats['matrix_pearson']
    print >>log, 'matrix correlation, Spearman:', comparison_stats['matrix_spearman']
    print >>log, 'num nans in matrix 1:', comparison_stats['n_nans_1']
    print >>log, 'num nans in matrix 2:', comparison_stats['n_nans_2']
    print >>log, 'num nans in matrix 1 but not 2:', comparison_stats['n_nans_in_1_not_2']
    print >>log, 'num nans in matrix 2 but not 1:', comparison_stats['n_nans_in_2_not_1']
    
    return aligned_dataset_1, aligned_dataset_2, comparison_stats, align_stats

def export_histogram(x, xlabel, ylabel, output, prefix, fname, bins = 400, x_range = None):
    freq, endpoints = np.histogram(x[np.isfinite(x)], bins = bins, range = x_range)
    with open(os.path.join(output, '{}_{}.txt'.format(prefix, fname)), 'wt') as txt_f:
        txt_f.write('bin_lower\tbin_upper\tcounts\n')
        for i in range(len(freq)):
            txt_f.write('{}\t{}\t{}\n'.format(endpoints[i], endpoints[i+1], freq[i]))
    fig = plt.figure()
    plt.bar((endpoints[:-1] + endpoints[1:]) / 2, freq,
            align = 'center', width = (endpoints[1] - endpoints[0]))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(prefix)
    fig.savefig(os.path.join(output, '{}_{}.png'.format(prefix, fname)))
    plt.close()



## Compares two 2d matrices
#def compare_2d(mat1, mat2, log_file, output, count=1):
#
#    # Opens log to write to and creates it if necessary
#    log = None
#    if os.path.exists(log_file):
#        log = file(log_file, "a")
#    else:
#        log = file(log_file, "w+")
#
#    # Statistics to be recorded
#    col_pearson_corr = []
#    col_diffs = []
#    col_diff_counts = []
#    no_diffs = 0
#
#    # Iterates through all columns
#    for i in range(mat1.shape[1]):
#        m_col = mat1[:, i]
#        n_col = mat2[:, i]
#        col_pearson_corr.append(stats.pearsonr(m_col, n_col)[0])
#        diff_counts = 0
#        for j in range(0, len(m_col)):
#            if not np.isclose(m_col[j], n_col[j]):
#                diff_counts += 1
#                col_diffs.append(abs(m_col[j] - n_col[j]))
#            else: no_diffs += 1
#        col_diff_counts.append(diff_counts)
#
#    # Plots pearson correlations in a histogram
#    fig = plt.figure()
#    plt.hist(col_pearson_corr, range=(0, 1), bins=200)
#    plt.xlabel("Pearson correlation coefficient b/w cols")
#    fig_name = "pearson_col_corr_" + str(count) + ".png"
#    fig.savefig(os.path.join(output, fig_name))
#
#    # Plots value differences in a histogram
#    fig = plt.figure()
#    plt.hist(col_diffs, range=(0, 10), bins=200)
#    plt.xlabel("Differences between individual values")
#    fig_name = "raw_diffs_" + str(count) + ".png"
#    fig.savefig(os.path.join(output, fig_name))
#
#    # Writes summary statistics to file
#    print >>log, "Matrix  " + str(count)
#    print >>log, ""
#    print >>log, "Total number of equivalent values: " + str(no_diffs)
#    print >>log, "Average number of non-close values in all columns: " + str(np.mean(col_diff_counts))
#    print >>log, "Average non-close value differences in all columns: " + str(np.mean(col_diffs))
#    print >>log, "Number of differences greater than 0.0001: " + str(len([d for d in col_diffs if d > 0.0001]))
#    print >>log, "Number of differences greater than 0.001: " + str(len([d for d in col_diffs if d > 0.001]))
#    print >>log, "Number of differences greater than 0.01: " + str(len([d for d in col_diffs if d > 0.01]))
#    print >>log, "Number of differences greater than 0.1: " + str(len([d for d in col_diffs if d > 0.1]))
#    print >>log, "Number of differences greater than 1: " + str(len([d for d in col_diffs if d > 1]))
#    print >>log, ""
#
#    # Closes log file
#    log.close()
#
#
#if __name__ == '__main__':
#
#    # Gets user's args
#    parser = make_arg_parser()
#    args = parser.parse_args()
#    args = vars(args)
#
#    # Makes output directory if it doesn't exist
#    if not os.path.exists(args['output']):
#        os.makedirs(args['output'])
#
#    # Deletes log file if it already exists
#    log_file = os.path.join(args['output'], "log.txt")
#    if os.path.exists(log_file):
#        os.remove(log_file)
#
#    # Loads matrices
#    m = None
#    n = None
#    with gzip.open(args['matrix1'], 'rb') as mat1:
#        m = pickle.load(mat1)
#    with gzip.open(args['matrix2'], 'rb') as mat2:
#        n = pickle.load(mat2)
#    # f = open(args['matrix1'])
#    # m = cPickle.load(f)
#    # f.close()
#    # f = open(args['matrix1'])
#    # n = cPickle.load(f)
#    # f.close()
#
#    # Only examines non-empty arrays
#    m_scores = np.asarray([np.nan_to_num(arr) for arr in m[3] if np.sum(~np.isnan(arr)) > 2])
#    n_scores = np.asarray([np.nan_to_num(arr) for arr in n[3] if np.sum(~np.isnan(arr)) > 2])
#
#    # Iterates through all 2d matrices
#    for i in range(m[3].shape[0]):
#        compare_2d(m[3][i, :, :], n[3][i, :, :], log_file, args['output'], i+1)
#

#    # Calculates Pearson p-values
#    p_vals = []
#    for i in range(len(m_scores)):
#        p = stats.ttest_ind(m_scores[i], p_scores[i])[1]
#        p_vals.append(p) 
#    print("Average Pearson p-value: " + str(np.nanmean(p_vals)))

#    # Calculates Spearman p-values
#    p_vals = []
#    for i in range(len(m_scores)):
#        p = stats.spearmanr(m_scores[i], p_scores[i])[1]
#        p_vals.append(p) 
#    print("Average Spearman p-value: " + str(np.nanmean(p_vals)))
