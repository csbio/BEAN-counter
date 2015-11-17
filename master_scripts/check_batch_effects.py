#!/usr/bin/env python
# This script takes in matrices, on which LDA has been performed and from which
# LDA components have been removed, and plots PR curves and historgrams
# visualizing the extent of the multiplex tag effect as each successive LDA
# component is removed.
import pandas as pd
import numpy as np
import scipy
from scipy.stats import rankdata
import sys, os, gzip
import jellyfish as jf
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import itertools as it
import cPickle
import argparse

from sklearn.metrics import precision_recall_curve

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append(os.path.join(barseq_path, 'lib'))

import config_file_parser as cfp
import compressed_file_opener as cfo
import cg_file_tools as cg_file
import correlation_functions as cf
from cg_common_functions import *
import plotting_tools as pt

def read_sample_table(tab_filename):

    # Read everything in as a string, to prevent vexing
    # number interpretation problems! Methods further down
    # can coerce to different types.
    tab = pd.read_table(tab_filename, dtype = 'S')
    return tab

def filter_sample_table(sample_table, final_conditions):

    all_conditions = np.array([tuple(row[['screen_name', 'expt_id']]) for i, row in sample_table.iterrows()])
    rows_to_keep = [i for i, val in enumerate(all_conditions) if a_is_row_in_b(val, final_conditions)]
    return sample_table.iloc[rows_to_keep]

def a_is_row_in_b(a, b):

    return np.any(np.all(a == b, axis = 1))

def get_lda_pr_curve_filename(folder, ymax):

    return os.path.join(folder, 'PR_batch-effect_{:.1f}.pdf'.format(ymax))

def get_batch_effect_histogram_filename(folder, ncomps):

    return os.path.join(folder, 'Histogram_batch-effect_{}-components-removed.pdf'.format(ncomps))

def load_3d_dataset(data_filename):

    f = gzip.open(data_filename, 'rb')
    ncomps, barcodes, conditions, matrix = cPickle.load(f)
    dataset = [np.array(ncomps), np.array(barcodes), np.array(conditions), matrix]
    f.close()
    return dataset

def get_labels_from_conditions(conditions, sample_table, label_name):

    sample_table = sample_table.set_index(['screen_name', 'expt_id'])
    conditions_tuple = [tuple(x) for x in conditions]
    return np.array(sample_table[label_name].ix[conditions_tuple])

def compute_PR_vectors(corr_mat, batches, verbosity):

    triu_rows, triu_cols = np.triu_indices_from(corr_mat, k = 1)

    corr_vec = corr_mat[triu_rows, triu_cols]

    batches_1 = batches[triu_rows]
    batches_2 = batches[triu_cols]
    batches_match = np.array(batches_1 == batches_2, dtype = np.int)
    # Get the number of true positives so we can use that as recall instead of
    # a fraction from 0 to 1.
    num_true_positives = np.sum(batches_match)

    if verbosity >= 1:
        print '\t\tnumber of NaN correlations: {}'.format(np.sum(np.isnan(corr_vec)))
    
    precision, recall, thresholds = precision_recall_curve(y_true = batches_match, probas_pred = corr_vec)

    # Reverse the orders so that they can easily be printed to file
    # with the highest recalls first
    return precision[::-1], recall[::-1] * num_true_positives, thresholds[::-1]

def plot_pr_curves_diff_scales(comps_removed, precision_vectors, recall_vectors, batch_column, pr_folder):

    ymax_list = [1, 0.1]
    filenames = [get_lda_pr_curve_filename(pr_folder, y) for y in ymax_list]

    for i, ymax in enumerate(ymax_list):
        plot_one_pr_curve(comps_removed, precision_vectors, recall_vectors, ymax, batch_column, filenames[i])

def plot_one_pr_curve(comps_removed, precision_vectors, recall_vectors, ymax, batch_column, filename):

    colors_lines = pt.get_style_combos(['lines', 'colors'])
    fig = plt.figure(figsize = (7, 7))
    ax = fig.add_subplot(1,1,1)
    for i, ncomps in enumerate(comps_removed):
        style = colors_lines.next()
        ax.plot(recall_vectors[i], precision_vectors[i], style, label = ncomps)
    ax.set_xscale('log')
    ax.set_ylim([0, ymax])
    ax.set_ylabel('Precision')
    ax.set_xlabel('Recall')
    lgd = ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, title = "Number of components\nremoved")
    plt.title('Precision-recall analysis of\n"{}" batch effects'.format(batch_column))
    plt.savefig(filename, bbox_extra_artists = ([lgd]), bbox_inches = 'tight')

def plot_histograms(within_batch_corrs, between_batch_corrs, ncomps, batch_column, hist_folder):

    filename = get_batch_effect_histogram_filename(hist_folder, ncomps)

    hist_bins = np.linspace(-1, 1, 41)
    bin_centers = 0.5 * (hist_bins[1:] + hist_bins[:-1])
    y_within_batch, bins = np.histogram(within_batch_corrs, bins = hist_bins, density = True)
    y_between_batch, bins = np.histogram(between_batch_corrs, bins = hist_bins, density = True)

    mean_within_batch_cor = np.nanmean(within_batch_corrs)
    mean_between_batch_cor = np.nanmean(between_batch_corrs)

    colors = pt.color_cycle()
    fig = plt.figure(figsize = (7, 7))
    ax = fig.add_subplot(1,1,1)
    color = colors.next()
    ax.plot(bin_centers, y_within_batch, color + '-', label = 'Within batch (mean = {:.2f})'.format(mean_within_batch_cor))
    color = colors.next()
    ax.plot(bin_centers, y_between_batch, color + '-', label = 'Between batch (mean = {:.2f})'.format(mean_between_batch_cor))
    ax.set_xlim([-1, 1])
    ax.set_xlabel('Average correlation within/between batch(es)')
    ax.set_ylabel('Density')
    lgd = ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, title = "Correlation")
    plt.title('{} batch correlations,\n{} components removed'.format(batch_column, ncomps))
    plt.savefig(filename, bbox_extra_artists = ([lgd]), bbox_inches = 'tight')

def print_pr_values(precision, recall, threshold, ncomps, pr_folder):

    pr_data_folder = os.path.join(pr_folder, 'PR_data')
    if not os.path.isdir(pr_data_folder):
        os.makedirs(pr_data_folder)
    filename = os.path.join(pr_data_folder, 'PR_{}-components-removed.txt.gz'.format(ncomps))

    # Make sure sure threshold vector is same length as recall and precision vectors
    threshold = np.append(threshold, np.nan)

    f = gzip.open(filename, 'wb')
    f.write('Recall\tPrecision\tThreshold\n')
    for i, rec in enumerate(recall):
        f.write('{}\t{}\t{}\n'.format(rec, precision[i], threshold[i]))
    f.close()

def print_corrs(groups, corrs, n, string, ncomps, hist_folder):

    corr_folder = os.path.join(hist_folder, 'correlation_data')
    if not os.path.isdir(corr_folder):
        os.makedirs(corr_folder)
    filename = os.path.join(corr_folder, '{}_correlations_{}-components-removed.txt'.format(string, ncomps))

    # Sort the correlations and the groups in descending order!
    # Sorting using the negative of corrs forces NaNs to stay
    # at the end of the vector despite the otherwise reversed
    # sorting order (descending)
    sort_inds = np.argsort(-corrs)
    corrs_sorted = corrs[sort_inds]
    groups_sorted = [tuple(x) for x in np.array(groups)[sort_inds]]
    n_sorted = n[sort_inds]

    f = open(filename, 'wt')
    f.write('batch_1\tbatch_2\tavg_correlation\tnum_non-NaN-correlations\n')
    
    for i, group in enumerate(groups_sorted):
        f.write('{}\t{}\t{}\t{}\n'.format(group[0], group[1], corrs_sorted[i], n_sorted[i]))

    f.close()
    return None

def main(dataset_3d, sample_table, batch_column, output_folder, verbosity, num_test_batches):

    # Load in the entire 3-d stacked matrix and dimension names!
    components_removed, barcodes, conditions, matrix_3d = dataset

    # Filter the sample table, so that no samples that have been eliminated
    # are used in further processing
    sample_table = filter_sample_table(sample_table, conditions)
   
    if verbosity >= 2:
        print sample_table

    # Get a vector of the batches that lines up with the vectors of screen_name, expt_id, name, etc
    batches = get_labels_from_conditions(conditions, sample_table, batch_column)

    # For each LDA matrix, get precision and recall!
    precisions = []
    recalls = []
    
    # Use "range(2)" only when needed to test on a small number of removed components!
    # for ncomps in range(2):
    for ncomps in range(len(components_removed)):
    
        if verbosity >= 1:
            print "Evaluating {} batch effect for matrix with {} components removed".format(batch_column, ncomps)
        
        # Get the matrix with 'ncomps' components removed
        matrix = matrix_3d[ncomps]
        
        if num_test_batches > -1:
            # Create a small matrix and set of batch classes for faster testing!
            small_batches_uniq = np.unique(batches)[0:num_test_batches]
            small_cond_inds = np.array([i for i, batch in enumerate(batches) if batch in small_batches_uniq])
            batches = batches[small_cond_inds]
            matrix = matrix[:, small_cond_inds]

        # Compute the correlation matrix here, since it is used for both the
        # PR curves AND the histograms!
        if verbosity >= 1:
            print "\tComputing correlation matrix"
        corr_matrix = np.corrcoef(matrix, rowvar = 0)
        
        # Compute the precision and recall, based on ranking Pearson
        # correlation coefficients of the profiles and asking if the
        # positive pairs tend to share the same batch
        if verbosity >= 1:
            print "\tCalculating precision and recall"
        precision, recall, threshold = compute_PR_vectors(corr_matrix, batches, verbosity)

        # Dump the PR values out to files for reference
        pr_folder = os.path.join(output_folder, 'PR_analysis')
        if not os.path.isdir(pr_folder):
            os.makedirs(pr_folder)
        print_pr_values(precision, recall, threshold, ncomps, pr_folder)

        # Dump the PR values into lists so I can plot everything together!
        precisions.append(precision)
        recalls.append(recall)

        # Since I'm in the loop, might as well generate batch effect
        # histograms as well...and print the values out to files.
        if verbosity >= 1:
            print "\tCalculating batch correlations"
            print "\t\tWithin batch"
        within_batch_groups, within_batch_corrs, within_batch_n = cf.get_group_correlations(corr_matrix, batches, within_group = True)
        if verbosity >= 1:
            print "\t\tBetween batches"
        between_batch_groups, between_batch_corrs, between_batch_n = cf.get_group_correlations(corr_matrix, batches, within_group = False)
       
        hist_folder = os.path.join(output_folder, 'histograms')
        if not os.path.isdir(hist_folder):
            os.makedirs(hist_folder)
        print_corrs(within_batch_groups, within_batch_corrs, within_batch_n, 'within-batch', ncomps, hist_folder)
        print_corrs(between_batch_groups, between_batch_corrs, between_batch_n, 'between-batch', ncomps, hist_folder)
        if verbosity >= 1:
            print "\t\tPlotting histograms"
        plot_histograms(within_batch_corrs, between_batch_corrs, ncomps, batch_column, hist_folder)


    # Plot the precision and recall vectors!
    plot_pr_curves_diff_scales(components_removed, precisions, recalls, batch_column, pr_folder)


# If this is run from the command line (which is the intention, then run the rest of this script!
if __name__ == '__main__':

    # Get arguments

    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_file', help = 'The dataset on which to perform batch correction.')
    parser.add_argument('sample_table', help = 'The sample table corresponding to the dataset.')
    parser.add_argument('batch_column', help = 'The column from the sample table that defines the batches.')
    parser.add_argument('-v', '--verbosity', help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')
    parser.add_argument('--num_test_batches', help = 'The number of unique batches to use. ONLY for testing purposes. To be used to reduce time and benchmark, not to generate accurate results')

    args = parser.parse_args()
    
    # Take care of verbosity right away
    if args.verbosity is None:
        verbosity = 1
    elif args.verbosity.isdigit():
        verbosity = int(args.verbosity)
    else:
        verbosity = 1

    # Get the data ready to rumble!
    dataset_file = os.path.abspath(args.dataset_file)
    assert os.path.isfile(dataset_file), "Dataset file does not exist."
    dataset = load_3d_dataset(args.dataset_file)

    if verbosity >= 2:
        print dataset

    assert os.path.isfile(args.sample_table), "Sample table file does not exist."
    sample_table = read_sample_table(args.sample_table)

    batch_column = args.batch_column

    input_folder = os.path.dirname(dataset_file)
    output_folder = os.path.join(input_folder, '{}_effect_eval'.format(batch_column))
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    if args.num_test_batches is not None:
        assert args.num_test_batches.isdigit(), "--num_test_batches option must be an integer"
        num_test_batches = int(args.num_test_batches)
    else:
        num_test_batches = -1

    main(dataset, sample_table, batch_column, output_folder, verbosity, num_test_batches)

