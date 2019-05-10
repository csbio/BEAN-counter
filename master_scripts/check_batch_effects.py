#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.6.1'

# This script takes in matrices, on which LDA has been performed and from which
# LDA components have been removed, and plots PR curves and historgrams
# visualizing the extent of the multiplex tag effect as each successive LDA
# component is removed.
import pandas as pd
import numpy as np
import scipy
from scipy.stats import rankdata
import scipy.stats as stats
import sys, os, gzip
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import itertools as it
import cPickle
import argparse

#from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve, auc

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'lib'))

from pr import precision_recall_curve
from version_printing import update_version_file

import compressed_file_opener as cfo
import cg_file_tools as cg_file
import correlation_functions as cf
from cg_common_functions import read_sample_table, bool_dict
import plotting_tools as pt

#def read_sample_table(tab_filename):
#
#    # Read everything in as a string, to prevent vexing
#    # number interpretation problems! Methods further down
#    # can coerce to different types.
#    tab = pd.read_table(tab_filename, dtype = 'S')
#    return tab

def filter_sample_table(sample_table, final_conditions):

    all_conditions = np.array([tuple(row[['screen_name', 'expt_id']]) for i, row in sample_table.iterrows()])
    rows_to_keep = [i for i, val in enumerate(all_conditions) if a_is_row_in_b(val, final_conditions)]
    return sample_table.iloc[rows_to_keep]

def a_is_row_in_b(a, b):

    return np.any(np.all(a == b, axis = 1))

def filter_dataset_by_conditions(conditions, matrix, conds_to_keep):

    inds_to_keep = np.array([i for i, val in enumerate(conditions) if a_is_row_in_b(val, conds_to_keep)])
    # print inds_to_keep
    filtered_conditions = conditions[inds_to_keep]
    filtered_matrix = matrix[:, inds_to_keep]

    return [filtered_conditions, filtered_matrix]

def get_include_col_conditions(sample_table, include_col):

    # Allow a leading exclamation point to negate a column. Default is false
    negate = False
    if include_col.startswith('!'):
        negate = True
        include_col = include_col[1:]
    
    # Do some checks to make sure the columns are in the sample table
    assert include_col in sample_table.columns, "Specified include column '{}' not in the sample table".format(include_col)
    
    include_vals = [bool_dict[x] for x in sample_table[include_col]]
    if negate:
        include_vals = np.invert(include_vals)
    conditions = np.array(sample_table.loc[include_vals, ['screen_name', 'expt_id']])
    return conditions

def get_lda_pr_curve_filename(folder, ymax):

    return os.path.join(folder, 'PR_batch-effect_{:.1f}.pdf'.format(ymax))

def get_lda_roc_curve_filename(folder):

    return os.path.join(folder, 'ROC_batch-effect.pdf')

def get_batch_effect_histogram_filename(folder, ncomps, total_comp):

    if total_comp > 1:
        return os.path.join(folder, 'Histogram_batch-effect_{}-components-removed.pdf'.format(ncomps))
    else:
        return os.path.join(folder, 'Histogram_batch-effect.pdf'.format(ncomps))

def get_batch_effect_histogram_ttest_filename(folder):

    return os.path.join(folder, 't-stat_vs_components-removed.pdf')

def get_batch_effect_PR_AUC_filename(folder):

    return os.path.join(folder, 'MW_AUC_vs_components-removed.pdf')

def get_batch_effect_PR_pval_filename(folder):

    return os.path.join(folder, 'MW-U_MW-pval_vs_components-removed.txt')

def get_batch_effect_histogram_pval_filename(folder):

    return os.path.join(folder, 't-stats_p-vals_vs_components-removed.txt')

def load_3d_dataset(data_filename):

    f = gzip.open(data_filename, 'rb')
    ncomps, barcodes, conditions, matrix = cPickle.load(f)
    dataset = [np.array(ncomps), np.array(barcodes), np.array(conditions), matrix]
    f.close()
    return dataset

def load_2d_or_3d_dataset(data_filename):

    f = gzip.open(data_filename, 'rb')
    dataset_raw = cPickle.load(f)
    dataset = [np.array(x) if (i+1) < len(dataset_raw) else x for i,x in enumerate(dataset_raw)]
    f.close()
    # If it's a 2D dataset, add a "0 components removed" array as the first element
    # and add another dimension to the data matrix to turn it into a 3d array. This
    # makes it compatible with all remaining code.
    if len(dataset) == 3:
        new_mat_dims = (1, dataset[2].shape[0], dataset[2].shape[1])
        dataset[2] = dataset[2].reshape(new_mat_dims)
        dataset.insert(0, np.array([0]))
    #print len(dataset_raw)
    #print len(dataset)
    return dataset

def get_labels_from_conditions(conditions, sample_table, label_name):

    sample_table = sample_table.set_index(['screen_name', 'expt_id'])
    conditions_tuple = [tuple(x) for x in conditions]
    return np.array(sample_table[label_name].ix[conditions_tuple])

def get_nonreplicating_conditions(sample_table, batch_col, nondup_col_list):

    # Do some checks to make sure the columns are in the sample table
    assert batch_col in sample_table.columns, "Specified batch column '{}' not in the sample table".format(batch_col)
    for col in nondup_col_list:
        assert col in sample_table.columns, "Specified column '{}' to prevent within-batch duplication is not in the sample table".format(col)

    # Create an empty dictionary of dictionaries. The first layer of dictionaries is for the name of the
    # property that should not be duplicated within the same batch (for example, the basic name of the
    # condition, or the chemical ID of the drug, etc). For each of these non-duplicating properties, a
    # dictionary is constructed with one key for each batch. Each of the values for each property --> tag
    # combination is a set containing the actual property values that must not be duplicated. This dictionary
    # is built up and then checked to see if any new property/tag combinations have already occurred.
    nondup_dict = {}
    for col in nondup_col_list:
        nondup_dict[col] = {}
        for tag in set(sample_table[batch_col]):
            nondup_dict[col][tag] = set()

    # Iterate over the rows of the sample table, build up the list of conditions
    # that will be used to get the components for removing batch effects.
    inds = []
    for i, row in sample_table.iterrows():
        # First, if the row has been slated to not be included, then don't include! Move on to the next row
        # print row
        if not bool_dict[row['include?']]:
            continue
        # Accept all control conditions, because they should not show strong correlations with batch 
        # identity unless there is a batch effect
        if bool_dict[row['control?']]:
            inds.append(i)
        else:
            accept = True
            tag = row[batch_col]
            for col in nondup_col_list:
                if row[col] in nondup_dict[col][tag]:
                    # If any property (for example, the condition name) has already occurred for a certain
                    # tag, flag that condition as unacceptable to accept as a nonduplicating condition
                    accept = False
            if accept:
                inds.append(i)
                for col in nondup_col_list:
                    nondup_dict[col][tag].add(row[col])

    nonreplicating_conditions = np.array([row[['screen_name', 'expt_id']] for i, row in sample_table.iterrows() if i in inds])

    return nonreplicating_conditions

def filter_3d_dataset_by_conditions(conditions, matrix_3d, conds_to_keep):

    inds_to_keep = np.array([i for i, val in enumerate(conditions) if a_is_row_in_b(val, conds_to_keep)])
    #print inds_to_keep
    #print len(inds_to_keep)
    filtered_conditions = conditions[inds_to_keep]
    filtered_matrix_3d = matrix_3d[:, :, inds_to_keep]

    return [filtered_conditions, filtered_matrix_3d]

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

    # This calculates the "normalized" AUC, since the recall values have not been multiplied by the
    # number of true positives yet.
    norm_AUPR = auc(recall, precision)

    # Reverse the orders so that they can easily be printed to file
    # with the highest recalls first
    return precision[::-1], recall[::-1] * num_true_positives, thresholds[::-1], norm_AUPR

def compute_ROC_vectors(corr_mat, batches, verbosity):

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
    
    fpr, tpr, thresholds = roc_curve(y_true = batches_match, y_score = corr_vec)
    roc_auc = auc(fpr, tpr)

    # Reverse the orders so that they can easily be printed to file
    # with the highest recalls first
    return fpr, tpr, thresholds, roc_auc

def plot_pr_curves_diff_scales(comps_removed, precision_vectors, recall_vectors, auprc_list, batch_column, total_comp, pr_folder):

    ymax_list = [1, 0.1]
    filenames = [get_lda_pr_curve_filename(pr_folder, y) for y in ymax_list]

    for i, ymax in enumerate(ymax_list):
        plot_one_pr_curve(comps_removed, precision_vectors, recall_vectors, auprc_list, ymax, batch_column, total_comp, filenames[i])

def plot_one_pr_curve(comps_removed, precision_vectors, recall_vectors, auprc_list, ymax, batch_column, total_comp, filename):

    colors_lines = pt.get_style_combos(['lines', 'colors'])
    fig = plt.figure(figsize = (7, 7))
    ax = fig.add_subplot(1,1,1)
    for i, ncomps in enumerate(comps_removed):
        # Since sklearn sets the first precision at recall=0 to 1,
        # I remove this value and any subsequent values at recall=0
        # so that plotting only starts with the precision obtained
        # at the threshold of the first true positive observation
        # (in case multiple true positives have the same associated
        # score).
        recall_vec = recall_vectors[i]
        precision_vec = precision_vectors[i]
        first_nonzero_recall_ind = np.nanmin(np.where(recall_vec > 0)[0])
        first_nonzero_recall_val = recall_vec[first_nonzero_recall_ind]
        last_first_nonzero_recall_ind = np.nanmax(np.where(recall_vec == first_nonzero_recall_val)[0])
        recall_vec = recall_vec[last_first_nonzero_recall_ind:]
        precision_vec = precision_vec[last_first_nonzero_recall_ind:]
        # On to plotting
        style = colors_lines.next()
        ax.plot(recall_vec, precision_vec, style, label = '{} (normalized AUPRC = {:.3f})'.format(ncomps, auprc_list[i]))
    ax.set_xscale('log')
    ax.set_ylim([0, ymax])
    ax.set_ylabel('Precision')
    ax.set_xlabel('Recall')
    plt.title('Precision-recall analysis of\n"{}" batch effects'.format(batch_column))
    if total_comp > 1:
        plt.title('Precision-recall analysis of\n"{}" batch effects'.format(batch_column))
        lgd = ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, title = "Number of components\nremoved")
        plt.savefig(filename, bbox_extra_artists = ([lgd]), bbox_inches = 'tight')
    else:
        plt.title('Precision-recall analysis of "{}"\nbatch effects (normalized AUPRC = {:.3f})'.format(batch_column, auprc_list[i]))
        plt.savefig(filename)

def plot_roc_curves(comps_removed, fpr_vecs, tpr_vecs, auc_list, batch_column, total_comp, pr_folder):

    filename = get_lda_roc_curve_filename(pr_folder)

    colors_lines = pt.get_style_combos(['lines', 'colors'])
    fig = plt.figure(figsize = (7, 7))
    ax = fig.add_subplot(1,1,1)
    for i, ncomps in enumerate(comps_removed):
        style = colors_lines.next()
        ax.plot(fpr_vecs[i], tpr_vecs[i], style, label = '{} (AUC = {:.3f})'.format(ncomps, auc_list[i]))
    ax.set_ylim([0, 1])
    ax.set_xlabel('False Positive Rate (1 - Specificity)')
    ax.set_ylabel('True Positive Rate (Sensitivity)')
    if total_comp > 1:
        plt.title('ROC analysis of\n"{}" batch effects'.format(batch_column))
        lgd = ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, title = "Number of components\nremoved")
        plt.savefig(filename, bbox_extra_artists = ([lgd]), bbox_inches = 'tight')
    else:
        plt.title('ROC analysis of "{}"\nbatch effects (AUC = {:.3f})'.format(batch_column, auc_list[i]))
        plt.savefig(filename)

def plot_histograms(within_batch_corrs, between_batch_corrs, ncomps, batch_column, hist_folder, total_comp, verbosity):

    filename = get_batch_effect_histogram_filename(hist_folder, ncomps, total_comp)

    hist_bins = np.linspace(-1, 1, 41)
    bin_centers = 0.5 * (hist_bins[1:] + hist_bins[:-1])
    # Due to a change in numpy, I must add "range=(bins.min(),bins.max())"
    # https://github.com/numpy/numpy/issues/7503
    y_within_batch, bins = np.histogram(within_batch_corrs, bins = hist_bins, range=(hist_bins.min(), hist_bins.max()), density = True)
    y_between_batch, bins = np.histogram(between_batch_corrs, bins = hist_bins, range=(hist_bins.min(), hist_bins.max()), density = True)

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
    if total_comp > 1:
        plt.title('"{}" batch correlations,\n{} components removed'.format(batch_column, ncomps))
    else:
        plt.title('"{}" batch correlations'.format(batch_column, ncomps))
    lgd = ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, title = "Correlation origin")
    plt.savefig(filename, bbox_extra_artists = ([lgd]), bbox_inches = 'tight')

    # Here I compute the t-statistic for each distribution
    # where the null hypothesis is that the mean is 0.
    # Also, take the sample standard deviation, not the population
    # standard deviation (ddof = 1).
    #if verbosity >= 2:
    #    print 'within batch t-statistic calculation:'
    #    print '{} / ({} / sqrt({}))'.format(mean_within_batch_cor, np.nanstd(within_batch_corrs, ddof = 1), np.size(within_batch_corrs))
    #    print ''
    #    print 'between batch t-statistic calculation:'
    #    print '{} / ({} / sqrt({}))'.format(mean_between_batch_cor, np.nanstd(between_batch_corrs, ddof = 1), np.size(between_batch_corrs))

    #within_batch_t = mean_within_batch_cor / (np.nanstd(within_batch_corrs, ddof = 1) / np.sqrt(np.size(within_batch_corrs)))
    #between_batch_t = mean_between_batch_cor / (np.nanstd(between_batch_corrs, ddof = 1) / np.sqrt(np.size(between_batch_corrs)))

    ## Take the t-stats and convert to pvalues.
    #within_batch_p = stats.t.sf(np.abs(within_batch_t), np.size(within_batch_corrs) - 1) * 2
    #between_batch_p = stats.t.sf(np.abs(between_batch_t), np.size(between_batch_corrs) - 1) * 2

    ## Return t-statstics and p-values
    #return [np.array([[within_batch_t, between_batch_t]]), np.array([[within_batch_p, between_batch_p]])]

    # Compute a Mann-Whitney U statistic and p-value comparing the distributions
    # of within-batch and between-batch correlations (Mann-Whitney U is equivalent
    # to the AUC for a ROC curve comparing the two classes.
    within_batch_corrs_nanfree = within_batch_corrs[np.invert(np.isnan(within_batch_corrs))]
    between_batch_corrs_nanfree = between_batch_corrs[np.invert(np.isnan(between_batch_corrs))]
    MW_U, MW_p = stats.mannwhitneyu(within_batch_corrs_nanfree, between_batch_corrs_nanfree)
    if verbosity >= 2 :
        print within_batch_corrs_nanfree
        print between_batch_corrs_nanfree
    MW_AUC = MW_U / (np.size(within_batch_corrs_nanfree) * np.size(between_batch_corrs_nanfree))
    return MW_AUC, MW_p

def plot_MW_AUCs(comps_removed, AUC_vals, batch_column, pr_folder):

    filename = get_batch_effect_PR_AUC_filename(pr_folder)

    fig = plt.figure(figsize = (7, 7))
    ax = fig.add_subplot(1,1,1)
    ax.plot(comps_removed, AUC_vals, 'k-')
    ax.set_xlabel('Number of components removed')
    ax.set_ylabel(r'Mann-Whitney $U$' '\n' r'(Area under ROC curve)')
    plt.title('Comparison of within- and between-batch \ncorrelation distributions (area under ROC curve)'.format(batch_column))
    plt.savefig(filename)

def plot_t_stats(comps_removed, t_stats, batch_column, hist_folder):

    # This function plots the change in t-statistic for both
    # within- and between-batch correlations, as components
    # (either LDA or SVD) are removed.
    filename = get_batch_effect_histogram_ttest_filename(hist_folder) 
   
    colors = pt.color_cycle()
    fig = plt.figure(figsize = (7, 7))
    ax = fig.add_subplot(1,1,1)
    color = colors.next()
    ax.plot(comps_removed, t_stats[:, 0], color + '-', label = 'Within batch')
    color = colors.next()
    ax.plot(comps_removed, t_stats[:, 1], color + '-', label = 'Between batch')
    ax.set_xlabel('Number of components removed')
    ax.set_ylabel(r'$t$-statistic' '\n' r'($H_0: \bar \mu = 0$)')
    lgd = ax.legend(bbox_to_anchor = (1.05, 1), loc = 2, title = 'Correlation origin')
    plt.title(r'$t$-statistics of {}' '\nbatch correlation distributions'.format(batch_column))
    plt.savefig(filename, bbox_extra_artists = ([lgd]), bbox_inches = 'tight')

def write_MW_pvals(comps_removed, AUC_vals, pvals, batch_column, pr_folder):

    filename = get_batch_effect_PR_pval_filename(pr_folder)

    f = open(filename, 'wt')
    f.write('# components removed\tMann-Whitney U statistic\tMann-Whitney p-value\n')
    for i, comp in enumerate(comps_removed):
        f.write('{}\t{}\t{}\n'.format(comps_removed[i], AUC_vals[i], pvals[i]))

    f.close()

def write_p_vals(comps_removed, t_stats, p_vals, batch_column, hist_folder):

    filename = get_batch_effect_histogram_pval_filename(hist_folder)

    f = open(filename, 'wt')
    f.write('# components removed\twithin-batch t-stat\twithin-batch p-value\tbetween-batch t-stat\tbetween-batch p-value\n')
    for i, comp in enumerate(comps_removed):
        f.write('{}\t{}\t{}\t{}\t{}\n'.format(comps_removed[i], t_stats[i, 0], p_vals[i, 0], t_stats[i, 1], p_vals[i, 1]))

    f.close()

def print_pr_values(precision, recall, threshold, ncomps, total_comp, pr_folder):

    pr_data_folder = os.path.join(pr_folder, 'PR_data')
    if not os.path.isdir(pr_data_folder):
        os.makedirs(pr_data_folder)
    if total_comp > 1:
        filename = os.path.join(pr_data_folder, 'PR_{}-components-removed.txt.gz'.format(ncomps))
    else:
        filename = os.path.join(pr_data_folder, 'PR.txt.gz')

    # Make sure sure threshold vector is same length as recall and precision vectors
    threshold = np.append(threshold, np.nan)

    f = gzip.open(filename, 'wb')
    f.write('Recall\tPrecision\tThreshold\n')
    for i, rec in enumerate(recall):
        f.write('{}\t{}\t{}\n'.format(rec, precision[i], threshold[i]))
    f.close()

def print_roc_values(fpr, tpr, threshold, ncomps, total_comp, pr_folder):

    roc_data_folder = os.path.join(pr_folder, 'ROC_data')
    if not os.path.isdir(roc_data_folder):
        os.makedirs(roc_data_folder)
    if total_comp > 1:
        filename = os.path.join(roc_data_folder, 'ROC_{}-components-removed.txt.gz'.format(ncomps))
    else:
        filename = os.path.join(roc_data_folder, 'ROC.txt.gz')

    f = gzip.open(filename, 'wb')
    f.write('TPR\tFPR\tThreshold\n')
    for i, tp in enumerate(tpr):
        f.write('{}\t{}\t{}\n'.format(tp, fpr[i], threshold[i]))
    f.close()

def print_corrs(groups, corrs, n, string, ncomps, total_comp, hist_folder):

    corr_folder = os.path.join(hist_folder, 'correlation_data')
    if not os.path.isdir(corr_folder):
        os.makedirs(corr_folder)
    if total_comp > 1:
        filename = os.path.join(corr_folder, '{}_correlations_{}-components-removed.txt.gz'.format(string, ncomps))
    else:
        filename = os.path.join(corr_folder, '{}_correlations.txt.gz'.format(string))

    # Sort the correlations and the groups in descending order!
    # Sorting using the negative of corrs forces NaNs to stay
    # at the end of the vector despite the otherwise reversed
    # sorting order (descending)
    sort_inds = np.argsort(-corrs)
    corrs_sorted = corrs[sort_inds]
    groups_sorted = [tuple(x) for x in np.array(groups)[sort_inds]]
    n_sorted = n[sort_inds]

    f = gzip.open(filename, 'wt')
    f.write('batch_1\tbatch_2\tavg_correlation\tnum_non-NaN-correlations\n')
    
    for i, group in enumerate(groups_sorted):
        f.write('{}\t{}\t{}\t{}\n'.format(group[0], group[1], corrs_sorted[i], n_sorted[i]))

    f.close()
    return None

def main(dataset_3d, sample_table, batch_column, nondup_col_list, include_column, output_folder, verbosity, num_test_batches):

    # Load in the entire 3-d stacked matrix and dimension names!
    components_removed, barcodes, conditions, matrix_3d = dataset
    total_comp = len(components_removed)

    if verbosity >= 2:
        print 'matrix dimensions should be: {}, {}, {}'.format(len(components_removed), len(barcodes), len(conditions))
        print 'matrix dimensions are: {}'.format(matrix_3d.shape)

    # If an include column was specified:
    # include only conditions that were specified in the optional "include_column" argument
    if include_column is not None:
        batch_viz_include_cond = get_include_col_conditions(sample_table, include_column)
        if verbosity >= 2:
            print batch_viz_include_cond[0:10]
            print conditions[0:10]
        sample_table = filter_sample_table(sample_table, batch_viz_include_cond)
        conditions, matrix_3d = filter_3d_dataset_by_conditions(conditions, matrix_3d, batch_viz_include_cond)
    
    if verbosity >= 2:
        print 'matrix dimensions should be: {}, {}, {}'.format(len(components_removed), len(barcodes), len(conditions))
        print 'matrix dimensions are: {}'.format(matrix_3d.shape)

    # If columns were specified for values that should not be duplicated in the same batch...
    if nondup_col_list != ['none']:
        # Get the indices of the nonreplicating conditions
        nonrep_cond_inds = get_nonreplicating_conditions(sample_table, batch_column, nondup_col_list)

        # Filter the 3d matrix such that it only contains the nonreplicating conditions
        conditions, matrix_3d = filter_3d_dataset_by_conditions(conditions, matrix_3d, nonrep_cond_inds)
        if verbosity >= 2:
            print conditions[0:10]

    if verbosity >= 2:
        print 'matrix dimensions should be: {}, {}, {}'.format(len(components_removed), len(barcodes), len(conditions))
        print 'matrix dimensions are: {}'.format(matrix_3d.shape)
    
    # Filter the sample table, so that no samples that have been eliminated
    # are used in further processing
    sample_table = filter_sample_table(sample_table, conditions)
   
    if verbosity >= 2:
        print sample_table[0:10]
        print sample_table.shape

    # Get a vector of the batches that lines up with the vectors of screen_name, expt_id, name, etc
    batches = get_labels_from_conditions(conditions, sample_table, batch_column)
    if verbosity >= 2:
        print len(conditions)
        print len(batches)
        print matrix_3d.shape

    # For each LDA matrix, get precision and recall, as well as tpr/fpr!
    precisions = []
    recalls = []
    auprs = []
    tprs = []
    fprs = []
    aucs = []

    ## Initialize empty data structures to hold the t-statistics and p-values
    #t_stats = np.zeros([0, 2])
    #p_vals = np.zeros([0, 2])
    AUC_values = []
    MW_pvals = []
    
    # Set up the output PR and ROC results folder
    pr_folder = os.path.join(output_folder, 'PR_ROC_analysis')
    if not os.path.isdir(pr_folder):
        os.makedirs(pr_folder)
    
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
        # Here I turn NaNs into zeros so profiles with some NaNs don't break everything
        matrix[np.isnan(matrix)] = 0.
        corr_matrix = np.corrcoef(matrix, rowvar = 0)
        if verbosity >= 3:
            print "Number of NaNs in correlation matrix: {}".format(np.sum(np.isnan(corr_matrix)))
        
        # Compute the precision and recall, based on ranking Pearson
        # correlation coefficients of the profiles and asking if the
        # positive pairs tend to share the same batch
        if verbosity >= 1:
            print "\tCalculating precision and recall"
        precision, recall, threshold_pr, aupr = compute_PR_vectors(corr_matrix, batches, verbosity)
        # Since there are still some bugs appearing here, add in
        # a verbose option.
        if verbosity >= 2:
            print "\t\tshape of recall vector: {}".format(recall.shape)
            print "\t\tshape of precision vector: {}".format(precision.shape)
            print "\t\tshape of threshold vector: {}".format(threshold_pr.shape)

        # Generate ROC curves based on the same values used to generate
        # the PR curves. I am particularly interested in the AUC values
        # from this analysis in determining when to stop removing
        # components
        if verbosity >= 1:
            print '\tCalculating true and false positive rates for ROC'
        fpr, tpr, threshold_roc, auc = compute_ROC_vectors(corr_matrix, batches, verbosity)
        if verbosity >= 2:
            print "\t\tshape of FPR vector: {}".format(fpr.shape)
            print "\t\tshape of TPR vector: {}".format(tpr.shape)
            print "\t\tshape of threshold vector: {}".format(threshold_roc.shape)

        # Dump the PR and ROC values out to files for reference
        print_pr_values(precision, recall, threshold_pr, ncomps, total_comp, pr_folder)
        print_roc_values(fpr, tpr, threshold_roc, ncomps, total_comp, pr_folder)

        # Dump the PR values into lists so I can plot everything together!
        precisions.append(precision)
        recalls.append(recall)
        auprs.append(aupr)
        
        # Dump the ROC values into lists so I can plot everything together!
        fprs.append(fpr)
        tprs.append(tpr)
        aucs.append(auc)

        # Since I'm in the loop, might as well generate batch effect
        # histograms as well...and print the values out to files.
        if verbosity >= 1:
            print "\tCalculating batch correlations"
            print "\t\tWithin batch"
        within_batch_groups, within_batch_corrs, within_batch_n = cf.get_group_correlations(corr_matrix, batches, within_group = True, verbosity = verbosity)
        if verbosity >= 1:
            print "\t\tBetween batches"
        between_batch_groups, between_batch_corrs, between_batch_n = cf.get_group_correlations(corr_matrix, batches, within_group = False, verbosity = verbosity)
       
        hist_folder = os.path.join(output_folder, 'histograms')
        if not os.path.isdir(hist_folder):
            os.makedirs(hist_folder)
        print_corrs(within_batch_groups, within_batch_corrs, within_batch_n, 'within-batch', ncomps, total_comp, hist_folder)
        print_corrs(between_batch_groups, between_batch_corrs, between_batch_n, 'between-batch', ncomps, total_comp, hist_folder)
        if verbosity >= 1:
            print "\t\tPlotting histograms"
        #between_within_t_stats, between_within_p_vals = plot_histograms(within_batch_corrs, between_batch_corrs, ncomps, batch_column, hist_folder, verbosity)
        plot_histograms(within_batch_corrs, between_batch_corrs, ncomps, batch_column, hist_folder, total_comp, verbosity)
        #AUC_values.append(MW_AUC_val)
        #MW_pvals.append(MW_pval)
        #t_stats = np.vstack([t_stats, between_within_t_stats])
        #p_vals = np.vstack([p_vals, between_within_p_vals])
    #AUC_values = np.array(AUC_values)
    #MW_pvals = np.array(MW_pvals)

    # Plot the precision and recall vectors!
    plot_pr_curves_diff_scales(components_removed, precisions, recalls, auprs, batch_column, total_comp, pr_folder)
    plot_roc_curves(components_removed, fprs, tprs, aucs, batch_column, total_comp, pr_folder)

    # Plot t-statistics for the histograms (and write out associated p-values)
    # as they evolve over the course of component removal
    #plot_MW_AUCs(components_removed, AUC_values, batch_column, pr_folder)
    #write_MW_pvals(components_removed, AUC_values, MW_pvals, batch_column, pr_folder)
    update_version_file(output_folder, VERSION)

# If this is run from the command line (which is the intention, then run the rest of this script!
if __name__ == '__main__':

    # Get arguments

    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_file', help = 'The dataset on which to evaluate batch effects. Can be either a 2D matrix dataset or a 3D stacked matrix dataset (the latter being the immediate output of batch correction).')
    parser.add_argument('sample_table', help = 'The sample table corresponding to the dataset.')
    parser.add_argument('batch_column', help = 'The column from the sample table that defines the batches.')
    parser.add_argument('nondup_columns', help = 'Comma delimited. The columns that contain condition attributes that should not be duplicated within the same batch.')
    parser.add_argument('-incl', '--include_column', help = 'Column in the sample table that specifies True/False whether or not the condition in that row should be included in the batch effect evaluation/visualization.')
    parser.add_argument('-v', '--verbosity', type = int, default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')
    parser.add_argument('--num_test_batches', type = int, default = -1, help = 'The number of unique batches to use. ONLY for testing purposes. To be used to reduce time and benchmark, not to generate accurate results')

    args = parser.parse_args()
    
    # Get the data ready to rumble!
    dataset_file = os.path.abspath(args.dataset_file)
    assert os.path.isfile(dataset_file), "Dataset file does not exist."
    dataset = load_2d_or_3d_dataset(args.dataset_file)

    if args.verbosity >= 2:
        print dataset

    assert os.path.isfile(args.sample_table), "Sample table file does not exist."
    sample_table = read_sample_table(args.sample_table)

    batch_column = args.batch_column

    # If the user specifies "none", then no columns are used to prevent duplicates of
    # certain values within a particular batch
    if args.nondup_columns is 'none':
        nondup_col_list = ['none']
    else:
        nondup_col_list = args.nondup_columns.split(',')

    # Take care of the include column optional argument
    if args.include_column is not None:
        include_column = args.include_column
    else:
        include_column = None

    input_folder = os.path.dirname(dataset_file)
    output_folder = os.path.join(input_folder, '{}_effect_eval'.format(batch_column))
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    main(dataset, sample_table, batch_column, nondup_col_list, include_column, output_folder, args.verbosity, args.num_test_batches)

