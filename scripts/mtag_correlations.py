# This script reads in all of the per-lane condition-strain interaction
# files (z-scores) and computes index tag correlations. It outputs
# graphs and a sorted list with the most correlated index tags at the
# top, so that the "worst" offenders can be removed during the count
# matrix filtering step.

import pandas as pd
import numpy as np
import scipy
import sys, os, gzip
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import itertools as it
import cPickle

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append('./lib')

import config_file_parser as cfp
import compressed_file_opener as cfo
import cg_file_tools as cg_file


def get_sample_table(config_params):

    filename = config_params['sample_table_file']

    # Read everything in as a string, to prevent vexing
    # number interpretation problems! Methods further down
    # can coerce to different types.
    tab = pd.read_table(filename, dtype = 'S')
    return tab

def get_lane_interactions_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'interactions', lane_id)

def get_index_tag_correlation_path(config_params):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'index_tag_correlations')

def get_all_lane_ids(sample_table):

    return np.unique(np.array(sample_table['lane']))

def combine_zscore_matrices(config_params):

    sample_table = get_sample_table(config_params)
    all_lane_ids = get_all_lane_ids(sample_table)
    
    zscore_matrix_list = []
    condition_id_list = []

    for lane_id in all_lane_ids:
        print lane_id
        lane_interactions_path = get_lane_interactions_path(config_params, lane_id)
        lane_interactions_filename = os.path.join(lane_interactions_path, '{}_scaled_dev.dump.gz'.format(lane_id))
        f = gzip.open(lane_interactions_filename)

        gene_barcode_ids, condition_ids, zscore_matrix = cPickle.load(f)
        f.close()
        zscore_matrix_list.append(zscore_matrix)
        condition_id_list.append(condition_ids)

    all_zscore_matrix = np.hstack(zscore_matrix_list)
    all_condition_ids = np.hstack(condition_id_list)

    # Note: for all matrices, the gene_barcode_ids, should be the same. When generating the count matrix,
    # I do not let any gene_barcode_ids disappear, even if they have no counts whatsoever. This should
    # simplify initial processing steps, and the filtering does not need to happen until later.
    return gene_barcode_ids, all_condition_ids, all_zscore_matrix

def dump_dataset(dataset, filename):

    f = gzip.open(filename, 'wb')
    cPickle.dump(dataset, f)
    f.close()

def get_control_condition_ids(dataset, sample_table):

    [barcode_gene_ids, condition_ids, matrix] = dataset

    bool_dict = {'True': True, 'TRUE': True, 'False': False, 'FALSE': False}
    control_bool_ind = np.array([bool_dict[x] for x in sample_table['control?']])
    control_table = sample_table[control_bool_ind]
    control_screen_names = control_table['screen_name']
    control_expt_ids = control_table['expt_id']
    control_condition_ids = ['{0}-{1}'.format(*x) for x in it.izip(control_screen_names, control_expt_ids)]

    control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if cond_id in control_condition_ids])
    final_control_condition_ids = condition_ids[control_condition_indices]

    return final_control_condition_ids

def get_control_dataset(dataset, control_condition_ids):

    [barcode_gene_ids, condition_ids, matrix] = dataset

    control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if cond_id in control_condition_ids])

    control_condition_ids = condition_ids[control_condition_indices]
    control_matrix = matrix[:, control_condition_indices]

    return [barcode_gene_ids, control_condition_ids, control_matrix]

def get_control_index_tag_correlations(control_dataset, sample_table):

    gene_barcode_ids, condition_ids, matrix = control_dataset

    # Must replace all NaNs in the matrix with zeroes for corrcoef to return valid numbers
    matrix = matrix.astype(np.float)
    matrix[np.isnan(matrix)] = 0

    # Remove rows (strains) and columns (conditions) that do not have
    # any values - they kill the clustering process!
    good_rows = np.nansum(matrix, axis = 1).astype(np.bool)
    good_cols = np.nansum(matrix, axis = 0).astype(np.bool)

    gene_barcode_ids = np.array(gene_barcode_ids)[good_rows]
    condition_ids = np.array(condition_ids)[good_cols]
    matrix = matrix[np.ix_(good_rows, good_cols)]

    # rowvar = 0 forces correlation computation on the columns of the matrix
    corr_mat = np.corrcoef(matrix, rowvar = 0)
    print corr_mat

    condition_id_to_index_tag = get_condition_id_to_index_tag(sample_table)
    index_tags = np.unique([condition_id_to_index_tag[condition_id] for condition_id in condition_ids])
    print index_tags[0:10]

    index_tag_to_conditions = {}
    for mapping in condition_id_to_index_tag.iteritems():
        condition_id = mapping[0]
        index_tag = mapping[1]

        if not index_tag_to_conditions.has_key(index_tag):
            index_tag_to_conditions[index_tag] = []
        index_tag_to_conditions[index_tag].append(condition_id)
    print index_tag_to_conditions.items()[0:10]

    # For each index tag, get the correlation matrix that corresponds only to that tag.
    # Then, take the nanmean of the upper triangular of that matrix to get that index tag's
    # correlation with itself across control conditions
    mean_index_tag_corrs = []
    for index_tag in index_tags:
        # Get indices of correlation matrix (symmetric)
        index_tag_condition_ids = index_tag_to_conditions[index_tag]
        print index_tag_condition_ids
        index_tag_inds = np.array([i for i, cond_id in enumerate(condition_ids) if cond_id in index_tag_condition_ids], dtype = np.int)
        print index_tag_inds
        # This bugs out, probably when there is one are zero columns
        # in the correlation matrix with a particular index tag.
        # Solution? Append a big, fat "0.000" to the array of index
        # tag correlations, because we cannot compute the correlation!
        # Also, because of code higher up within this function, no
        # correlations should have nans anymore. Every column and row
        # vector contains no nans and has a nonzero absolute degree.
        if index_tag_inds.size > 1:
            index_tag_ind_1, index_tag_ind_2 = np.hsplit(np.vstack(it.combinations(index_tag_inds, 2)), 2)
            one_index_tag_corr = corr_mat[np.squeeze(index_tag_ind_1), np.squeeze(index_tag_ind_2)]
            # print index_tag_corr_mat
            mean_index_tag_corrs.append(np.nanmean(one_index_tag_corr))
        else:
            # I believe appending a 'nan' is the most appropriate
            # thing to do here, and I do not believe it will screw
            # anything up in future steps that depend on these
            # correlations (not many - just filtering out the most
            # correlated index tags). And np.nan > float(x) returns
            # False, not 'nan'. This is required for the filtering.
            mean_index_tag_corrs.append(np.nan)
       
    mean_index_tag_corrs = np.array(mean_index_tag_corrs)

    corrs_sort_indices = np.argsort(-mean_index_tag_corrs)
    index_tags_sorted = index_tags[corrs_sort_indices]
    mean_index_tag_corrs_sorted = mean_index_tag_corrs[corrs_sort_indices]

    return index_tags_sorted, mean_index_tag_corrs_sorted


def get_condition_id_to_index_tag(sample_table):

    index_tags = list(sample_table['index_tag'])
    screen_names = list(sample_table['screen_name'])
    expt_ids = list(sample_table['expt_id'])

    condition_id_to_index_tag = {}

    for i in range(len(index_tags)):
        index_tag = index_tags[i]
        screen_name = screen_names[i]
        expt_id = expt_ids[i]
       
        condition_id = '{0}-{1}'.format(screen_name, expt_id)    
        condition_id_to_index_tag[condition_id] = index_tag        

    return condition_id_to_index_tag

def write_index_tag_corrs(index_tags, index_tag_corrs, index_tag_path):
    
    txt_filename = os.path.join(index_tag_path, 'control_index_tag_correlations.txt')
    dump_filename = os.path.join(index_tag_path, 'control_index_tag_correlations.dump')

    txt_f = open(txt_filename, 'wt')
    txt_f.write('index_tag\tcorrelation\n')
    for i in range(len(index_tags)):
        txt_f.write('{0}\t{1}\n'.format(index_tags[i], index_tag_corrs[i]))
    txt_f.close()

    dump_f = open(dump_filename, 'wt')
    cPickle.dump([index_tags, index_tag_corrs], dump_f)
    dump_f.close()

def plot_control_index_tag_correlations(index_tag_corrs, index_tag_path):

    plt_filename = os.path.join(index_tag_path, 'control_index_tag_correlations_hist.png')
    
    plt.figure()
    plt.hist(index_tag_corrs, 20, (-1, 1))
    plt.xlabel('Within index tag correlation')
    plt.ylabel('Number of occurrences')
    plt.savefig(plt_filename)


def main(config_file):

    # Read in the config params
    print 'parsing parameters...'
    config_params = cfp.parse(config_file)
    sample_table = get_sample_table(config_params)

    # Read in all of the z-score matrices and combine into one matrix
    dataset = combine_zscore_matrices(config_params)

    # Get directory for index_tag_correlation analysis
    index_tag_path = get_index_tag_correlation_path(config_params)
    if not os.path.isdir(index_tag_path):
        os.makedirs(index_tag_path)

    # Export the initial combined z-score matrix
    per_lane_zscore_dataset_filename = os.path.join(index_tag_path, 'combined_per_lane_zscore_dataset.dump.gz')
    dump_dataset(dataset, per_lane_zscore_dataset_filename)

    # Get just the control dataset, and dump that out too
    control_condition_ids = get_control_condition_ids(dataset, sample_table)
    control_dataset = get_control_dataset(dataset, control_condition_ids)
    per_lane_control_zscore_dataset_filename = os.path.join(index_tag_path, 'combined_per_lane_control_zscore_dataset.dump.gz')
    dump_dataset(dataset, per_lane_control_zscore_dataset_filename)

    # Get the sorted index tag correlations for control conditions
    index_tags_sorted, control_index_tag_correlations_sorted = get_control_index_tag_correlations(control_dataset, sample_table)
    
    # Export the sorted index tag correlations to dump and text files
    write_index_tag_corrs(index_tags_sorted, control_index_tag_correlations_sorted, index_tag_path)   

    # Plot a histogram of the index tag correlations
    plot_control_index_tag_correlations(control_index_tag_correlations_sorted, index_tag_path)

# call: python mtag_correlations.py <config_file>
if '__name__' == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: python mtag_correlations.py <config_file>'
    else:
        config_file = sys.argv[1]
        main(config_file)
