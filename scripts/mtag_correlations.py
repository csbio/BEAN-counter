#!/usr/bin/env python

VERSION='2.6.0'

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
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append('./lib')

import compressed_file_opener as cfo
import cg_file_tools as cg_file
from cg_common_functions import get_verbosity, get_sample_table, get_amplicon_struct_params, get_barcode_table, parse_yaml, bool_dict
from read_type_functions import determine_read_type, get_seq_params
from version_printing import update_version_file

# import pdb

#def get_sample_table(config_params):
#
#    filename = config_params['sample_table_file']
#
#    # Read everything in as a string, to prevent vexing
#    # number interpretation problems! Methods further down
#    # can coerce to different types.
#    tab = pd.read_table(filename, dtype = 'S')
#    return tab

def get_lane_interactions_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'interactions', lane_id)

def get_index_tag_correlation_path(config_params):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'index_tag_qc')

def get_all_lane_ids(sample_table):

    return np.unique(np.array(sample_table['lane']))

def a_is_row_in_b(a, b):

    return np.any(np.all(a == b, axis = 1))

def combine_zscore_matrices(config_params):

    sample_table = get_sample_table(config_params)
    all_lane_ids = get_all_lane_ids(sample_table)
    
    zscore_matrix_list = []
    condition_id_list = []
    gene_barcode_id_list = []
    all_gene_barcode_ids = set()

    for lane_id in all_lane_ids:
        if get_verbosity(config_params) >= 2:
            print lane_id
        lane_interactions_path = get_lane_interactions_path(config_params, lane_id)
        lane_interactions_filename = os.path.join(lane_interactions_path, '{}_scaled_dev.dump.gz'.format(lane_id))
        f = gzip.open(lane_interactions_filename)

        gene_barcode_ids, condition_ids, zscore_matrix = cPickle.load(f)
        # print "number of nonzero barcodes: {}".format(len(gene_barcode_ids))
        f.close()
        zscore_matrix_list.append(zscore_matrix)
        condition_id_list.append(condition_ids)
        gene_barcode_id_list.append(gene_barcode_ids)
        all_gene_barcode_ids.update(gene_barcode_ids) 

    # Because these matrices may not share all of the same rows
    # in the same order, here I must align the rows of the
    # matrices and add missing rows if needed. Using NaNs to
    # fill the rows in some matrices that didn't have any counts
    # because they are handled later in this script.
    all_gene_barcode_ids = list(all_gene_barcode_ids)
    all_gene_barcode_indices = {x:i for i,x in enumerate(all_gene_barcode_ids)}
    aligned_zscore_matrix_list = []
    for i in range(len(zscore_matrix_list)):
        aligned_zscore_matrix_list.append(np.empty((len(all_gene_barcode_ids), len(condition_id_list[i]))))
        aligned_zscore_matrix_list[i][:] = np.nan
        orig_barcodes = gene_barcode_id_list[i]
        for j, barcode in enumerate(orig_barcodes):
            aligned_zscore_matrix_list[i][all_gene_barcode_indices[barcode], :] = zscore_matrix_list[i][j, :]

    all_zscore_matrix = np.hstack(aligned_zscore_matrix_list)
    all_condition_ids = np.vstack(condition_id_list)

    # print all_zscore_matrix.shape
    assert len(all_gene_barcode_ids) == all_zscore_matrix.shape[0], "Number of barcodes: {}; Number of rows in matrix: {}; these should match!".format(len(all_gene_barcode_ids), len(all_zscore_matrix.shape[0]))
    assert len(all_condition_ids) == all_zscore_matrix.shape[1], "Number of conditions: {}; Number of columns in matrix: {}; these should match!".format(len(all_condition_ids), len(all_zscore_matrix.shape[1]))

    return all_gene_barcode_ids, all_condition_ids, all_zscore_matrix

def generate_barcode_specific_template_profiles(gene_barcodes):
    
    profile_ids = ['A', 'C', 'G', 'T']
    profiles = []
    for base in profile_ids:
        profiles.append(np.array([1 if x[0] == base else 0 for x in gene_barcodes]))
    profile_mat = np.vstack(profiles).transpose()
    
    #pdb.set_trace()

    return profile_ids, profile_mat

def compute_max_correlation_barcode_specific_offenders(template_profile_ids, template_profile_mat, condition_ids, matrix, config_params):

    # First, fill NaNs in matrix with very small random noise
    # Get the indices of NaNs in the matrix
    nan_inds = np.isnan(matrix)
    num_nans = np.sum(nan_inds)

    # Replace the NaN values with extremely small noise values
    # (noise has standard deviation of 1 million times less than the data)
    np.random.seed(15263748)
    data_sd = np.nanstd(matrix)
    noise = np.random.randn(num_nans) * data_sd / 1e6
    matrix[nan_inds] = noise

    # Compute the correlation coefficients between the template profiles and
    # the observed profiles.
    corr_mat = pearson_between_two_mats_columns(matrix, template_profile_mat)

    #pdb.set_trace()

    # Get the max value for each condition (row in this correlation matrix)
    max_corrs = np.max(corr_mat, axis = 1)

    # Match the template profiles to the condition profiles that match them the best
    corr_mat_inds = corr_mat == max_corrs[:, None]
    if get_verbosity(config_params) >= 2:
        print corr_mat_inds[0:10]
    template_profile_ids = np.array(template_profile_ids)
    most_correlated_template_profiles = np.array([','.join(template_profile_ids[bool_inds]) for bool_inds in corr_mat_inds])

    # Sort by the max_correlations, descending order
    sort_inds = np.argsort(-max_corrs)
    max_corrs_sorted = max_corrs[sort_inds]
    condition_ids_sorted = condition_ids[sort_inds]
    most_correlated_template_profiles_sorted = most_correlated_template_profiles[sort_inds]

    # These max_corrs should line up with the condition_ids for printing out to file
    return condition_ids_sorted, max_corrs_sorted, most_correlated_template_profiles_sorted

def pearson_between_two_mats_columns(A, B):
    '''
    Computes Pearson correlations between columns of two matrices.
    The column names of A will be the row names of the resulting correlation matrix.
    Adapted from Divakar on:
    http://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays
    '''
    
    # Subtract the column means from each column
    A_mA = A - A.mean(axis = 0)[None, :]
    B_mB = B - B.mean(axis = 0)[None, :]

    # Get the sum of squares of errors for each column
    ssA = (A_mA**2).sum(0)
    ssB = (B_mB**2).sum(0)

    # And compute the correlation coefficient!
    return np.dot(A_mA.T, B_mB) / np.sqrt(np.dot(ssA[:,None],ssB[None]))

def dump_dataset(dataset, filename):

    f = gzip.open(filename, 'wb')
    cPickle.dump(dataset, f)
    f.close()

def get_control_condition_ids(dataset, sample_table):

    [barcode_gene_ids, condition_ids, matrix] = dataset

    control_bool_ind = np.array([bool_dict[x] for x in sample_table['control?']])
    control_table = sample_table[control_bool_ind]
    control_screen_names = control_table['screen_name']
    control_expt_ids = control_table['expt_id']
    control_condition_ids = np.array(list(it.izip(control_screen_names, control_expt_ids)))

    control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, control_condition_ids)])
    final_control_condition_ids = condition_ids[control_condition_indices]

    return final_control_condition_ids

def get_control_dataset(dataset, control_condition_ids):

    [barcode_gene_ids, condition_ids, matrix] = dataset

    control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, control_condition_ids)])

    control_condition_ids = condition_ids[control_condition_indices]
    control_matrix = matrix[:, control_condition_indices]

    return [barcode_gene_ids, control_condition_ids, control_matrix]

def get_control_index_tag_correlations(control_dataset, sample_table, config_params):

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

    condition_id_to_index_tag = get_condition_id_to_index_tag(sample_table)
    index_tags = np.unique([condition_id_to_index_tag[tuple(condition_id)] for condition_id in condition_ids])
    if get_verbosity(config_params) >= 2:
        print corr_mat
        print index_tags[0:10]

    index_tag_to_conditions = {}
    for mapping in condition_id_to_index_tag.iteritems():
        condition_id = mapping[0]
        index_tag = mapping[1]

        # Reversing the mapping now, each key is an index tag, each value is
        # a list of condition id tuples (should be fine)
        if not index_tag_to_conditions.has_key(index_tag):
            index_tag_to_conditions[index_tag] = []
        index_tag_to_conditions[index_tag].append(condition_id)
    if get_verbosity(config_params) >= 2:
        print index_tag_to_conditions.items()[0:10]

    # For each index tag, get the correlation matrix that corresponds only to that tag.
    # Then, take the nanmean of the upper triangular of that matrix to get that index tag's
    # correlation with itself across control conditions
    mean_index_tag_corrs = []
    for index_tag in index_tags:
        # Get indices of correlation matrix (symmetric)
        index_tag_condition_ids = np.vstack(index_tag_to_conditions[index_tag])
        index_tag_inds = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, index_tag_condition_ids)], dtype = np.int)
        if get_verbosity(config_params) >= 3:
            print index_tag_condition_ids
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
            if get_verbosity(config_params) >= 3:
                print index_tag_corr_mat
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
      
        # Here the condition id is a tuple instead of a ndarray - watch out!!!
        condition_id = (screen_name, expt_id)    
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

def write_barcode_specific_template_corrs(condition_ids, template_correlations, template_ids, index_tag_path):
    
    txt_filename = os.path.join(index_tag_path, 'barcode-specific_template_correlations.txt')
    dump_filename = os.path.join(index_tag_path, 'barcode-specific_template_correlations.dump')

    txt_f = open(txt_filename, 'wt')
    txt_f.write('condition_id\tcorrelation with template\tfirst barcode base\n')
    for i in range(len(condition_ids)):
        txt_f.write('{0}\t{1}\t{2}\n'.format(condition_ids[i], template_correlations[i], template_ids[i]))
    txt_f.close()

    dump_f = open(dump_filename, 'wt')
    cPickle.dump([condition_ids, template_correlations, template_ids], dump_f)
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
    config_params = parse_yaml(config_file)
    amplicon_struct_params = get_amplicon_struct_params(config_params)
    sample_table = get_sample_table(config_params)
    barcode_table = get_barcode_table(config_params)

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
    dump_dataset(control_dataset, per_lane_control_zscore_dataset_filename)

    # Get the sorted index tag correlations for control conditions
    index_tags_sorted, control_index_tag_correlations_sorted = get_control_index_tag_correlations(control_dataset, sample_table, config_params)
    
    # Export the sorted index tag correlations to dump and text files
    write_index_tag_corrs(index_tags_sorted, control_index_tag_correlations_sorted, index_tag_path)   

    # Determine what type of barcoding/indexing scheme was used
    read_params = [get_seq_params(amplicon_struct_params, 'read_1'), get_seq_params(amplicon_struct_params, 'read_2')]
    read_type_dict, read_params_clean = determine_read_type(read_params[0], read_params[1])
   
    #pdb.set_trace()

    # If a single barcoding scheme was used, get the correlations of each
    # profile to the barcode-specific template profiles
    if len(read_type_dict['barcode']) == 1:
        single_read = read_type_dict['barcode'][0]
        barcode_column = amplicon_struct_params[single_read]['genetic_barcode']['barcode_file_column']
        barcode_table.set_index('Strain_ID', drop = False, inplace = True)
        barcodes = barcode_table.loc[dataset[0], barcode_column].values
        template_profile_ids, template_profile_mat = generate_barcode_specific_template_profiles(barcodes)
        condition_ids_sorted, barcode_specific_template_correlations_sorted, template_profile_ids_sorted = compute_max_correlation_barcode_specific_offenders(template_profile_ids, template_profile_mat, dataset[1], dataset[2], config_params)

        # Export the sorted correlations of profiles to the barcode-specific template profiles
        write_barcode_specific_template_corrs(condition_ids_sorted, barcode_specific_template_correlations_sorted, template_profile_ids_sorted, index_tag_path)
    else:
        write_barcode_specific_template_corrs([], [], [], index_tag_path)
    
    ## Plot a histogram of the index tag correlations
    plot_control_index_tag_correlations(control_index_tag_correlations_sorted, index_tag_path)

    update_version_file(index_tag_path, VERSION)
    update_version_file(config_params['output_folder'], VERSION)

# call: python mtag_correlations.py <config_file>
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: python mtag_correlations.py <config_file>'
    else:
        config_file = sys.argv[1]
        main(config_file)
