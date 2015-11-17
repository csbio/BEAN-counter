#!/usr/bin/env python
# This script reads in the combined count matrix and filters out all of the
# conditions and strains that either do not meet specified quality measures
# and/or are specifically slated for removal

import pandas as pd
import numpy as np
import scipy
import sys, os, gzip
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import cPickle
import gzip
import itertools as it

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append(os.path.join(barseq_path, 'lib'))

import config_file_parser as cfp
import compressed_file_opener as cfo
import cg_file_tools as cg_file
from cg_common_functions import *

def get_sample_table(config_params):

    filename = config_params['sample_table_file']

    # Read everything in as a string, to prevent vexing
    # number interpretation problems! Methods further down
    # can coerce to different types.
    tab = pd.read_table(filename, dtype = 'S')
    return tab

def get_lane_data_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'intermediate', lane_id)

def get_filtered_strains_conditions_path(config_params):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'filtered_strains_conditions')

def get_dumped_count_matrix_filename(config_params, lane_id):

    lane_folder = get_lane_data_path(config_params, lane_id)
    return os.path.join(lane_folder, '{}_barseq_matrix.dump.gz'.format(lane_id))

def load_dumped_count_matrix(config_params, lane_id):

    filename = get_dumped_count_matrix_filename(config_params, lane_id)
    f = gzip.open(filename, 'rb')
    barcodes, conditions, matrix = cPickle.load(f)
    dataset = [np.array(barcodes), np.array(conditions), matrix]
    f.close()

    return dataset

def dump_dataset(dataset, filename):

    f = gzip.open(filename, 'wb')
    cPickle.dump(dataset, f)
    f.close()

def a_is_row_in_b(a, b):

    return np.any(np.all(a == b, axis = 1))

def filter_dataset_for_include_2(dataset, sample_table):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
   
    bool_dict = {'True': True, 'TRUE': True, 'False': False, 'FALSE': False}
    include_bool_ind = np.array([bool_dict[x] for x in sample_table['include?']])
    include_table = sample_table[include_bool_ind]
    not_include_table = sample_table[np.invert(include_bool_ind)]

    include_screen_names = include_table['screen_name']
    include_expt_ids = include_table['expt_id']
    include_condition_ids = np.array(list(it.izip(include_screen_names, include_expt_ids)))
    
    include_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, include_condition_ids)])
    
    filtered_condition_ids = condition_ids[include_condition_indices]
    filtered_matrix = matrix[:, include_condition_indices]

    # ALSO return the sample table of excluded conditions!
    return [barcode_gene_ids, filtered_condition_ids, filtered_matrix], not_include_table

def filter_dataset_for_barcodes(dataset, config_params):

    barcodes_to_remove_file = config_params['barcodes_to_remove_file']
    if barcodes_to_remove_file.strip() == '':
        return dataset, []
    if not os.path.exists(barcodes_to_remove_file):
        return dataset, []
    with open(barcodes_to_remove_file, 'rt') as f:
        # Throw away the header
        f.readline()
        gene_barcode_ids_to_remove = np.array([line.rstrip().split('\t')[::-1] for line in f])
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
    if get_verbosity(config_params) >= 2:
        print 'barcodes 1-10:'
        print barcode_gene_ids[0:10]
        print '\n\n'
        print 'barcodes to remove:'
        print gene_barcode_ids_to_remove
    inds_to_keep = np.array([i for i, bg_id in enumerate(barcode_gene_ids) if not a_is_row_in_b(bg_id, gene_barcode_ids_to_remove)])

    filtered_barcode_gene_ids = barcode_gene_ids[inds_to_keep]
    filtered_matrix = matrix[inds_to_keep, :]
   
    return [filtered_barcode_gene_ids, condition_ids, filtered_matrix], gene_barcode_ids_to_remove

def get_index_tag_correlation_path(config_params):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'index_tag_qc')

def filter_dataset_for_index_tags(dataset, config_params):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
    
    sample_table = get_sample_table(config_params)

    index_tag_corr_path = get_index_tag_correlation_path(config_params)
    index_tag_corr_filename = os.path.join(index_tag_corr_path, 'control_index_tag_correlations.dump')
   
    index_tag_correlation_cutoff = float(config_params['index_tag_correlation_cutoff'])

    f = open(index_tag_corr_filename, 'rt')
    index_tags, correlations = cPickle.load(f)
    index_tags_to_remove = np.array([index_tags[i] for i,corr in enumerate(correlations) if corr >= index_tag_correlation_cutoff])
    # print 'index_tags_to_remove:'
    # print '\n'.join(index_tags_to_remove) + '\n'
   
    to_remove_idx = sample_table['index_tag'].isin(index_tags_to_remove)
    to_remove_table = sample_table[to_remove_idx]
    to_keep_table = sample_table[~to_remove_idx]

    to_keep_screen_names = to_keep_table['screen_name']
    to_keep_expt_ids = to_keep_table['expt_id']
    to_keep_condition_ids = np.array(list(it.izip(to_keep_screen_names, to_keep_expt_ids)))

    to_keep_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, to_keep_condition_ids)])
    
    filtered_condition_ids = condition_ids[to_keep_condition_indices]
    filtered_matrix = matrix[:, to_keep_condition_indices]

    # ALSO return the sample table of filtered out conditions!
    to_remove_table = to_remove_table.reset_index(drop = True)
    return [barcode_gene_ids, filtered_condition_ids, filtered_matrix], to_remove_table

def filter_dataset_for_barcode_specific_patterns(dataset, config_params):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
    
    sample_table = get_sample_table(config_params)

    index_tag_corr_path = get_index_tag_correlation_path(config_params)
    barcode_spec_corr_filename = os.path.join(index_tag_corr_path, 'barcode-specific_template_correlations.dump')
   
    barcode_specific_correlation_cutoff = float(config_params['barcode_specific_template_correlation_cutoff'])

    f = open(barcode_spec_corr_filename, 'rt')
    barcode_spec_condition_ids, barcode_spec_correlations, start_base = cPickle.load(f)
    cond_ids_to_remove = np.array([barcode_spec_condition_ids[i] for i,corr in enumerate(barcode_spec_correlations) if corr >= barcode_specific_correlation_cutoff])
    if get_verbosity(config_params) >= 2:
        print 'condition_ids_to_remove:'
        print cond_ids_to_remove
   
    to_remove_idx = np.array([a_is_row_in_b(tuple(x[1]), cond_ids_to_remove) for x in sample_table[['screen_name', 'expt_id']].iterrows()])
    to_remove_table = sample_table[to_remove_idx]
    to_keep_table = sample_table[~to_remove_idx]

    to_keep_screen_names = to_keep_table['screen_name']
    to_keep_expt_ids = to_keep_table['expt_id']
    to_keep_condition_ids = np.array(list(it.izip(to_keep_screen_names, to_keep_expt_ids)))

    to_keep_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, to_keep_condition_ids)])
    
    filtered_condition_ids = condition_ids[to_keep_condition_indices]
    filtered_matrix = matrix[:, to_keep_condition_indices]

    ### ALSO return the sample table of filtered out conditions!
    ### But first, add the correlation for each condition to the sample table
    ### and which of the barcode first bases ('A', 'C', 'G', 'T') it matched to
    to_remove_table = to_remove_table.reset_index(drop = True)
    # First, create a table of the correlation information
    cor_tab = pd.DataFrame({'cond_ids': [tuple(x) for x in barcode_spec_condition_ids], 'corrs': barcode_spec_correlations, 'bases': start_base})
    cor_tab = cor_tab.set_index('cond_ids')
    # Now rearrange the table based on the "to_remove_table"
    cor_tab_cond_ids_idx = [tuple(x[1]) for x in to_remove_table[['screen_name', 'expt_id']].iterrows()]
    cor_tab_ordered = cor_tab.ix[cor_tab_cond_ids_idx]
    if get_verbosity(config_params) >= 3:
        print cor_tab_cond_ids_idx
        print cor_tab_ordered
    to_remove_table['barcode-specific_effect_correlation'] = np.array(cor_tab_ordered['corrs'])
    to_remove_table['barcode-specific_effect_start_base'] = np.array(cor_tab_ordered['bases'])

    # And, finally return everything!
    return [barcode_gene_ids, filtered_condition_ids, filtered_matrix], to_remove_table

def filter_dataset_for_count_degree(dataset, config_params, sample_table):

    strain_pass_read_count = float(config_params['strain_pass_read_count'])
    strain_pass_fraction = float(config_params['strain_pass_fraction'])
    condition_pass_read_count = float(config_params['condition_pass_read_count'])
    condition_pass_fraction = float(config_params['condition_pass_fraction'])

    [barcode_gene_ids, condition_ids, matrix] = dataset

    # Get the indices of the strains to keep (and those to eliminate)
    passing_strain_read_count_matrix = matrix >= strain_pass_read_count
    passing_strain_read_count_sums = np.nansum(passing_strain_read_count_matrix, axis = 1)
    passing_strain_fraction = passing_strain_read_count_sums / float(passing_strain_read_count_matrix.shape[1])
    # print passing_strain_fraction
    strains_to_keep_inds = passing_strain_fraction >= strain_pass_fraction
    strains_to_remove_inds = np.invert(strains_to_keep_inds)
    strains_to_remove = barcode_gene_ids[strains_to_remove_inds]

    # Get the indices of the conditions to keep (and those to eliminate)
    passing_condition_read_count_matrix = matrix >= condition_pass_read_count
    passing_condition_read_count_sums = np.nansum(passing_condition_read_count_matrix, axis = 0)
    passing_condition_fraction = passing_condition_read_count_sums / float(passing_condition_read_count_matrix.shape[0])
    if get_verbosity(config_params) >= 3:
        print passing_condition_fraction
    conditions_to_keep_inds = passing_condition_fraction >= condition_pass_fraction
    conditions_to_remove_inds = np.invert(conditions_to_keep_inds)
    conditions_to_remove = condition_ids[conditions_to_remove_inds]
    conditions_to_remove_tuple = [tuple(x) for x in conditions_to_remove]

    sample_table = sample_table.set_index(['screen_name', 'expt_id'])
    conditions_to_remove_table = sample_table.loc[conditions_to_remove_tuple]

    # Return conditions and strains rejected by the count degree criteria 
    filtered_barcode_gene_ids = barcode_gene_ids[strains_to_keep_inds]
    filtered_condition_ids = condition_ids[conditions_to_keep_inds]
    filtered_matrix = matrix[strains_to_keep_inds, :]
    filtered_matrix = filtered_matrix[:, conditions_to_keep_inds]

    return [filtered_barcode_gene_ids, filtered_condition_ids, filtered_matrix], conditions_to_remove_table.reset_index(), strains_to_remove

def create_filtered_output_folder(config_params):
    filtered_path = get_filtered_strains_conditions_path(config_params)
    if not os.path.isdir(filtered_path):
        os.makedirs(filtered_path)

def write_filtered_include_table(tab, config_params):
    filtered_path = get_filtered_strains_conditions_path(config_params)
    # Write out conditions that were specified as "do not include"
    tab_file = os.path.join(filtered_path, 'prespecified_excluded_conditions.txt')
    tab.to_csv(tab_file, sep = '\t', index = False)

def write_count_degree_excluded_conditions(tab, config_params):
    filtered_path = get_filtered_strains_conditions_path(config_params)
    # Write out the conditions excluded because of their count degree
    cond_degree_file = os.path.join(filtered_path, 'count_degree_excluded_conditions.txt')
    tab.to_csv(cond_degree_file, sep = '\t', index = False)

def write_correlated_index_tags_excluded_conditions(tab, config_params):
    filtered_path = get_filtered_strains_conditions_path(config_params)
    # Write out the conditions excluded because they were tagged with index tags
    # that had within-tag correlations above the specified threshold.
    filtered_index_tag_file = os.path.join(filtered_path, 'index_tag_correlation_excluded_conditions.txt')
    tab.to_csv(filtered_index_tag_file, sep = '\t', index = False)

def write_barcode_specific_excluded_conditions(tab, config_params):
    filtered_path = get_filtered_strains_conditions_path(config_params)
    # Write out the conditions excluded because they correlated with any of the
    # 'A', 'C', 'G', or 'T'-specific profiles above the specified threshold.
    filtered_barcode_specific_file = os.path.join(filtered_path, 'barcode-specific_correlation_excluded_conditions.txt')
    tab.to_csv(filtered_barcode_specific_file, sep = '\t', index = False)

def write_filtered_strain_file(barcodes, config_params):
    filtered_path = get_filtered_strains_conditions_path(config_params)
    # Write out the barcodes that were specifed as "do not include"
    filtered_barcode_file = os.path.join(filtered_path, 'prespecified_excluded_strains.txt')
    with open(filtered_barcode_file, 'wt') as f:
        f.write('Barcode\tStrain ID\n')
        for barcode in barcodes:
            f.write('\t'.join(barcode[1::-1]) + '\n')

def write_count_degree_excluded_strains(barcodes, config_params):
    filtered_path = get_filtered_strains_conditions_path(config_params)
    # Write out the barcodes excluded because of their count degree
    ###### Perhaps change this in the future so it exports in the same format as the barcode table
    strain_degree_file = os.path.join(filtered_path, 'count_degree_excluded_strains.txt')
    with open(strain_degree_file, 'wt') as f:
        f.write('Barcode\tStrain ID\n')
        for barcode in barcodes:
            f.write('\t'.join(barcode[1::-1]) + '\n')

def dump_filtered_count_matrix(config_params, dataset):
   
    lane_id = 'all_lanes_filtered'
    dataset_path = get_lane_data_path(config_params, lane_id)
    if not os.path.isdir(dataset_path):
        os.makedirs(dataset_path)

    dataset_filename = get_dumped_count_matrix_filename(config_params, lane_id)
    dump_dataset(dataset, dataset_filename)

def main(config_file):

    # Read in the config params
    config_params = cfp.parse(config_file)
    sample_table = get_sample_table(config_params)

    # Create the folder where all the filtered condition/strain info is dumped
    create_filtered_output_folder(config_params)

    # Get parameters that specify if some steps should be run or not
    bool_dict = {'True': True, 'TRUE': True, 'False': False, 'FALSE': False}
    remove_mtag_offenders = bool_dict[config_params['remove_correlated_index_tags']]
    remove_barcode_specific_conds = bool_dict[config_params['remove_barcode_specific_conditions']]

    dataset = load_dumped_count_matrix(config_params, 'all_lanes')
    if get_verbosity(config_params) >= 2:
        print dataset[2].shape

    if get_verbosity(config_params) >= 1:
        print 'Filtering out...'
        print '\tPrespecified conditions to exclude'
    dataset, filtered_include_tab = filter_dataset_for_include_2(dataset, sample_table)
    if get_verbosity(config_params) >= 2:
        print dataset[2].shape
    write_filtered_include_table(filtered_include_tab, config_params)

    if get_verbosity(config_params) >= 1:
        print '\tPrespecified barcodes to exclude'
    dataset, filtered_barcodes = filter_dataset_for_barcodes(dataset, config_params)
    if get_verbosity(config_params) >= 2:
        print dataset[2].shape
    
    write_filtered_strain_file(filtered_barcodes, config_params)
    
    if remove_mtag_offenders:
        if get_verbosity(config_params) >= 1:
            print '\tHighly-correlated index tags'
        dataset, filtered_index_tag_condition_table = filter_dataset_for_index_tags(dataset, config_params)
        if get_verbosity(config_params) >= 2:
            print dataset[2].shape
        write_correlated_index_tags_excluded_conditions(filtered_index_tag_condition_table, config_params)

    if remove_barcode_specific_conds:
        if get_verbosity(config_params) >= 1:
            print '\tConditions with barcode-specific signatures'
        dataset, filtered_barcode_specific_condition_table = filter_dataset_for_barcode_specific_patterns(dataset, config_params)
        if get_verbosity(config_params) >= 2:
            print dataset[2].shape
        write_barcode_specific_excluded_conditions(filtered_barcode_specific_condition_table, config_params)

    if get_verbosity(config_params) >= 1:
        print '\tConditions and strains with low counts'
    dataset, filtered_degree_condition_table, filtered_degree_barcodes = filter_dataset_for_count_degree(dataset, config_params, sample_table)
    if get_verbosity(config_params) >= 2:
        print dataset[2].shape
    write_count_degree_excluded_conditions(filtered_degree_condition_table, config_params)
    write_count_degree_excluded_strains(filtered_degree_barcodes, config_params)

    # Dump the dataset out to file
    dump_filtered_count_matrix(config_params, dataset)

# call: python filter_final_count_matrix.py <config_file>
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: python filter_final_count_matrix.py <config_file>'
    else:
        config_file = sys.argv[1]
        main(config_file)
