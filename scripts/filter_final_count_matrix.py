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

def filter_dataset_for_include_2(dataset, sample_table):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
    
    include_table = sample_table[sample_table['include?'].astype(np.bool)]
    not_include_table = sample_table[np.invert(sample_table['include?'].astype(np.bool))]
    include_screen_names = include_table['screen_name']
    include_expt_ids = include_table['expt_id']
    include_condition_ids = ['{0}-{1}'.format(*x) for x in it.izip(include_screen_names, include_expt_ids)]
    
    include_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if cond_id in include_condition_ids])
    
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
        gene_barcode_ids_to_remove = [line.rstrip().replace('\t', '_') for line in f]
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
    inds_to_keep = np.array([i for i, bg_id in enumerate(barcode_gene_ids) if bg_id not in gene_barcode_ids_to_remove])

    filtered_barcode_gene_ids = barcode_gene_ids[inds_to_keep]
    filtered_matrix = matrix[inds_to_keep, :]
   
    return [filtered_barcode_gene_ids, condition_ids, filtered_matrix], gene_barcode_ids_to_remove

def get_index_tag_correlation_path(config_params):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'index_tag_correlations')

def filter_dataset_for_index_tags(dataset, config_params):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset

    index_tag_corr_path = get_index_tag_correlation_path(config_params)
    index_tag_corr_filename = os.path.join(index_tag_corr_path, 'control_index_tag_correlations.dump')
    
    index_tags, correlations = cPickle.load(index_tag_corr_filename)
    index_tags_to_remove = [tag[i] for i,corr in enumerate(correlations) if corr >= index_tag_correlation_threshold]
    index_tags_to_keep = [tag[i] for i,corr in enumerate(correlations) if corr < index_tag_correlation_threshold]
   
    sample_table = sample_table.set_index('index_tag')
    to_remove_table = sample_table.ix[index_tags_to_remove]
    to_keep_table = sample_table.ix[index_tags_to_keep]

    to_keep_screen_names = to_keep_table['screen_name']
    to_keep_expt_ids = to_keep_table['expt_id']
    to_keep_condition_ids = ['{0}-{1}'.format(*x) for x in it.izip(to_keep_screen_names, to_keep_expt_ids)]

    to_keep_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if cond_id in to_keep_condition_ids])
    
    filtered_condition_ids = condition_ids[to_keep_condition_indices]
    filtered_matrix = matrix[:, to_keep_condition_indices]

    # ALSO return the sample table of filtered out conditions!
    to_remove_table = to_remove_table.reset_index()
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
    print passing_strain_fraction
    strains_to_keep_inds = passing_strain_fraction >= strain_pass_fraction
    strains_to_remove_inds = np.invert(strains_to_keep_inds)
    strains_to_remove = barcode_gene_ids[strains_to_remove_inds]

    # Get the indices of the conditions to keep (and those to eliminate)
    passing_condition_read_count_matrix = matrix >= condition_pass_read_count
    passing_condition_read_count_sums = np.nansum(passing_condition_read_count_matrix, axis = 0)
    passing_condition_fraction = passing_condition_read_count_sums / float(passing_condition_read_count_matrix.shape[0])
    print passing_condition_fraction
    conditions_to_keep_inds = passing_condition_fraction >= condition_pass_fraction
    conditions_to_remove_inds = np.invert(conditions_to_keep_inds)
    conditions_to_remove = condition_ids[conditions_to_remove_inds]
    conditions_to_remove_split = [tuple(x.split('-')) for x in conditions_to_remove]

    sample_table = sample_table.set_index(['screen_name', 'expt_id'])
    conditions_to_remove_table = sample_table.loc[conditions_to_remove_split]

    # Return conditions and strains rejected by the count degree criteria 
    filtered_barcode_gene_ids = barcode_gene_ids[strains_to_keep_inds]
    filtered_condition_ids = condition_ids[conditions_to_keep_inds]
    filtered_matrix = matrix[strains_to_keep_inds, :]
    filtered_matrix = filtered_matrix[:, conditions_to_keep_inds]

    return [filtered_barcode_gene_ids, filtered_condition_ids, filtered_matrix], conditions_to_remove_table.reset_index(), strains_to_remove

def write_filtered_strains_conditions(config_params, filtered_include_tab, filtered_barcodes, filtered_degree_condition_table, filtered_degree_barcodes):
    filtered_path = get_filtered_strains_conditions_path(config_params)
    if not os.path.isdir(filtered_path):
        os.makedirs(filtered_path)

    # Write out conditions that were specified as "do not include"
    filtered_include_tab_file = os.path.join(filtered_path, 'prespecified_excluded_conditions.txt')
    filtered_include_tab.to_csv(filtered_include_tab_file, sep = '\t', index = False)

    # Write out the conditions excluded because of their count degree
    cond_degree_file = os.path.join(filtered_path, 'count_degree_excluded_conditions.txt')
    filtered_degree_condition_table.to_csv(cond_degree_file, sep = '\t', index = False)

    # Write out the barcodes that were specifed as "do not include"
    filtered_barcode_file = os.path.join(filtered_path, 'prespecified_excluded_strains.txt')
    with open(filtered_barcode_file, 'wt') as f:
        for barcode in filtered_barcodes:
            f.write(barcode.replace('_', '\t') + '\n')

    # Write out the barcodes excluded because of their count degree
    strain_degree_file = os.path.join(filtered_path, 'count_degree_excluded_strains.txt')
    with open(strain_degree_file, 'wt') as f1:
        for barcode in filtered_degree_barcodes:
            f1.write(barcode.replace('_', '\t') + '\n')

def dump_filtered_count_matrix(config_params, dataset):
   
    lane_id = 'all_lanes_filtered'
    dataset_path = get_lane_data_path(config_params, lane_id)
    if not os.path.isdir(dataset_path):
        os.makedirs(dataset_path)

    dataset_filename = get_dumped_count_matrix_filename(config_params, lane_id)
    dump_dataset(dataset, dataset_filename)

def main(config_file):

    # Read in the config params
    print 'parsing parameters...'
    config_params = cfp.parse(config_file)
    sample_table = get_sample_table(config_params)

    dataset = load_dumped_count_matrix(config_params, 'all_lanes')

    dataset, filtered_include_tab = filter_dataset_for_include_2(dataset, sample_table)
    dataset, filtered_barcodes = filter_dataset_for_barcodes(dataset, config_params)
    dataset, filtered_index_tag_condition_table = filter_dataset_for_index_tags(dataset, config_params)
    dataset, filtered_degree_condition_table, filtered_degree_barcodes = filter_dataset_for_count_degree(dataset, config_params, sample_table)

    # Write filtered conditions/barcodes out to file
    write_filtered_strains_conditions(config_params, filtered_include_tab, filtered_barcodes, filtered_degree_condition_table, filtered_degree_barcodes)

    # Dump the dataset out to file
    dump_filtered_count_matrix(config_params, dataset)

# call: python filter_final_count_matrix.py <config_file>
if '__name__' == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: python filter_final_count_matrix.py <config_file>'
    else:
        config_file = sys.argv[1]
        main(config_file)
