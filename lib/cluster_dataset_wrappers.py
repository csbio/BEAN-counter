# Here I wrap the core dataset clustering function with functions
# that send specific datasets to it and define specific output
# locations for the cluster files.
import os,sys
import numpy as np

import config_file_parser as cfp

import cluster_dataset as clus


def get_sample_table(config_params):

    filename = config_params['sample_table_file']

    # Read everything in as a string, to prevent vexing
    # number interpretation problems! Methods further down
    # can coerce to different types.
    tab = pd.read_table(filename, dtype = 'S')
    return tab
    
def get_barcode_table(config_params):

    species_config_params = get_species_config_params(config_params)

    filename = species_config_params['gene_barcode_file']
    full_path = os.path.join('./data/barcodes/{0}'.format(filename))

    tab = pd.read_table(full_path, dtype = 'S')
    return tab

def get_species_config_params(config_params):
    species_config_file = './data/species_config_file.txt'
    all_species_params = cfp.parse_species_config(species_config_file)
    species_id = config_params['species_ID']
    species_params = all_species_params[species_id]
    return species_params

def get_lane_data_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'intermediate', lane_id)

def get_lane_interactions_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'interactions', lane_id)

def get_detection_limits(config_params):

    sample_detection_limit = float(config_params['sample_detection_limit'])
    control_detection_limit = float(config_params['control_detection_limit'])

    return sample_detection_limit, control_detection_limit

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

def get_clustered_count_matrix_filename(config_params, lane_id):

    lane_folder = get_lane_data_path(config_params, lane_id)
    return os.path.join(lane_folder, '{}_barseq_matrix'.format(lane_id))

def get_dumped_zscore_matrix_filename(config_params, lane_id):

    lane_interactions_path = get_lane_interactions_path(config_params, lane_id)
    return os.path.join(lane_interactions_path, '{}_scaled_dev.dump.gz'.format(lane_id))
    
def load_dumped_zscore_matrix(config_params, lane_id):

    filename = get_dumped_zscore_matrix_filename(config_params, lane_id)
    f = gzip.open(filename, 'rb')
    barcodes, conditions, matrix = cPickle.load(f)
    dataset = [np.array(barcodes), np.array(conditions), matrix]
    f.close()

    return dataset

def get_clustered_zscore_matrix_filename(config_params, lane_id):

    lane_interactions_path = get_lane_interactions_path(config_params, lane_id)
    return os.path.join(lane_interactions_path, '{}_scaled_dev'.format(lane_id))

def customize_strains(strains, config_params, fmt_string):

    barcode_table = get_barcode_table(config_params)
    barcode_table_cols = barcode_table.columns.values
    
    cols_with_commas = [',' in x for x in barcode_table_cols]
    if np.any(cols_with_commas):
        num_with_commas = np.sum(cols_with_commas)
        raise ColumnError('{} barcode table column names contain commas,\nwhich must be addressed before visualizing'.format(num_with_commas))

    custom_columns = fmt_string.split(',')
    custom_barcode_table = barcode_table[custom_columns].set_index(['Barcode', 'Strain_ID'])
    strain_indices = [tuple(strain.split('_')) for strain in strains]
    
    custom_strains = []
    for strain in strain_indices:
        custom_strain = custom_barcode_table.ix[strain].values
        custom_strain_strings = [str(x) for x in custom_strain]
        custom_strain_final = '_'.join(custom_strain_strings)
        custom_strains.append(custom_strain_final)

    return np.array(custom_strains)

def customize_conditions(conditions, config_params, fmt_string):

    sample_table = get_sample_table(config_params)
    sample_table_cols = sample_table.columns.values
    
    cols_with_commas = [',' in x for x in sample_table_cols]
    if np.any(cols_with_commas):
        num_with_commas = np.sum(cols_with_commas)
        raise ColumnError('{} sample table column names contain commas,\nwhich must be addressed before visualizing'.format(num_with_commas))

    custom_columns = fmt_string.split(',')
    custom_sample_table = sample_table[custom_columns].set_index(['screen_name', 'expt_id'])
    condition_indices = [tuple(cond.split('-')) for cond in conditions]
    
    custom_conditions = []
    for cond in condition_indices:
        custom_cond = custom_sample_table.ix[cond].values
        custom_cond_strings = [str(x) for x in custom_cond]
        custom_cond_final = '_'.join(custom_cond_strings)
        custom_conditions.append(custom_cond_final)

    return np.array(custom_conditions)


# Define an error class to handle when the sample table columns
# have commas
class ColumnError(ColError):
    '''
    Raise when a table's columns contain commas
    and the code is trying to reference those
    columns to customize row/column names in
    the clustering output.
    '''

# In the future, these functions will take in inputs to customize
# the output (the CDT labels)
# In the future, these functions will take in inputs to customize
# the output (the CDT labels)
def cluster_count_matrix(config_file, lane_id, strain_fmt_string, cond_fmt_string):

    config_params = cfp.parse(config_file)

    sample_detection_limit, control_detection_limit = get_detection_limits(config_params)

    # If the file does not exist, then do not attempt to cluster it!
    try:
        genes, conditions, matrix = load_dumped_count_matrix(config_params, lane_id)
    except IOError:
        return None

    thresholded_matrix = matrix
    
    thresholded_matrix[thresholded_matrix < sample_detection_limit] = sample_detection_limit
    logged_matrix = np.log2(thresholded_matrix)

    # Customize the strain and condition names for interpretable visualization!
    custom_genes = customize_strains(genes, config_params, strain_fmt_string)
    custom_conditions = customize_conditions(conditions, config_params, condition_fmt_string)

    dataset = [custom_genes, custom_conditions, logged_matrix]

    record, rows_tree, cols_tree = clus.cluster(dataset)

    f = get_clustered_count_matrix_filename(config_params, lane_id)
    record.save(f, rows_tree, cols_tree)


def cluster_zscore_matrix(config_file, lane_id, strain_fmt_string, cond_fmt_string):

    config_params = cfp.parse(config_file)

    # If the file does not exist, then do not attempt to cluster it!
    try:
        genes, conditions, matrix = load_dumped_zscore_matrix(config_params, lane_id)
    except: IOError:
        return None
    
    # Customize the strain and condition names for interpretable visualization!
    custom_genes = customize_strains(genes, config_params, strain_fmt_string)
    custom_conditions = customize_conditions(conditions, config_params, condition_fmt_string)

    dataset = [custom_genes, custom_conditions, matrix]
    
    record, rows_tree, cols_tree = clus.cluster(dataset)

    f = get_clustered_zscore_matrix_filename(config_params, lane_id)
    record.save(f, rows_tree, cols_tree)




