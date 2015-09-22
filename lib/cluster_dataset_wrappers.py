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



# In the future, these functions will take in inputs to customize
# the output (the CDT labels)
# In the future, these functions will take in inputs to customize
# the output (the CDT labels)
def cluster_count_matrix(config_file, lane_id):

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

    dataset = [genes, conditions, logged_matrix]

    record, rows_tree, cols_tree = clus.cluster(dataset)

    f = get_clustered_count_matrix_filename(config_params, lane_id)
    record.save(f, rows_tree, cols_tree)


def cluster_zscore_matrix(config_file, lane_id):

    config_params = cfp.parse(config_file)

    # If the file does not exist, then do not attempt to cluster it!
    try:
        dataset = load_dumped_zscore_matrix(config_params, lane_id)
    except: IOError:
        return None
    
    record, rows_tree, cols_tree = clus.cluster(dataset)

    f = get_clustered_zscore_matrix_filename(config_params, lane_id)
    record.save(f, rows_tree, cols_tree)




