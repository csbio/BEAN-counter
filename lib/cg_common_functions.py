#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.4.0'

# This library contains functions that are used across
# the barseq_counter pipeline.

import os
import pandas as pd, numpy as np
from datetime import datetime
import yaml
from multiprocessing import cpu_count

bool_dict = {'True': True, 'TRUE': True, 'False': False, 'FALSE': False,
        'T': True, 'F': False, '1': True, '0': False, True: True, False: False}

def get_verbosity(config_params):
 
    v = config_params['verbosity']
    if isinstance(v, int):
        return v
    else:
        return 0

#yaml = YAML()

def parse_yaml(filename):

    with open(filename, 'rt') as f:
        try:
            params = yaml.load(f)
        except yaml.YAMLError as e:
            assert False, 'yaml screen config file did not load properly.\nThe file is: {}\nOriginal error:\n{}'.format(filename, e)

    return params

def get_sample_table(config_params, required_columns = ['screen_name', 'expt_id', 'include?', 'control?']):
    return read_sample_table(config_params['sample_table_file'], required_columns)

def read_sample_table(tab_filename, required_columns = ['screen_name', 'expt_id']):
    tab = pd.read_table(tab_filename, dtype = 'S', comment = '#')
    assert all([x in tab.columns for x in required_columns]), 'One or more of the required columns were not found in the sample table. If it looks like all of the columns are there, check for unwanted spaces in the column names.\nThe missing required columns are: {}\nThe sample table is here: {}'.format([x for x in required_columns if x not in tab.columns], tab_filename)
    assert not any(tab[['screen_name', 'expt_id']].duplicated()), 'Duplicated combinations of "screen_name" and "expt_id" were found in the sample table. Please ensure that each pair of "screen_name" and "expt_id" are unique. The sample table is here: {}'.format(tab_filename)
    return tab

def get_amplicon_struct_params(config_params):
    return read_amplicon_struct_params(config_params['amplicon_struct_file'])

def read_amplicon_struct_params(filename):
    amplicon_struct_params = parse_yaml(filename)
    #amplicon_struct_params['gene_barcode_file'] = os.path.join(os.path.dirname(filename),
    #        amplicon_struct_params['gene_barcode_file'])
    return amplicon_struct_params

def get_barcode_table(config_params):
    tab = read_barcode_table(config_params['gene_barcode_file'], required_columns = ['Strain_ID', 'include?'])
    return tab

def read_barcode_table(tab_filename, required_columns = ['Strain_ID'], comment = '#'):
    tab = pd.read_table(tab_filename, dtype = 'S', comment = comment)
    assert all([x in tab.columns for x in required_columns]), 'One or more of the required columns were not found in the barcode table. If it looks like all of the columns are there, check for unwanted spaces in the column names.\nThe missing required columns are: {}\nThe barcode table is here: {}'.format([x for x in required_columns if x not in tab.columns], tab_filename)
    # This will happen at some point, but not ready for prime time yet
    # until I get the strain identifier thing worked out (aka not
    # dependent on barcodes anymore).
    assert not any(tab['Strain_ID'].duplicated()), 'Duplicated values were found for within the "Strain_ID" column of the barcode table. Please ensure that all of these values are unique. The barcode table is here: {}'.format(tab_filename)
    return tab

def get_num_cores(config_params):
    '''
    Use integers to specify the exact number of cores.
    Use floats to specify the proportion of the
    processor's cores to use. Use "all" to specify all
    cores. Defaults to the number of detected cores
    divided by 2.
    '''
    n_raw = config_params.get('num_cores')
    max_cores = cpu_count()
    if n_raw is None:
        n = int(np.ceil(max_cores / 2.))
    elif n_raw == 'all':
        n = max_cores
    elif n_raw == 1:
        n = 1
    elif isinstance(n_raw, int):
        n = n_raw
    elif isinstance(n_raw, float):
        n = int(np.ceil(max_cores * n_raw))
    else:
        print '\nWarning: the "num_cores" parameter specified in the config file ' \
                '("{}") is not valid. Defaulting to the number of detected cores ' \
                '/ 2'.format(n_raw)
        n = int(np.ceil(max_cores / 2.))
    
    return n

def get_temp_clustergram_name(output_folder, name):

    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    return os.path.join(output_folder,'{}_{}'.format(timestamp, name))

def a_is_row_in_b(a, b):

    return np.any(np.all(a == b, axis = 1))

# def get_distribution_dir(output_folder):
#
#    return os.path.join(output_folder, 'for_distribution')
#
#def get_clustergram_distribution_dir(output_folder):
#
#    return os.path.join(get_distribution_dir(output_folder), 'clustergrams')

