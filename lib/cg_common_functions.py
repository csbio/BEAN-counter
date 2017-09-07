#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This library contains functions that are used across
# the barseq_counter pipeline.

import os
import pandas as pd
from datetime import datetime

def get_verbosity(config_params):
 
    v = config_params['verbosity']
    if v.isdigit():
        v = int(v)
    else:
        v = 0
    return v

def get_temp_clustergram_name(output_folder, name):

    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    return os.path.join(output_folder,'{}_{}'.format(timestamp, name))

def get_sample_table(config_params, required_columns = ['screen_name', 'expt_id', 'include?', 'control?']):
    return read_sample_table(config_params['sample_table_file'], required_columns)

def read_sample_table(tab_filename, required_columns = ['screen_name', 'expt_id']):
    tab = pd.read_table(tab_filename, dtype = 'S')
    assert all([x in tab.columns for x in required_columns]), 'One or more of the required columns were not found in the sample table. If it looks like all of the columns are there, check for unwanted spaces in the columns names.\nThe missing required columns are: {}\nThe sample table is here: {}'.format([x for x in required_columns if x not in tab.columns], tab_filename)
    assert not any(tab[['screen_name', 'expt_id']].duplicated()), 'Duplicated combinations of "screen_name" and "expt_id" were found in the sample table. Please ensure that each pair of "screen_name" and "expt_id" are unique. The sample table is here: {}'.format(tab_filename)

    return tab

# def get_distribution_dir(output_folder):
#
#    return os.path.join(output_folder, 'for_distribution')
#
#def get_clustergram_distribution_dir(output_folder):
#
#    return os.path.join(get_distribution_dir(output_folder), 'clustergrams')

