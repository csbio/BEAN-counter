# This library contains functions that are used across
# the barseq_counter pipeline.

import os
from datetime import datetime

def get_verbosity(config_params):
 
    v = config_params['verbosity']
    if v.isdigit():
        v = int(v)
    else:
        v = 0
    return v

def get_temp_clustergram_name(config_params, name):

    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    clustergram_dist_dir = get_clustergram_distribution_dir(config_params)
    return os.path.join(clustergram_dist_dir,'{}_{}'.format(timestamp, name))

def get_distribution_dir(config_params):

    return os.path.join(config_params['output_folder'], 'for_distribution')

def get_clustergram_distribution_dir(config_params):

    return os.path.join(get_distribution_dir(config_params), 'clustergrams')

