#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

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

def get_temp_clustergram_name(output_folder, name):

    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    return os.path.join(output_folder,'{}_{}'.format(timestamp, name))

# def get_distribution_dir(output_folder):
#
#    return os.path.join(output_folder, 'for_distribution')
#
#def get_clustergram_distribution_dir(output_folder):
#
#    return os.path.join(get_distribution_dir(output_folder), 'clustergrams')

