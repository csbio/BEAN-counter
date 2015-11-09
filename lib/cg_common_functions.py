# This library contains functions that are used across
# the barseq_counter pipeline.

def get_verbosity(config_params):
 
    v = config_params['verbosity']
    if v.isdigit():
        v = int(v)
    else:
        v = 0
    return v

