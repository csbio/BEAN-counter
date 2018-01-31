#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.1.0'

# Helper functions for all testing scripts

import os, sys
from datetime import datetime
import inspect
import shutil

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."

def get_test_dirname(depth = 1):
    '''
    Assuming this module is imported directly into the testing script for a
    particular test case, this function returns the name of that script's
    parent directory. This is then used to determine where to put the test
    results. `depth` argument can be changed to that this function can be
    called from inside other functions in this library.
    '''

    testing_script_path = inspect.stack()[depth][1]
    testing_dirname_full = os.path.dirname(os.path.realpath(testing_script_path))
    return os.path.split(testing_dirname_full)[1]

def get_test_results_dirname():

    '''
    Assuming this module is imported directly into the testing script for a
    particular test case, this function returns the directory to which the
    test results are exported. This combines the directory that reflects the
    name of the test plus a timestamped subfolder.
    '''
    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    test_dirname = get_test_dirname(depth = 2)
    path = os.path.join(barseq_path, 'test_results', test_dirname, timestamp)
    if not os.path.isdir(path):
        os.makedirs(path)

    return path

def get_data_dir():

    '''
    This function should be called from within a testing script but not
    from inside another function in that testing script.
    '''
    testing_script_path = inspect.stack()[1][1]
    return os.path.join(os.path.dirname(os.path.realpath(testing_script_path)), 'data')

def get_results_dir():

    '''
    This function should be called from within a testing script but not
    from inside another function in that testing script.
    '''
    timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    testing_script_path = inspect.stack()[1][1]
    testing_dirname_full = os.path.dirname(os.path.realpath(testing_script_path))
    test_dirname = os.path.split(testing_dirname_full)[1]
    return os.path.join(barseq_path, 'test_results', test_dirname, timestamp)


def copy_data_to_results_dir(data_dir, results_dir):

    shutil.copytree(data_dir, results_dir)

def get_file():

    print __file__

#def set

