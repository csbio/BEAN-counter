#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

import os, sys

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."

sys.path.append(os.path.join(barseq_path, 'master_scripts'))
sys.path.append(os.path.join(barseq_path, 'lib'))

from testing_lib import get_data_dir, get_results_dir, copy_data_to_results_dir
import process_screen

# How do I make sure I'm running from the correct directory?!
# Or do I find a way to manipulate the config file...
#def test_filename():
#    print "running test_filename"
#    testing_lib.get_testing_dirname()
#    print __file__

data_dir = get_data_dir()
results_dir = get_results_dir()

copy_data_to_results_dir(data_dir, results_dir)

os.chdir(results_dir)

config_file = os.path.join(results_dir, 'config_files', 'config_file.yaml')

process_screen.main(config_file, 1, 1)

def test_folder():
    print data_dir
    print results_dir



