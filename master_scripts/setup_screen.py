#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script provides both interactive and programmable interfaces
# for setting up configuration files in the default directory
# structure and linking all files to each other appropriately.



## Set up the Arg class
#class Arg:
#
#    def __init__(self, name, value, type, help, options):
#        self.name = name
#        self.value = value
#        self.default = default
#        self.type = type
#        self.help = help
#        self.options = options
#
#    # Add function to control input behavior here!

## Define the file/folder locations
#loc_list = []
#loc_list.append(
#        Arg(name = 'output_directory',
#            value = 'output',
#            type = str,
#            help = 'Directory to which output is written',
#            options = None))


# Minimal imports
import os, sys

# Set up path
barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'lib'))

import parameters as p

# Parse arguments
import argparse

#arg_names = [x.name for x in arg_list]
#arg_name_idx = {x:i for i,x in enumerate(arg_names)}
#assert len(arg_name_idx) == len(arg_list), "One or more arguments have the same name."

# Define some parameter sets
loc_list = ['config_file', 'output_directory', 'lane_location_file', 'sample_table_file', 'screen_config_folder']

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--interactive', action = 'store_true', help = 'Activates interactive mode.')
parser.add_argument('--clobber', action = 'store_true', help = 'Overwrite existing files.')
loc_group = parser.add_argument_group('file/folder locations')
for param in loc_list:
    par_obj = getattr(p, param)
    loc_group.add_argument('--{}'.format(par_obj.name), default = par_obj.value, type = par_obj.type, help = par_obj.help)

args = parser.parse_args()


# Remainder of imports
import textwrap


# Interactive mode!
if args.interactive:

    for param in loc_list:
        par_obj = getattr(p, param)

        par_obj.get_input()


# Make sure I check that each parameter has a valid value

# Deal with modifying paths of all config_file locations

# Deal with moving screen config folder to the working dir

# Check if location files already exist, and quit if --clobber is not set

# Config file-writing function
def write_config_file(params):

    working_dir = os.getcwd()
    fname = params.config_file.value
    with open(fname, 'wt') as f:

        f.write('---\n\n')
        f.write('##########   BEAN-counter configuration file   #########\n')
        f.write('\n')
        f.write('####  File/folder locations  ####\n')
        f.write('\n')

        for param in ['output_directory', 'lane_location_file', 'sample_table_file', 'screen_config_folder']:
            getattr(params, param).write_config(f)
            f.write('\n')

        f.write('...\n')

# Lane location file-writing function
def write_lane_location_file(fname):
    pass


######  Final directory creation and file writing steps!!!  ######

# Create directory structure and generate necessary files
if not os.path.isdir('config_files'):
    os.makedirs('config_files')

write_config_file(p)


