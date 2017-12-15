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
parser.add_argument('--{}'.format(p.clobber), default = p.clobber.value, type = p.clobber.type, help = p.clobber.help)
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

    while True:
        print '\n'
        print 'Would you like to specify advanced options?'
        adv = raw_input('(y/n): ')
        if adv in ['y', 'n']:
            break
        else:
            print ''
            print 'Value must be in ("y", "n"). Please try again.'
        
    if adv == 'y':
        for param in ['verbosity', 'remove_barcode_specific_conditions', 'barcode_specific_template_correlation_cutoff',
                'remove_correlated_index_tags', 'index_tag_correlation_cutoff', 'common_primer_tolerance',
                'barcode_tolerance', 'control_detection_limit', 'sample_detection_limit', 'strain_pass_read_count',
                'strain_pass_fraction', 'condition_pass_read_count', 'condition_pass_fraction']:
            par_obj = getattr(p, param)
            par_obj.get_input()

# If not interactive mode, set all params via their command-line arguments
else:
    pass
    

# Make sure I check that each parameter has a valid value

# Deal with modifying paths of all config_file locations
working_dir = os.getcwd()
for param in ['output_directory', 'lane_location_file', 'sample_table_file', 'screen_config_folder']:
    par_obj = getattr(p, param)
    par_obj.value = os.path.join(working_dir, par_obj.value)

# Check if location files/folders already exist, and quit if --clobber is False
print '\n\n'
p.clobber.get_input()
print '\n\n'

f_exists_strings = []
for param in ['config_file', 'output_directory', 'lane_location_file', 'sample_table_file', 'screen_config_folder']:
    par_obj = getattr(p, param)
    if os.path.isdir(par_obj.value) or os.path.isfile(par_obj.value):
        f_exists_strings.append('{}: {}'.format(par_obj.name, par_obj.value))

if not p.clobber.value:
    if len(f_exists_strings) > 0:
        assert False, '\n\nThe following files/folders exist and cannot be overwritten unless '\
                '"--clobber" is specified: {}'.format('\n' + '\n'.join(f_exists_strings) + '\n')

# Deal with moving screen config folder to the working dir

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

        f.write('\n\n')
        f.write('####  Advanced parameters  ####\n\n')
        f.write('# Do not change unless you are comfortable doing so :-)\n')

        for param in ['verbosity', 'remove_barcode_specific_conditions', 'barcode_specific_template_correlation_cutoff',
                'remove_correlated_index_tags', 'index_tag_correlation_cutoff', 'common_primer_tolerance',
                'barcode_tolerance', 'control_detection_limit', 'sample_detection_limit', 'strain_pass_read_count',
                'strain_pass_fraction', 'condition_pass_read_count', 'condition_pass_fraction']:
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


