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
loc_list_noconfig = ['output_directory', 'lane_location_file', 'sample_table_file', 'screen_config_folder']
sample_tab_list = ['screen_name', 'plate_size', 'plates_per_lane', 'num_lanes', 'extra_columns']
bas_list = ['verbosity', 'sub_screen_column']
adv_list = ['remove_barcode_specific_conditions', 'barcode_specific_template_correlation_cutoff',
        'remove_correlated_index_tags', 'index_tag_correlation_cutoff', 'common_primer_tolerance',
        'barcode_tolerance', 'control_detection_limit', 'sample_detection_limit', 'strain_pass_read_count',
        'strain_pass_fraction', 'condition_pass_read_count', 'condition_pass_fraction']

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--interactive', action = 'store_true', help = 'Activates interactive mode.')
parser.add_argument('--{}'.format(p.clobber.name), action = 'store_true', help = p.clobber.help)

def add_arg_group(prsr, params, param_list, name):
    grp = prsr.add_argument_group(name)
    for param in param_list:
        par_obj = getattr(params, param)
        grp.add_argument('--{}'.format(par_obj.name), default = par_obj.value, type = par_obj.type, help = par_obj.help)
    return None

add_arg_group(parser, p, loc_list, 'file/folder locations')
add_arg_group(parser, p, sample_tab_list, 'sample table parameters')
add_arg_group(parser, p, bas_list, 'basic parameters')
add_arg_group(parser, p, adv_list, 'advanced parameters')

args = parser.parse_args()


# Remainder of imports
import shutil

# Interactive mode!
if args.interactive:

    # Location parameters
    print '\n\n\n'
    print '#####################################################'
    print '######           Location parameters           ######'
    print '#####################################################'
    print '\n'
    for param in loc_list:
        par_obj = getattr(p, param)
        par_obj.get_input()

    # Sample table parameters
    print '\n\n\n'
    print '#####################################################'
    print '######         Sample table parameters         ######'
    print '#####################################################'
    print '\n'
    for param in sample_tab_list:
        par_obj = getattr(p, param)
        par_obj.get_input()

    print '\n\n\n'
    print '#####################################################'
    print '######            Basic parameters             ######'
    print '#####################################################'
    print '\n'
    for param in bas_list:
        par_obj = getattr(p, param)
        par_obj.get_input()


    # Advanced parameters
    while True:
        print '\n'
        print 'Would you like to specify advanced parameters?'
        adv = raw_input('(y/n): ')
        if adv in ['y', 'n']:
            break
        else:
            print ''
            print 'Value must be in ("y", "n"). Please try again.'
        
    print '\n\n\n'
    print '#####################################################'
    print '######           Advanced parameters           ######'
    print '#####################################################'
    print '\n'
    if adv == 'y':
        for param in adv_list:
            par_obj = getattr(p, param)
            par_obj.get_input()

# If not interactive mode, set all params via their command-line arguments
# (This is a bit circular, but it should work out).
else:
    def params_from_args(params, arguments, param_list):
        for param in param_list:
            par_obj = getattr(params, param)
            par_obj.value = getattr(arguments, param)

    params_from_args(p, args, loc_list + sample_tab_list + bas_list + adv_list)
        

# Check that each parameter has a valid value
# I **could** also add these types of checks to the argparse configuration, but
# this is probably not the best use of my time. Error messages here give the same
# if not better information to the user.
def get_invalid_params(params, param_list):
    invalid_param_list = []
    for param in param_list:
        par_obj = getattr(params, param)
        if not par_obj.has_valid_value():
            invalid_param_list.append(par_obj)
    return invalid_param_list

## Adding some bad params as a test
#p.verbosity.value = 4
#p.num_lanes.value = 'a'

invalid_param_list = get_invalid_params(p, loc_list + sample_tab_list + bas_list + adv_list)
if len(invalid_param_list) > 0:
    string_list = []
    for x in invalid_param_list:
        string_list.append('{}: {}'.format(x.name, x.value))
        if (x.options is not None) and (x.options != '_any_'):
            string_list.append('    options: {}'.format('\n             '.join([str(y) for y in x.options])))
        string_list.append('')
    string = '\n'.join(string_list)
    assert False, '\n\nThe following parameters do not possess valid values '\
            '(parameter: value; options below if pre-specified):\n\n{}\n\n'.format(string)


# Deal with modifying paths of all config_file locations
working_dir = os.getcwd()
for param in loc_list:
    par_obj = getattr(p, param)
    par_obj.value = os.path.join(working_dir, par_obj.value)

# Check if location files already exist, and quit if they do
# and --clobber is False
if args.interactive:
    print '\n\n\n'
    print '#####################################################'
    print '######         File writing parameters         ######'
    print '#####################################################'
    print '\n'
    p.clobber.get_input()
    print '\n\n'
else:
    params_from_args(p, args, ['clobber'])

f_exists_strings = []
for param in loc_list:
    par_obj = getattr(p, param)
    if os.path.isfile(par_obj.value) or (os.path.isdir(par_obj.value) and len(os.listdir(par_obj.value)) > 0):
        f_exists_strings.append('{}: {}'.format(par_obj.name, par_obj.value))

if not p.clobber.value:
    if len(f_exists_strings) > 0:
        assert False, '\n\nThe following file(s)/folder(s) exist and cannot be overwritten unless '\
                '"--clobber" is specified: {}'.format('\n' + '\n'.join(f_exists_strings) + '\n')

# Config file-writing function
def write_config_file(params, location_list, basic_list, advanced_list):

    working_dir = os.getcwd()
    fname = params.config_file.value
    with open(fname, 'wt') as f:

        f.write('---\n\n')
        f.write('##########   BEAN-counter configuration file   #########\n')
        f.write('\n')
        f.write('####  File/folder locations  ####\n')
        f.write('\n')

        for param in location_list:
            getattr(params, param).write_config(f)
            f.write('\n')

        f.write('\n\n')
        f.write('####  Basic parameters  ####\n\n')

        for param in basic_list:
            getattr(params, param).write_config(f)
            f.write('\n')

        f.write('\n\n')
        f.write('####  Advanced parameters  ####\n\n')
        f.write('# Do not change unless you are comfortable doing so :-)\n')

        for param in advanced_list:
            getattr(params, param).write_config(f)
            f.write('\n')

        f.write('...\n')

# Creating directories for each lane
def get_formatted_lanes(num_lanes):
    max_digits = len(str(num_lanes))
    return ['lane{:0{dig}}'.format(i, dig = max_digits) for i in range(1, num_lanes + 1)]

def get_raw_dir():
    working_dir = os.getcwd()
    return os.path.join(working_dir, 'raw')

def write_raw_dirs(params):
    lanes = get_formatted_lanes(params.num_lanes.value)
    raw_dir = get_raw_dir()
    if not os.path.isdir(raw_dir):
        os.makedirs(raw_dir)
    for lane in lanes:
        lane_dir = os.path.join(raw_dir, lane)
        if not os.path.isdir(lane_dir):
            os.makedirs(lane_dir)

# Lane location file-writing function
def write_lane_location_file(params):
    lanes = get_formatted_lanes(params.num_lanes.value)
    fname = params.lane_location_file.value
    raw_dir = get_raw_dir()
    with open(fname, 'wt') as f:
        f.write('lane\tlocation\n')
        for lane in lanes:
            f.write('{}\t{}\n'.format(lane, os.path.join(raw_dir, lane)))
    return None

def write_sample_table(params):
    
    n_lns = params.num_lanes.value
    plts_p_ln = params.plates_per_lane.value
    plt_sz = params.plate_size.value

    lns = get_formatted_lanes(n_lns)

    n_smpls = plt_sz * plts_p_ln * n_lns
    max_digits = len(str(n_smpls)) + 2

    columns = ['screen_name', 'expt_id', 'name', 'include?', 'control?', 'lane']
    # First properly format the list of extra columns
    n_ex_cols = 0
    if p.extra_columns.value is not None:
        ex_cols_raw = p.extra_columns.value.split(',')
        # Remove flanking whitespace, and replace any internal spaces
        # with underscores instead of skipping or throwing an error
        ex_cols = [x.strip().replace(' ', '_') for x in ex_cols_raw]
    else:
        ex_cols = []

    # Add sub_screen_column to set of extra columns (replace any spaces)
    print 'sub_screen_column: ', p.sub_screen_column.value
    if p.sub_screen_column.value is not None:
        ex_cols.append(p.sub_screen_column.value.replace(' ', '_'))

    # Now add extra columns that are unique
    for col in ex_cols:
        if col not in columns:
            columns.append(col)

    print columns

    fname = params.sample_table_file.value
    fname_parent = os.path.dirname(fname)
    if not os.path.isdir(fname_parent):
        os.makedirs(fname_parent)
    with open(fname, 'wt') as f:
        f.write('\t'.join(columns) + '\n')
        n = 0
        for i in range(n_lns):
            for j in range(plts_p_ln):
                for k in range(plt_sz):
                    n += 1
                    line = [params.screen_name.value,
                            '1{:0{dig}}'.format(n, dig = max_digits),
                            '',
                            'True',
                            'False',
                            lns[i]] + [''] * len(ex_cols)
                    f.write('\t'.join(line) + '\n')
    
    print params.screen_name.value
    print n
    
    return None


def copy_screen_config(params):
    barseq_path = os.getenv('BARSEQ_PATH')
    assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
    sc_root_dir = os.path.join(barseq_path, 'data', 'screen_configs')
    final_dir = params.screen_config_folder.value
    dir_basename = os.path.basename(final_dir)
    orig_dir = os.path.join(sc_root_dir, dir_basename)

    # Assuming that if the script has gotten this far,
    # it's okay to overwrite the existing screen config
    # directory...
    if os.path.isdir(final_dir):
        shutil.rmtree(final_dir)
    
    shutil.copytree(orig_dir, final_dir)

    return None

######  Final directory creation and file writing steps!!!  ######

# Create directory structure and generate necessary files
if not os.path.isdir('config_files'):
    os.makedirs('config_files')

# Write config file
write_config_file(p, loc_list_noconfig, bas_list, adv_list)

# Create output directory
if not os.path.isdir(p.output_directory.value):
    os.makedirs(p.output_directory.value)

# Write lane location file and raw data directory
write_lane_location_file(p)
write_raw_dirs(p)

# Write sample table file
write_sample_table(p)

# Copy screen config directory over
copy_screen_config(p)

# Print congratulatory message with further instructions
#final_string = 'Congratulations! You have successfully set up your screen for processing '\
#        'with BEAN-counter. Please see below for further instructions:\n\n'\
#        'Move, copy, or symlink your sequencing data to the respective lane folders in the raw data folder:\n{raw}\n\n'\
#        'A template sample information table has been provided here, containing the maximum '\
#        'possible number of conditions given the number of conditions per plate, plates per lane, '\
#        'and total lanes: {sam}\n'\
#        'If you have already prepared a sample information table'

