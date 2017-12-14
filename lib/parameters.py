#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

import os
import textwrap
from cg_common_functions import read_screen_config_params, read_barcode_table

# Contains all default parameter values for filling in the
# config file.

# First define the Param class
class Param:

    def __init__(self, name, value, type, help, options):
        self.name = name
        self.value = value
        self.type = type
        self.help = help
        if isinstance(options, (basestring, int, float)):
            self.options = [options]
        elif isinstance(options, (list, tuple)):
            self.options = list(options)
        elif options is None:
            self.options = options
        else:
            assert False, 'Not a valid value for the "options" argument.'

    # Function to standardize the input value to its
    # intended value (aka: if it's an integer, grab the
    # value at that index in the list of options, BUT
    # give priority to matching on the real values).
    def _standardize_value(self, val):
        if self.options is not None:
            if val in self.options:
                return val
            elif val.isdigit():
                if int(val) in range(len(self.options)):
                    return self.options[int(val)]
            else:
                return val

    # Function to check for a valid "value" value
    def _is_valid_value(self, val):
        if self.options is not None:
            return val in self.options
        else:
            return val not in [None, '']

    def has_valid_value(self):
        return self._is_valid_value(self.value)

    def _write_help_string(self):
        # Setup
        opt_string = 'Type "o" to see options and "h" to print help again.'
        #help_string = '{} {} (default: {})'.format(self.help, opt_string, self.value)
        # Printing
        print '\n' * 4
        print '#' * (len(self.name) + 8)
        print '##  {}  ##'.format(self.name)
        print '#' * (len(self.name) + 8)
        print ''
        for x in textwrap.wrap(self.help, width = 60):
            print x
        print 'default: {}'.format(self.value)
        print ''
        print opt_string
        print ''
    
    # Add function to control input behavior here!
    def get_input(self):

        self._write_help_string()

        while True:
            y = raw_input('{}: '.format(self.name))
            if y == 'o':
                if self.options is not None:
                    for i, opt in enumerate(self.options):
                        print ''
                        print 'Options for {}:'.format(self.name)
                        print '(specify via either name or number)'
                        print '{}: {}'.format(i, opt)
                        print
                else:
                    print ''
                    print 'No pre-defined options for "{}".'.format(self.name)
                    print ''
            elif y == 'h':
                self._write_help_string()
            else:
                y = self._standardize_value(y)
                if self._is_valid_value(y):
                    self.value = y
                if self.has_valid_value():
                    break
                else:
                    print ''
                    if not self._is_valid_value(y):
                        print '"{}" is not a valid value for "{}".'.format(y, self.name)
                    print 'No default value exists for "{}", so the user must supply one!'.format(self.name)
                    print ''

    # Function to spit param out to config file
    def write_config(self, f):
        f.writelines(['# {}\n'.format(x) for x in textwrap.wrap(self.help, width = 60)])
        f.write('{}: {}\n'.format(self.name, self.value))



# Detection of valid screen_config directories in the
# 'data/screen_config/' directory.
def is_valid_screen_config_dir(folder):

    files = os.listdir(folder)
    print 'files:', files
    if 'screen_config.yaml' in files:
        sc_params = read_screen_config_params(os.path.join(folder, 'screen_config.yaml'))
        if 'gene_barcode_file' in sc_params:
            print 'gene barcode file:', sc_params['gene_barcode_file']
            try:
                read_barcode_table(sc_params['gene_barcode_file'])
            except Exception as e:
                print 'Assertion failed'
            else:
                print 'Assertion passed'

            # If no exception was raised, then return True!
            if 'e' not in locals():
                print 'returning True...'
                return True

    # If anything fails, return False
    return False

def list_valid_screen_config_dirs():
    barseq_path = os.getenv('BARSEQ_PATH')
    assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
    sc_root_dir = os.path.join(barseq_path, 'data', 'screen_configs')
    print 'root dir:', sc_root_dir
    sc_dirs = [x for x in os.listdir(sc_root_dir) if os.path.isdir(os.path.join(sc_root_dir, x))]
    print 'screen config dirs:', sc_dirs
    valid_sc_dirs = [x for x in sc_dirs if is_valid_screen_config_dir(os.path.join(sc_root_dir, x))]
    print 'valid screen config dirs:', valid_sc_dirs
    return valid_sc_dirs

print list_valid_screen_config_dirs()

# Definitions of file/folder locations for the config file
# (and of the config file itself!)
#loc_list = []
config_file = Param(
        name = 'config_file',
        value = os.path.join('config_files', 'config_file.yaml'),
        type = str,
        help = 'File that coordinates all aspects of interaction scoring.',
        options = None)

output_directory = Param(
        name = 'output_directory',
        value = 'output',
        type = str,
        help = 'Directory to which output is written',
        options = None)

lane_location_file = Param(
        name = 'lane_location_file',
        value = os.path.join('config_files', 'lane_locations.txt'),
        type = str,
        help = 'File that maps the lane column in the sample table ' \
                'to the folder containing raw sequencing data for that lane.',
        options = None)

sample_table_file = Param(
        name = 'sample_table_file',
        value = os.path.join('sample_table_files', 'sample_table.txt'),
        type = str,
        help = 'File that contains all sample information, including unique IDs, '\
                'negative control status, the sequence ' \
                'of the index tag(s) each sample was labeled with and the '\
                'sequencing lane from which the data came.',
        options = None)

screen_config_folder = Param(
        name = 'screen_config_folder',
        value = None,
        type = str,
        help = 'Folder in "$BARSEQ_PATH/data/screen_configs/" that contains: ' \
                '1) a file that defines the structure of the '\
                'sequenced PCR product(s), which is required for '\
                'correctly parsing the sequencing data (screen_config.yaml), and 2) '\
                'the file that maps barcodes '\
                'to strains and their identifiers (barcodes.txt).',
        options = list_valid_screen_config_dirs())



