#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

import os
import textwrap

# Contains all default parameter values for filling in the
# config file.

# First define the Param class
class Param:

    def __init__(self, name, value, type, help, options):
        self.name = name
        self.value = value
        self.type = type
        self.help = help
        self.options = options

    # Function to check for a valid "value" value
    def _is_valid_value(self, val):
        return val not in [None, '']

    def has_valid_value(self):
        if self.options is not None:
            return self.value in self.options
        else:
            return self._is_valid_value(self.value)

    # Add function to control input behavior here!
    def get_input(self):
        opt_string = 'Type "o" to see options.'
        help_string = '{} {} (default: {})'.format(self.help, opt_string, self.value)
        for x in textwrap.wrap(help_string, width = 60):
            print x
        while True:
            y = raw_input('{}: '.format(self.name))
            if y == 'o':
                if self.options is not None:
                    if isinstance(self.options, basestring):
                        print self.options
                    elif isinstance(self.options, (list, tuple)):
                        for opt in self.options:
                            print opt
                else:
                    print 'No pre-defined options for "{}".'.format(self.name)
            else:
                if self._is_valid_value(y):
                    self.value = y
                if self.has_valid_value():
                    break
                else:
                    print 'No default value exists for "{}", so the user must supply one!'.format(self.name)


# Detection of valid screen_config directories in the
# 'data/screen_config/' directory.
#def is_valid_screen_config_dir(folder):
#
#    files = os.listdir(folder)
#    if 'screen_config.yaml' in files:
#        sc_params = parse_yaml(os.path.join(folder, 'screen_config.yaml'))
#        if 'gene_barcode_file' in sc_params:
#            try:
#                read_barcode_table(sc_params['gene_barcode_file'])
#                return True
#            except AssertionError as e:
#                pass
#
#    # If anything fails, return False
#    return False
#
#def list_valid_screen_config_dirs():
#    barseq_path = os.getenv('BARSEQ_PATH')
#    assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
#    sc_root_dir = os.path.join(barseq_path, 'data', 'screen_configs')
#    sc_dirs = [os.path.join(sc_root_dir, x) for x in os.listdir(sc_root_dir) if os.path.isdir(x)]
#    valid_sc_dirs = [x for x in sc_dirs if is_valid_screen_config_dir(x)]
#
#print list_valid_screen_config_dirs()

# Definitions of file/folder locations for the config file
# (and of the config file itself!)
loc_list = []
loc_list.append(
        Param(name = 'config_file',
            value = os.path.join('config_files', 'config_file.yaml'),
            type = str,
            help = 'File that coordinates all aspects of interaction scoring.',
            options = None))

loc_list.append(
        Param(name = 'output_directory',
            value = 'output',
            type = str,
            help = 'Directory to which output is written',
            options = None))

loc_list.append(
        Param(name = 'lane_location_file',
            value = os.path.join('config_files', 'lane_locations.txt'),
            type = str,
            help = 'File that maps the lane column in the sample table ' \
                    'to the folder containing raw sequencing data for that lane.',
            options = None))

loc_list.append(
        Param(name = 'sample_table_file',
            value = os.path.join('sample_table_files', 'sample_table.txt'),
            type = str,
            help = 'File that contains all sample information, including unique IDs, '\
                    'negative control status, the sequence ' \
                    'of the index tag(s) each sample was labeled with and the '\
                    'sequencing lane from which the data came.',
            options = None))

loc_list.append(
        Param(name = 'screen_config_folder',
            value = None,
            type = str,
            help = 'Folder in "$BARSEQ_PATH/data/screen_configs/" that contains: ' \
                    '1) a file that defines the structure of the '\
                    'sequenced PCR product(s), which is required for '\
                    'correctly parsing the sequencing data (screen_config.yaml), and 2) '\
                    'the file that maps barcodes '\
                    'to strains and their identifiers (barcodes.txt).',
            options = None))
#            options = list_valid_screen_config_dirs()))




