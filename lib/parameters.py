#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

import os

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

    # Add function to control input behavior here!
    # def get_input(self):


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
        Param(name = 'screen_config_file',
            value = None,
            type = str,
            help = 'File that contains the information on the structure of the '\
                    'sequenced PCR product(s), which is required for '\
                    'correctly parsing the sequencing data. '\
                    'Also points to the barcode file that maps barcodes '\
                    'to strains and their identifiers.',
            options = None))





