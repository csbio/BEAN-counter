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
        if isinstance(options, (basestring, int, float, bool)):
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

        # First, if it's blank, other functions will handle that
        if val in [None, '']:
            return val
        
        # See if the value matches the type it is required to be.
        # The coercion to string prevents things like coercing "2"
        # to boolean from giving the appearance of working.
        try:
            matches_type = str(self.type(val)) == str(val)
        except:
            matches_type = False

        # See if the value is coerceable to int and if it
        # matches the range of options given.
        matches_int = True
        try:
            int(val)
        except:
            matches_int = False


        # Put it all together here
        if self.options is None:
            if matches_type:
                return self.type(val)
            else:
                # This nonmatching value will be taken care of in subsequent steps
                return val
                #assert False, 'Value "{}" for parameter "{}" cannot be coerced to type {}, and no options are pre-defined for integer-based selection.'.format(val, self.name, self.type)
        else:
            if matches_type:
                if self.type(val) in self.options:
                    return self.type(val)
            if matches_int:
                if int(val) in range(len(self.options)):
                    return self.options[int(val)]
            # If the value could not be matched to the options either
            # via it's correctly-typed value or via integer indexing,
            # then it is not a valid value for that parameter. Return
            # it as-is to be dealt with by other validation functions.
            return val
            #assert False, 'Value "{}" could not be matched to parameter "{}" by name or by position in the list of options. It is required to be of type {}'.format(val, self.name, self.type)


    # Function to check for a valid "value" value
    def _is_valid_value(self, val):
      
        # This function must be performed on "standardized" values,
        # as it will say that integers are not valid inputs even
        # when they are used for selection from a list of parameters.

        # First, see if the value matches the required type.
        # If not, then it's not valid!
        if not isinstance(val, self.type):
            return False
        
        # If it matches type, then it must match the pre-defined options.
        if self.options is not None:
            return val in self.options
        # And if there are no pre-defined options, ensure it's not blank.
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

###########################
##  Advanced parameters  ##
###########################

verbosity = Param(
        name = 'verbosity',
        value = 1,
        type = int,
        help = 'Verbosity of messages printed to stdout. '\
                '0 is none, 1 is default (step-relevant announcements), '\
                '2 prints object inspection results, and 3 presents even '\
                'deeper (and more cryptic) inspection results.',
        options = range(4))

remove_barcode_specific_conditions = Param(
        name = 'remove_barcode_specific_conditions',
        value = True,
        type = bool,
        help = 'Remove conditions matching the "barcode-specific" pattern ' \
                '(positive interactions for barcodes beginning with one ' \
                'particular base) with correlation greater than that defined ' \
                'by "barcode_specific_template_correlation_cutoff".',
        options = [True, False])

barcode_specific_template_correlation_cutoff = Param(
        name = 'barcode_specific_template_correlation_cutoff',
        value = 0.3,
        type = float,
        help = 'See "remove_barcode_specific_conditions".',
        options = None)

remove_correlated_index_tags = Param(
        name = 'remove_correlated_index_tags',
        value = True,
        type = bool,
        help = 'Remove conditions associated with highly-self-correlating '\
                'index tags (in control conditions) with correlation greater ' \
                ' than that defined by "index_tag_correlation_cutoff". ' \
                'Intended for use with larger screens, as these are more likely ' \
                'to have coverage of multiple negative control conditions with ' \
                'the same index tag.',
        options = [True, False])

index_tag_correlation_cutoff = Param(
        name = 'index_tag_correlation_cutoff',
        value = 0.4,
        type = float,
        help = 'See "remove_correlated_index_tags".',
        options = None)

common_primer_tolerance = Param(
        name = 'common_primer_tolerance',
        value = 2,
        type = int,
        help = 'Number of substitution errors allowed in the common primer sequence.',
        options = None)

barcode_tolerance = Param(
        name = 'barcode_tolerance',
        value = 2,
        type = int,
        help = 'Number of substitution errors allowed in the barcode sequence.',
        options = None)

control_detection_limit = Param(
        name = 'control_detection_limit',
        value = 20,
        type = int,
        help = 'Read count below which strains in control conditions are '\
                'labeled as "not detected" and set to NA prior to calculation '\
                'of the mean control profile.',
        options = None)

sample_detection_limit = Param(
        name = 'sample_detection_limit',
        value = 20,
        type = int,
        help = 'Read count below which strains in non-control conditions are '\
                '"not detected." All counts below this number are set to this '\
                'number, which prevents NAs when log-transforming. Also applies '\
                'to control conditions when their normalized profiles are being computed.',
        options = None)

strain_pass_read_count = Param(
        name = 'strain_pass_read_count',
        value = 20,
        type = int,
        help = 'For a strain to pass, the fraction of conditions with read counts '\
                '>= "strain_pass_read_count" must be >= "strain_pass_fraction". '\
                'This filtering step occurs AFTER removing strains/conditions with ' \
                'zero counts.',
        options = None)

strain_pass_fraction = Param(
        name = 'strain_pass_fraction',
        value = 0.25,
        type = int,
        help = 'See "strain_pass_read_count".',
        options = None)

condition_pass_read_count = Param(
        name = 'condition_pass_read_count',
        value = 20,
        type = int,
        help = 'For a condition to pass, the fraction of conditions with read counts '\
                '>= "condition_pass_read_count" must be >= "condition_pass_fraction". '\
                'This filtering step occurs AFTER removing strains/conditions with ' \
                'zero counts.',
        options = None)

condition_pass_fraction = Param(
        name = 'condition_pass_fraction',
        value = 0.25,
        type = int,
        help = 'See "condition_pass_read_count".',
        options = None)





