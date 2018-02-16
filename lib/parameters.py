#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.2.4'

import os
import textwrap
from cg_common_functions import read_amplicon_struct_params, read_barcode_table

# Contains all default parameter values for filling in the
# config file.

# First define the Param class
class Param:

    def __init__(self, name, value, type, help, options, config_help = None):
        self.name = name
        self.value = value
        self.type = type
        self.help = help
        if config_help is None:
            self.config_help = self.help
        else:
            self.config_help = config_help
        if options == '_any_':
            self.options = options
        elif isinstance(options, (basestring, int, float, bool)):
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

        # First, if it's blank/None, there is nothing to standardize.
        # Other functions will handle that
        if val in [None, '']:
            return None
        
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
                # This type-nonmatching value will be taken care of in subsequent steps
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

        # First, if value is None, the only way this is okay is if
        # '_any_' was specified for the "options" field
        if val is None:
            return self.options == '_any_'

        # Then, see if the value matches the required type.
        # If not, then it's not valid!
        if not isinstance(val, self.type):
            return False
        
        # If it matches type, then it must match the pre-defined options
        # if they exist.
        if self.options not in [None, '_any_']:
            return val in self.options

        # Otherwise, you have a value that matches self.type that can be anything.
        return True

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
            y = y.strip()
            if y == 'o':
                print ''
                if self.options is not None:
                    print 'Options for {}:'.format(self.name)
                    print '(specify via either name or number)'
                    for i, opt in enumerate(self.options):
                        print '{}: {}'.format(i, opt)
                else:
                    print 'No pre-defined options for "{}".'.format(self.name)
                print ''
            elif y == 'h':
                self._write_help_string()
            else:
                #print 'raw val:', y
                y = self._standardize_value(y)
                #print 'standardized val:', y
                if self._is_valid_value(y):
                    #print 'is_valid_value'
                    self.value = y
                if self.has_valid_value():
                    #print 'has_valid_value'
                    break
                else:
                    print ''
                    if not self._is_valid_value(y):
                        print '"{}" is not a valid value for "{}".'.format(y, self.name)
                    print 'No default value exists for "{}", so the user must supply one!'.format(self.name)
                    print ''

    # Function to spit param out to config file
    def write_config(self, f):
        f.writelines(['# {}\n'.format(x) for x in textwrap.wrap(self.config_help, width = 60)])
        if self.value is None:
            f.write('{}: {}\n'.format(self.name, ''))
        else:
            f.write('{}: {}\n'.format(self.name, self.value))

def list_valid_amplicon_struct_files():
    barseq_path = os.getenv('BARSEQ_PATH')
    assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
    as_root_dir = os.path.join(barseq_path, 'data', 'amplicon_struct_files')
    as_files = [x for x in os.listdir(as_root_dir) if x.endswith('.yaml')]
    valid_as_files = []
    for as_file in as_files:
        try:
            read_amplicon_struct_params(os.path.join(as_root_dir, as_file))
        except Exception as e:
            pass
        else:
            valid_as_files.append(as_file)

    return valid_as_files

def list_valid_gene_barcode_files():
    barseq_path = os.getenv('BARSEQ_PATH')
    assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
    bc_root_dir = os.path.join(barseq_path, 'data', 'gene_barcode_files')
    bc_files = os.listdir(bc_root_dir)
    valid_bc_files = []
    for bc_file in bc_files:
        try:
            read_barcode_table(os.path.join(bc_root_dir, bc_file))
        except Exception as e:
            pass
        else:
            valid_bc_files.append(bc_file)

    return valid_bc_files


#print list_valid_amplicon_struct_dirs()


#########################
## Random command-line ##
## params              ##
#########################
clobber = Param(
        name = 'clobber',
        value = False,
        type = bool,
        help = 'If True, overwrites existing config and '\
                'condition/strain information files.',
        options = [False, True])

############################################################
# Definitions of file/folder locations for the config file #
# (and of the config file itself!)                         #
############################################################
config_file = Param(
        name = 'config_file',
        value = os.path.join('config_files', 'config_file.yaml'),
        type = str,
        help = 'File that coordinates all aspects of interaction scoring.',
        options = None)

output_folder = Param(
        name = 'output_folder',
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
        value = os.path.join('sample_table', 'sample_table.txt'),
        type = str,
        help = 'File that contains all sample information, including unique IDs, '\
                'negative control status, the sequence ' \
                'of the index tag(s) each sample was labeled with and the '\
                'sequencing lane from which the data came.',
        options = None)

gene_barcode_file = Param(
        name = 'gene_barcode_file',
        value = None,
        type = str,
        help = 'Tab-delimited text file in $BARSEQ_PATH/data/gene_barcode_files/ '\
                'that maps barcodes to strains and their identifiers. '\
                'Must contain unique "Strain_ID" column.',
        options = list_valid_gene_barcode_files(),
        config_help = 'Tab-delimited text file that maps barcodes to strains and ' \
                'their identifiers. Must contain unique "Strain_ID" column.')

amplicon_struct_file = Param(
        name = 'amplicon_struct_file',
        value = None,
        type = str,
        help = 'YAML-formatted file in $BARSEQ_PATH/data/amplicon_struct_files/ that ' \
                'defines the structure of the sequenced PCR amplicons. This is '\
                'required for correctly parsing the sequencing data.',
        options = list_valid_amplicon_struct_files(),
        config_help = 'YAML-formatted file that defines the structure of the sequenced ' \
                'PCR amplicons. This is required for correctly parsing the sequencing data.')
                
###############################
##  Sample table parameters  ##
###############################
new_sample_table = Param(
        name = 'new_sample_table',
        value = False,
        type = bool,
        help = 'If True, a sample table is created with the filename '\
                'specified in the "sample_table_file" parameter. The '\
                '"screen_name," "plate_size," "plates_per_lane," '\
                'and "num_lanes" parameters are used to create the '\
                'table; each parameter must be specified if it '\
                'does not possess an acceptable default value.',
        options = [False, True])

screen_name = Param(
        name = 'screen_name',
        value = None,
        type = str,
        help = 'Name of the screen that uniquely identifies '\
                'it from all other performed screens.',
        options = None)

plate_size = Param(
        name = 'plate_size',
        value = 96,
        type = int,
        help = 'Number of conditions per plate.',
        options = None)

plates_per_lane = Param(
        name = 'plates_per_lane',
        value = None,
        type = int,
        help = 'Number of full condition plates per sequencing lane.',
        options = None)

num_lanes = Param(
        name = 'num_lanes',
        value = None,
        type = int,
        help = 'Number of sequencing lanes',
        options = None)

extra_columns = Param(
        name = 'extra_columns',
        value = None,
        type = str,
        help = 'Extra columns to include in the sample information table.'\
                ' Must be separated by commas and not include spaces.',
        options = '_any_')


###########################
##  Basic parameters     ##
###########################

sub_screen_column = Param(
        name = 'sub_screen_column',
        value = None,
        type = str,
        help = 'Column in the sample information table that specifies how '\
                'to partition the raw count data prior to scoring interactions. '\
                'The user must manually add the values into this column that '\
                'identify the different partitions.',
        options = '_any_')

verbosity = Param(
        name = 'verbosity',
        value = 1,
        type = int,
        help = 'Verbosity of messages printed to stdout. '\
                '0 is none, 1 is default (step-relevant announcements), '\
                '2 prints object inspection results, and 3 presents even '\
                'deeper (and more cryptic) inspection results.',
        options = range(4))

###########################
##  Advanced parameters  ##
###########################

num_cores = Param(
        name = 'num_cores',
        value = None,
        type = str,
        help = 'Number of cores to use for processes spread across multiple cores. ' \
                'Defaults to the number of cores detected divided by 2. Can be ' \
                'specified as an integer number of cores (num_cores >= 1), as a fraction of the ' \
                'number of available cores (0 < num_cores < 1), or as "all" to run on all cores. '\
                'Note that "all" may not have much or any speed boost over 0.5.',
        options = '_any_')

remove_barcode_specific_conditions = Param(
        name = 'remove_barcode_specific_conditions',
        value = True,
        type = bool,
        help = 'Remove conditions matching the "barcode-specific" pattern ' \
                '(positive interactions for barcodes beginning with one ' \
                'particular base) with correlation greater than that defined ' \
                'by "barcode_specific_template_correlation_cutoff".',
        options = [False, True])

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
        options = [False, True])

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
        type = float,
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
        type = float,
        help = 'See "condition_pass_read_count".',
        options = None)





