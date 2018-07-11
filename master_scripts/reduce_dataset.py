#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.6.0'

# This script takes in a dataset and a corresponding sample table with
# a reduced number of rows. This script will reduce the dataset so that
# it contains only the samples corresponding to the rows in the
# dataset. It also exports a sample table that perfectly corresponds
# to the matrix (any extra rows in the sample table have been removed).
import argparse
import numpy as np, pandas as pd
import networkx as nx
import os, sys
import gzip
import cPickle

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'lib'))
from cg_common_functions import read_sample_table, read_barcode_table, bool_dict
from version_printing import update_version_file

#def read_sample_table(tab_filename):
#
#    # Read everything in as a string, to prevent vexing
#    # number interpretation problems! Methods further down
#    # can coerce to different types.
#    tab = pd.read_table(tab_filename, dtype = 'S')
#    return tab

def load_dataset(data_filename):

    f = gzip.open(data_filename, 'rb')
    barcodes, conditions, matrix = cPickle.load(f)
    dataset = [np.array(barcodes), np.array(conditions), matrix]
    f.close()
    return dataset

def dump_dataset(dataset, filename):

    f = gzip.open(filename, 'wb')
    cPickle.dump(dataset, f)
    f.close()

def filter_strain_table(strain_table, final_strains):
    all_strains = np.array(strain_table['Strain_ID'])
    rows_to_keep = [i for i, val in enumerate(all_strains) if val in all_strains]
    return strain_table.iloc[rows_to_keep]

def filter_sample_table(sample_table, final_conditions):

    all_conditions = np.array([tuple(row[['screen_name', 'expt_id']]) for i,row in sample_table.iterrows()])
    rows_to_keep = [i for i, val in enumerate(all_conditions) if a_is_row_in_b(val, final_conditions)]
    #print final_conditions
    #print all_conditions
    #print rows_to_keep

    return sample_table.iloc[rows_to_keep]

def filter_dataset(dataset, info_table, dim):

    strains, conditions, matrix = dataset

    if dim == 'strains':
        reduced_strains = np.array(info_table['Strain_ID'])
        matrix_strains_to_keep = np.array([i for i, val in enumerate(strains) if val in reduced_strains])
        return [strains[matrix_strains_to_keep], conditions, matrix[matrix_strains_to_keep, :]]
    elif dim == 'conditions':
        reduced_conditions = np.array([tuple(row[['screen_name', 'expt_id']]) for i,row in info_table.iterrows()])
        matrix_conds_to_keep = np.array([i for i, val in enumerate(conditions) if a_is_row_in_b(val, reduced_conditions)])
        return [strains, conditions[matrix_conds_to_keep], matrix[:, matrix_conds_to_keep]]
    else:
        assert False, '"dim" argument must be either "strains" or "conditions."'
        

def filter_info_table_on_column(info_table, column, invert):

    keep_bool = [bool_dict[x] for x in info_table[column]]
    if invert:
        keep_bool = np.invert(keep_bool)

    # Return subsetted info table
    return info_table[keep_bool]

#def filter_matrix(sample_table, matrix_conditions, matrix):
#
#    reduced_conditions = np.array([tuple(row[['screen_name', 'expt_id']]) for i,row in sample_table.iterrows()])
#    matrix_conds_to_keep = np.array([i for i, val in enumerate(matrix_conditions) if a_is_row_in_b(val, reduced_conditions)])
#    #print final_conditions
#    #print all_conditions
#    #print rows_to_keep
#
#    return matrix_conditions[matrix_conds_to_keep], matrix[:, matrix_conds_to_keep]

def a_is_row_in_b(a, b):

    return np.any(np.all(a == b, axis = 1))

def main(dataset, info_table, dim, filename, col, inv, verbosity):

    strains, conditions, matrix = dataset

    # First, filter the info table to make sure it doesn't have extra strains/conditions
    if dim == 'strains':
        info_table = filter_strain_table(info_table, strains)
    elif dim == 'conditions':
        info_table = filter_sample_table(info_table, conditions)
    else:
        assert False, '"dim" argument must be either "strains" or "conditions."'

    # Then, reduce the info table based on the column (if specified)
    if col is not None:
        info_table = filter_info_table_on_column(info_table, col, inv)

    # Then, reduce the matrix and condition array
    reduced_dataset = filter_dataset([strains, conditions, matrix], info_table, dim)

    #if verbosity >= 2:
    #    print sample_table
    #    print reduced_conditions.shape
    #    print reduced_matrix.shape

    dump_dataset(reduced_dataset, filename)
    
    update_version_file(os.path.dirname(os.path.abspath(filename)), VERSION)

# If this is run from the command line (which is the intention), then run the rest of this script!
if __name__ == '__main__':

    # Get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_file', help = 'The dataset to be reduced.')
    parser.add_argument('info_table', help = 'The strain/condition info table that contains fewer strains/conditions than are in the dataset.')
    parser.add_argument('output_file', help = 'Filename for the reduced dataset.')
    parser.add_argument('--strains', action = 'store_true', help = 'By default, filtering is performed on the condition side of the dataset (columns). Use this flag to filter the strain side of the dataset')
    parser.add_argument('--column', help = 'As an alternative to using a slimmed-down info table to reduce the dataset, only retain strains/conditions for which the value in this info table column is True')
    parser.add_argument('--invert', action = 'store_true', help = 'If --column is specified, invert the True/False values')    
    parser.add_argument('-v', '--verbosity', default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')

    args = parser.parse_args()

    # Get the data ready to rumble!
    dataset_file = os.path.abspath(args.dataset_file)
    assert os.path.isfile(dataset_file), "Dataset file does not exist."
    dataset = load_dataset(args.dataset_file)

    if args.verbosity >= 2:
        print dataset

    if args.strains:
        dim = 'strains'
    else:
        dim = 'conditions'
    
    assert os.path.isfile(args.info_table), "File {} does not exist.".format(args.info_table)
    if dim == 'strains':
        info_table = read_barcode_table(args.info_table)
    elif dim == 'conditions':
        info_table = read_sample_table(args.info_table)

    if args.verbosity >= 2:
        print info_table
        print dataset


    if args.column is not None:
        assert args.column in info_table.columns, 'Specified --column {} is not a column name in the info table.'.format(args.column)

    # If the filename has directories in it that do not exist, make them
    outdir = os.path.dirname(os.path.abspath(args.output_file))
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    main(dataset, info_table, dim, args.output_file, args.column, args.invert, args.verbosity)

