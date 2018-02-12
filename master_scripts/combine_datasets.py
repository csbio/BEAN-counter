#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.2.0'

# This script takes in mulitple datasets and their corresponding sample tables and
# spits out a combined dataset and a combined sample table. The only option
# for this script is whether or not strains present in only one of the
# datasets are kept as nans in the combined dataset or if they are excluded
# entirely. The default is to only keep the intersect of the strains across the
# entire dataset. If the user is to perform a batch correction to reduce effects
# specific to each dataset, then one column in the sample table should clearly
# define the dataset that each of the samples came from.
import argparse
import numpy as np, pandas as pd
import networkx as nx
import os, sys
import gzip
import cPickle

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'lib'))
from cg_common_functions import read_sample_table
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

def filter_sample_table(sample_table, final_conditions):

    all_conditions = np.array([tuple(row[['screen_name', 'expt_id']]) for i,row in sample_table.iterrows()])
    rows_to_keep = [i for i, val in enumerate(all_conditions) if a_is_row_in_b(val, final_conditions)]
    #print final_conditions
    #print all_conditions
    #print rows_to_keep

    return sample_table.iloc[rows_to_keep]

def a_is_row_in_b(a, b):

    return np.any(np.all(a == b, axis = 1))

def strain_match(a, b):

    a_inds = []
    b_inds = []
    for i, x in enumerate(a):
        # Changing from matching based on 2-D array to just lists
        # b/c "Strain_ID" is now just a single identifier string.
        #if a_is_row_in_b(x, b):
        if x in b:
            #b_match_ind = np.where(np.all(x == b, axis = 1))[0]
            b_match_ind = np.where(x == np.array(b))[0]
            assert np.size(b_match_ind) == 1, "Strain {} somehow is matching multiple rows in the final strains array".format(x)
            b_match_int = int(b_match_ind)
            a_inds.append(i)
            b_inds.append(b_match_int)

    return [a_inds, b_inds]

# Define a function to get unique rows from a 2-D array (for getting
# unique strains from a combined list). From user545424 on this
# stackoverflow page: http://stackoverflow.com/questions/8560440/removing-duplicate-columns-and-rows-from-a-numpy-2d-array/8567929#8567929
def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def combine_datasets(datasets, all_strains, verbosity):

    # NO MORE STRAIN TUPLE IDENTIFIERS
    # First, get all strains to include in the final matrix - this depends
    # on if the value of all_strains is True (union of all strains) or
    # False (intersect of all strains)
    #strains_tuple_list = [[tuple(z) for z in y] for y in [x[0] for x in datasets]]
    # But I do need to get the set of strains from each dataset
    strains_list = [x[0] for x in datasets]

    if verbosity >= 2:
        print strains_list
        #print strains_tuple_list

    final_strains = set(strains_list[0])
    if all_strains:
        for i in range(1, len(strains_list)):
            final_strains = final_strains.union(set(strains_list[i]))
    else:
        for i in range(1, len(strains_list)):
            final_strains = final_strains.intersection(set(strains_list[i]))
    final_strains = np.array(list(final_strains))

    if verbosity >= 2:
        print final_strains[0:10]
        print len(final_strains)

    # Now, get indices to reorder each matrix's rows to match the new order of the strains
    old_strain_indices = []
    new_strain_indices = []
    for strains in strains_list:
        strain_indices = strain_match(strains, final_strains)
        old_strain_indices.append(np.array(strain_indices[0]))
        new_strain_indices.append(np.array(strain_indices[1]))

    # Create a large matrix that will be filled in, dataset by dataset
    condition_lists = [x[1] for x in datasets]
    final_conditions = np.vstack(condition_lists)
  
    if verbosity >= 2:
        print "shapes of condition arrays before concatenation:"
        for x in condition_lists:
            print x.shape
        print "shape of concatenated condition array: {}".format(final_conditions.shape)
   
    cond_lengths = np.array([x.shape[0] for x in condition_lists])
    tot_conds = np.sum(cond_lengths)
    end_cols = np.cumsum(cond_lengths)
    start_cols = end_cols - cond_lengths
    final_matrix = np.zeros((final_strains.shape[0], tot_conds))
    final_matrix.fill(np.nan)

    # Fill in the final matrix!
    if verbosity >= 2:
        print 'final_matrix shape: {}'.format(final_matrix.shape)
    for i, dset in enumerate(datasets):
        if verbosity >= 2:
            print new_strain_indices[i].shape
            print start_cols[i]
            print end_cols[i]
            print dset[2].shape
            print old_strain_indices[i].shape
            print np.max(new_strain_indices[i])
        final_matrix[new_strain_indices[i], start_cols[i]:end_cols[i]] = dset[2][old_strain_indices[i]]

    return [final_strains, final_conditions, final_matrix]

def combine_sample_tables(tables):

    # Iterate over all of the sample tables and get the column names
    colname_list = [list(x.columns) for x in tables]

    # This should be pretty simple using the append method for data frames...
    combined_table = tables[0].append(tables[1:]).reset_index(drop = True)

    return combined_table

def main(dataset_list, sample_table_list, all_strains, output_folder, verbosity):

    combined_dataset = combine_datasets(dataset_list, all_strains, verbosity)

    if verbosity >= 2:
        print combined_dataset
        print combined_dataset[0].shape
        print combined_dataset[1].shape
        print combined_dataset[2].shape

    combined_sample_table = combine_sample_tables(sample_table_list)

    if verbosity >= 2:
        print combined_sample_table

    dataset_filename = os.path.join(output_folder, 'combined_dataset.dump.gz')
    table_filename = os.path.join(output_folder, 'combined_sample_table.txt')

    dump_dataset(combined_dataset, dataset_filename)
    combined_sample_table.to_csv(table_filename, sep = '\t', header = True, index = False)

    update_version_file(output_folder, VERSION)

# If this is run from the command line (which is the intention), then run the rest of this script!
if __name__ == '__main__':

    # Get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_and_sample_table_files', metavar = 'DATASET_1 SAMPLE_TABLE_1 DATASET_2 SAMPLE_TABLE_2 ...', nargs = '+', help = 'The datasets and sample tables that are to be combined, in alternating order.')
    parser.add_argument('--output_folder', help = 'The folder to which the resulting combined matrix and sample table are written.')
    parser.add_argument('--all_strains', action = 'store_true', help = 'Add this flag if you want to keep all of the strains in the combined dataset, instead of the interesect of the strains. Strains present in one dataset but not the other will be represented as NaNs in the combined dataset.')
    parser.add_argument('-v', '--verbosity', type = int, default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')

    args = parser.parse_args()

    # Get the data ready to rumble!
    dataset_list = []
    sample_table_list = []
    for i, fname in enumerate(args.dataset_and_sample_table_files):
        fname = os.path.abspath(fname)
        assert os.path.isfile(fname), "File {} does not exist.".format(fname)
        if i % 2 == 0:
            assert fname.endswith('dump.gz'), "File {} is not a gzipped dump file".format(fname)
            f = gzip.open(fname)
            dataset = cPickle.load(f)
            dataset_list.append(dataset)
            f.close()
        elif i % 2 == 1:
            assert fname.endswith('.txt'), "File {} is not a text file".format(fname)
            tab = read_sample_table(fname)
            sample_table_list.append(tab)
        else:
            assert False, "something is very wrong if 'i % 2' returns something other than 0 or 1."

    all_strains = args.all_strains

    if args.verbosity >= 2:
        print dataset_list
        print sample_table_list
        print args.all_strains

    output_folder = args.output_folder
    assert output_folder is not None, "User must specify a value for the --output_folder option"
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    main(dataset_list, sample_table_list, all_strains, output_folder, args.verbosity)

