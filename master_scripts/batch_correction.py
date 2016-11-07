#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script takes in a chemical genetic interaction score dataset (matrix and strain/
# condition ids) as well as a sample table containing information on each condition.
# Then, it associates each condition with a "batch" as defined by the entires in the
# specified column in the sample table. Using the batches, the script first performs a
# filtering step in order to remove duplicated conditions within the same batch that we
# would expect to show similarity, as this would create artificial batch effects. Then,
# using only this filtered set of conditions, multiclass LDA is performed to remove
# effects attributable to the batch id. The user specifies the maximum number of LDA
# components to remove.
import pandas as pd
import numpy as np
import scipy
import sys, os, gzip
import jellyfish as jf
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import itertools as it
import cPickle
import argparse

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append(os.path.join(barseq_path, 'lib'))

import config_file_parser as cfp
import compressed_file_opener as cfo
import cg_file_tools as cg_file

sys.path.append(os.path.join(barseq_path, 'lib/python2.7/site-packages'))

def read_sample_table(tab_filename):

    # Read everything in as a string, to prevent vexing
    # number interpretation problems! Methods further down
    # can coerce to different types.
    tab = pd.read_table(tab_filename, dtype = 'S')
    return tab

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

def write_corrected_data_info(dataset_3d, batch_col, nondup_col_list, filename, input_filename):

    f = open(filename, 'wt')
    f.write("- Batch effect correction, using each condition's '{}' value from the sample table\n\n".format(batch_col))
    f.write("- Within each batch, the condtions could not possess duplicate values in these sample table columns: {}\n\n".format(', '.join(nondup_col_list)))
    f.write("- Final dimensions of corrected dataset: {0} stacked matrices, each {1} strains X {2} conditions, representing 0 through {3} LDA components removed\n\n".format(len(dataset_3d[0]), len(dataset_3d[1]), len(dataset_3d[2]), len(dataset_3d) - 1))
    f.write("- Dataset on which batch correction was performed: {}\n".format(input_filename))

def filter_sample_table(sample_table, final_conditions):

    all_conditions = np.array([tuple(row[['screen_name', 'expt_id']]) for i,row in sample_table.iterrows()])
    rows_to_keep = [i for i, val in enumerate(all_conditions) if a_is_row_in_b(val, final_conditions)]
    #print final_conditions
    #print all_conditions
    #print rows_to_keep
    
    return sample_table.iloc[rows_to_keep]

def a_is_row_in_b(a, b):
      
    return np.any(np.all(a == b, axis = 1))

def get_nonreplicating_conditions(sample_table, batch_col, nondup_col_list):

    # Do some checks to make sure the columns are in the sample table
    assert batch_col in sample_table.columns, "Specified batch column '{}' not in the sample table".format(batch_col)
    for col in nondup_col_list:
        assert col in sample_table.columns, "Specified column '{}' to prevent within-batch duplication is not in the sample table".format(col)

    # Create an empty dictionary of dictionaries. The first layer of dictionaries is for the name of the
    # property that should not be duplicated within the same batch (for example, the basic name of the
    # condition, or the chemical ID of the drug, etc). For each of these non-duplicating properties, a
    # dictionary is constructed with one key for each batch. Each of the values for each property --> tag
    # combination is a set containing the actual property values that must not be duplicated. This dictionary
    # is built up and then checked to see if any new property/tag combinations have already occurred.
    nondup_dict = {}
    for col in nondup_col_list:
        nondup_dict[col] = {}
        for tag in set(sample_table[batch_col]):
            nondup_dict[col][tag] = set()

    # Iterate over the rows of the sample table, build up the list of conditions
    # that will be used to get the components for removing batch effects.
    bool_dict = {'True': True, 'TRUE': True, 'False': False, 'FALSE': False}
    inds = []
    for i, row in sample_table.iterrows():
        # First, if the row has been slated to not be included, then don't include! Move on to the next row
        # print row
        if not bool_dict[row['include?']]:
            continue
        # Accept all control conditions, because they should not show strong correlations with batch 
        # identity unless there is a batch effect
        if bool_dict[row['control?']]:
            inds.append(i)
        else:
            accept = True
            tag = row[batch_col]
            for col in nondup_col_list:
                if row[col] in nondup_dict[col][tag]:
                    # If any property (for example, the condition name) has already occurred for a certain
                    # tag, flag that condition as unacceptable to accept as a nonduplicating condition
                    accept = False
            if accept:
                inds.append(i)
                for col in nondup_col_list:
                    nondup_dict[col][tag].add(row[col])
    
    nonreplicating_conditions = np.array([row[['screen_name', 'expt_id']] for i, row in sample_table.iterrows() if i in inds])
    
    return nonreplicating_conditions

def filter_dataset_by_conditions(conditions, matrix, conds_to_keep):
    
    inds_to_keep = np.array([i for i, val in enumerate(conditions) if a_is_row_in_b(val, conds_to_keep)])
    filtered_conditions = conditions[inds_to_keep]
    filtered_matrix = matrix[:, inds_to_keep]

    return [filtered_conditions, filtered_matrix]

def get_batch_classes(dataset, sample_table, batch_col):

    # First, get a mapping of each tag to an integer class
    batches = set(sample_table[batch_col])
    batch_mapping = {}
    for i, batch in enumerate(batches):
        batch_mapping[batch] = i

    # Now, get a vector along the matrix columns indicating which conditions
    # are in which batch
    barcodes, conditions, matrix = dataset
    sample_table = sample_table.set_index(['screen_name', 'expt_id'])
    classes = [batch_mapping[sample_table[batch_col].ix[tuple(cond)]] for cond in conditions]
    return classes

def writeList(l, filename):
    with open(filename, 'wt') as f:
        for x in l:
            f.write(str(x) + '\n')

def get_corrected_data_filename(output_folder):
    return os.path.join(output_folder, 'batch_corrected_datasets.dump.gz')

def get_corrected_data_info_filename(output_folder):
    return os.path.join(output_folder, 'batch_corrected_datasets_info.txt')

def LDA_batch_normalization(dataset, sample_table, batch_col, output_folder, ncomps): # this is actually the batch normalization method
   
    tmp_output_folder = os.path.join(output_folder, 'tmp')

    if not os.path.isdir(tmp_output_folder):
        os.makedirs(tmp_output_folder)
    
    barcodes, filtered_conditions, filtered_matrix, conditions, matrix = dataset
    
    # Remove any remaining NaNs and Infs from the filtered matrix - they would screw
    # up the LDA. 
    filtered_matrix[scipy.isnan(filtered_matrix)] = 0
    filtered_matrix[scipy.isinf(filtered_matrix)] = 0

    # For full matrix, also eliminate NaNs and Infs, BUT preserve the indices and values
    # so they can be added back into the matrix later (not implemented yet, and may never
    # be - there should no longer be NaNs and Infs in the dataset)
    # The NaNs and Infs will mess up the final step of the MATLAB LDA script, which uses
    # matrix multiplication to remove the specified number of components!
    matrix_nan_inds = scipy.isnan(matrix)
    matrix_nan_vals = matrix[matrix_nan_inds]
    matrix_inf_inds = scipy.isinf(matrix)
    matrix_inf_vals = matrix[matrix_inf_inds]

    matrix[matrix_nan_inds] = 0
    matrix[matrix_inf_inds] = 0

    # Save both the small matrix (for determining the components to remove) and the 
    # full matrix for the matlab script
    filtered_matrix_tmp_filename = os.path.join(tmp_output_folder, 'nonreplicating_matrix.txt')
    full_matrix_tmp_filename = os.path.join(tmp_output_folder, 'full_matrix.txt')
    
    np.savetxt(filtered_matrix_tmp_filename, filtered_matrix)
    np.savetxt(full_matrix_tmp_filename, matrix)

    # Map the batch to integers for matlab, and write out to a file so matlab can read
    # Note that yes, the batch_classes should match up with the filtered matrix, not
    # the full matrix
    batch_classes = get_batch_classes(dataset = [barcodes, filtered_conditions, filtered_matrix], sample_table = sample_table, batch_col = batch_col)
    class_tmp_filename = os.path.join(tmp_output_folder, 'classes.txt')
    writeList(batch_classes, class_tmp_filename)
   
    output_tmp_filename = os.path.join(tmp_output_folder, 'full_matrix_lda_normalized.txt')
    runLDAMatlabFunc(filtered_matrix_filename = filtered_matrix_tmp_filename, \
            matrix_filename = full_matrix_tmp_filename, \
            class_filename = class_tmp_filename, \
            ncomps = ncomps, \
            output_filename = output_tmp_filename)
    # The X norm that is returned is the full matrix. In the future, we could add in
    # returning the components to remove so they can be visualized or applied to other
    # one-off datasets
    Xnorm =  scipy.genfromtxt(output_tmp_filename)

    ## Dump the dataset out!
    #output_filename = os.path.join(mtag_effect_folder, 'scaleddeviation_full_mtag_lda_{}.dump.gz'.format(ncomps))
    #of = gzip.open(output_filename, 'wb')
    #cPickle.dump([barcodes, conditions, Xnorm], of)
    #of.close()

    return [barcodes, conditions, Xnorm]

def runLDAMatlabFunc(filtered_matrix_filename, matrix_filename, class_filename, ncomps, output_filename):
    #print "cd /project/csbio/raamesh/projects/smallprojects/coveringYeastGIArray/After330screens/data/Jeff/nameddrugs/; matlab -nodisplay -nodesktop -nojvm -nosplash -r 'multi_class_lda_py('\\'%s\\'', '\\'%s\\'', %f, '\\'%s\\''); exit' > /dev/null;" % (matrixfilename, classfilename, perc, outputfilename)
    barseq_path = os.getenv('BARSEQ_PATH')
    barseq_lib_path = os.path.join(barseq_path, 'lib')
    os.system( "cd %s; matlab -nodisplay -nodesktop -nojvm -nosplash -r 'multi_class_lda_ncomps_full_dataset_py('\\'%s\\'', '\\'%s\\'', '\\'%s\\'', %d, '\\'%s\\''); exit'" % (barseq_lib_path, filtered_matrix_filename, matrix_filename, class_filename, ncomps, output_filename) )


def main(dataset, sample_table, batch_column, nondup_col_list, max_comps, output_folder, input_file, verbosity):

    # sample_table = read_sample_table(sample_table_filename)
    # print sample_table

    # Get the maximum number of LDA components to remove
    # Program will remove 0 through that number of components
    max_comps_remove = int(max_comps)

    # Load in the final zscore dataset
    # barcodes, conditions, matrix = load_dataset(dataset_filename)
    barcodes, conditions, matrix = dataset

    # filter the sample table, so that no samples that have been eliminated
    # are used in further processing
    sample_table = filter_sample_table(sample_table, conditions)
    if verbosity >= 2:
        print sample_table

    # Determine IDs of nonreplicating conditions
    # The user defines which columns of the sample table to use, such that
    # no two conditions with the same batch have the same value in
    # any of the other specified columns
    if nondup_col_list != ['none']:
        # Get the indices of the nonreplicating conditions
        nonrep_cond_inds = get_nonreplicating_conditions(sample_table, batch_column, nondup_col_list)

        # Filter the 3d matrix such that it only contains the nonreplicating conditions
        filtered_conditions, filtered_matrix = filter_3d_dataset_by_conditions(conditions, matrix, nonrep_cond_inds)
        if verbosity >= 2:
            print conditions[0:10]
    else:
        filtered_conditions = conditions
        filtered_matrix = matrix
    
    #nonreplicating_conditions = get_nonreplicating_conditions(sample_table, batch_column, nondup_col_list)
    #if verbosity >= 2:
    #    print nonreplicating_conditions
    #    print nonreplicating_conditions.shape

    # Filter the dataset to include only the nonreplicating conditions
    #filtered_conditions, filtered_matrix = filter_dataset_by_conditions(conditions, matrix, nonreplicating_conditions)

    # Create a list to capture all of the LDA-batch-corrected matrices
    all_mats = []

    # Perform the batch effect normalization (LDA), using a matlab callout
    # Remove 0 components up to and including the max number of comps specified by the user
    final_dataset = [barcodes, filtered_conditions, filtered_matrix, conditions, matrix]
    for i in range(max_comps_remove + 1):
        if verbosity >= 1:
            print "Performing batch effect correction and removing {} components".format(i)
        corrected_dataset = LDA_batch_normalization(final_dataset, sample_table, batch_column, output_folder, i)
        all_mats.append(corrected_dataset[2])

    # Turn the final dataset into a numpy array!
    components_removed = np.array(list(range(max_comps_remove + 1)))
    all_mats = [components_removed, barcodes, conditions, np.array(all_mats)]

    # Write out the 3-d array with different LDA components removed
    corrected_data_filename = get_corrected_data_filename(output_folder)
    dump_dataset(all_mats, corrected_data_filename)

    # Write info on the dumped stacked matrix to a file
    corrected_data_info_filename = get_corrected_data_info_filename(output_folder)
    write_corrected_data_info(all_mats, batch_column, nondup_col_list, corrected_data_info_filename, input_file)


# If this is run from the command line (which is the intention, then run the rest of this script!
if __name__ == '__main__':

    # Get arguments

    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_file', help = 'The dataset on which to perform batch correction.')
    parser.add_argument('sample_table', help = 'The sample table corresponding to the dataset.')
    parser.add_argument('batch_column', help = 'The column from the sample table that defines the batches.')
    parser.add_argument('nondup_columns', help = 'Comma delimited. The columns that contain condition identifiers that should not be duplicated within the same batch.')
    parser.add_argument('max_components', help = 'The maximum number of LDA components to remove.')
    parser.add_argument('output_folder', help = 'The folder to which results are exported.')
    parser.add_argument('-v', '--verbosity', help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')

    args = parser.parse_args()

    # Get the data ready to rumble!
    dataset_file = os.path.abspath(args.dataset_file)
    assert os.path.isfile(dataset_file), "Dataset file does not exist."
    dataset = load_dataset(args.dataset_file)

    #nondup_col_list = args.nondup_columns.split(',')
    # If the user specifies "none", then no columns are used to prevent duplicates of
    # certain values within a particular batch
    if args.nondup_columns is 'none':
        nondup_col_list = ['none']
    else:
        nondup_col_list = args.nondup_columns.split(',')

    assert os.path.isfile(args.sample_table), "Sample table file does not exist."
    sample_table = read_sample_table(args.sample_table)

    assert args.max_components.isdigit(), "The specified maximum number of components to remove is not an integer."
    max_components = int(args.max_components)

    output_folder = os.path.abspath(args.output_folder)
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    if args.verbosity.isdigit():
        verbosity = int(args.verbosity)
    else:
        verbosity = 1

    main(dataset, sample_table, args.batch_column, nondup_col_list, max_components, output_folder, dataset_file, verbosity)
