#!/usr/bin/env python
# This script takes in a chemical genetic interaction score matrix and information
# that maps each condition to their index tags, such that any effects in the data
# that can be explained by the index tag are removed.
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

def load_dataset(data_filename)

    f = gzip.open(data_filename, 'rb')
    barcodes, conditions, matrix = cPickle.load(f)
    dataset = [np.array(barcodes), np.array(conditions), matrix]
    f.close()
    return dataset

def dump_dataset(dataset, filename):

    f = gzip.open(filename, 'wb')
    cPickle.dump(dataset, f)
    f.close()

def write_corrected_data_info(dataset_3d, batch_col, nondup_cols, filename):

    f = open(filename, 'wt')
    f.write("- Batch effect correction, using each condition's '{}' value from the sample table\n\n".format(batch_col)
    f.write("- Within each batch, the condtions could not possess duplicate: {} values\n\n".format(nondup_cols.replace(',', ', ')))
    f.write("- Final dimensions of corrected dataset: {0} stacked matrices, each {1} strains X {2} conditions, representing 0 through {3} LDA components removed\n".format(len(dataset_3d[0]), len(dataset_3d[1]), len(dataset_3d[2]), len(dataset_3d) - 1)

def filter_sample_table(sample_table, final_conditions):

    all_conditions = np.array([tuple(row[['screen_name', 'expt_id']]) for i,row in sample_table.iterrows()])
    rows_to_keep = [i for i, val in enumerate(all_conditions) if a_is_row_in_b(val, final_conditions)]
    #print final_conditions
    #print all_conditions
    #print rows_to_keep
    
    return sample_table.iloc[rows_to_keep]

def a_is_row_in_b(a, b):
      
    return np.any(np.all(a == b, axis = 1))

def get_nonreplicating_conditions(sample_table, batch_col, nondup_cols):

    nondup_col_list = nondup_cols.split(',')

    # Do some checks to make sure the columns are in the sample table
    assert batch_col in sample_table.columns, "Specified batch column '{}' not in the sample table".format(batch_col)
    for col in nondup_col_list:
        assert col in sample_table.columns, "Specified column '{}' to prevent within-batch duplication is not in the sample table".format(col)

    # Create an empty dictionary of dictionaries. The first layer of dictionaries is for the name of the
    # property that should not be duplicated within the same index tag (for example, the basic name of the
    # condition, or the chemical ID of the drug, etc). For each of these non-duplicating properties, a
    # dictionary is constructed with one key for each index tag. Each of the values for each property --> tag
    # combination is a set containing the actual property values that must not be duplicated. This dictionary
    # is built up and then checked to see if any new property/tag combinations have already occurred.
    nondup_dict = {}
    for col in nondup_col_list:
        nondup_dict[col] = {}
        for tag in set(sample_table[batch_col]):
            nondup_dict[col][tag] = set()

    # Iterate over the rows of the sample table, build up the list of conditions
    # that will be used to get the components for removing index tag effects.
    bool_dict = {'True': True, 'TRUE': True, 'False': False, 'FALSE': False}
    inds = []
    for i, row in sample_table.iterrows():
        # First, if the row has been slated to not be included, then don't include! Move on to the next row
        # print row
        if not bool_dict[row['include?']]:
            continue
        # Accept all control conditions, because they should not show strong correlations with index tag
        # idendity unless there is an index tag effect
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

def filter_dataset_nonreplicating(conditions, matrix, conds_to_keep):
    
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
    # are in which index tag class
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

    # Map the index tags to integers for matlab, and write out to a file so matlab can read
    # Note that yes, the index_tag_classes should match up with the filtered matrix, not
    # the full matrix
    batch_classes = get_batch_classes(dataset = [barcodes, filtered_conditions, filtered_matrix], sample_table, batch_col)
    class_tmp_filename = os.path.join(tmp_output_folder, 'classes.txt')
    writeList(index_tag_classes, class_tmp_filename)
   
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


def main(dataset, sample_table, batch_column, nondup_columns, output_folder, max_comps):

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
    print sample_table

    # Determine IDs of nonreplicating conditions
    # The user defines which columns of the sample table to use, such that
    # no two conditions with the same index tag have the same value in
    # any of the other selected columns
    nonreplicating_conditions = get_nonreplicating_conditions(sample_table, batch_column, nondup_columns)
    print nonreplicating_conditions
    print nonreplicating_conditions.shape

    # Filter the dataset to include only the nonreplicating conditions
    filtered_conditions, filtered_matrix = filter_dataset_nonreplicating(conditions, matrix, nonreplicating_conditions)

    # Create a list to capture all of the LDA-batch-corrected matrices
    all_mats = []

    # Perform the index tag effect normalization (LDA), using a matlab callout
    # Remove 0 components up to and including the max number of comps specified by the user
    final_dataset = [barcodes, filtered_conditions, filtered_matrix, conditions, matrix]
    for i in range(max_comps_remove + 1):
        print "Performing index tag normalization and removing {} components".format(i)
        all_mats.append(LDA_batch_normalization(final_dataset, sample_table, batch_column, output_folder, i))

    # Turn the final dataset into a numpy array!
    components_removed = np.array(list(range(max_comps_remove + 1)))
    all_mats = [components_removed, barcodes, conditions, np.array(all_mats)]

    # Write out the 3-d array with different LDA components removed
    corrected_data_filename = get_corrected_data_filename(output_folder)
    dump_dataset(all_mats, corrected_data_filename)

    # Write info on the dumped stacked matrix to a file
    corrected_data_info_filename = get_corrected_data_info_filename(output_folder)
    write_corrected_data_info(all_mats, corrected_data_info_filename)

