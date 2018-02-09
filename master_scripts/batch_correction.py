#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.2.0'

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

import compressed_file_opener as cfo
import cg_file_tools as cg_file
from cg_common_functions import read_sample_table
from version_printing import update_version_file

sys.path.append(os.path.join(barseq_path, 'lib/python2.7/site-packages'))

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import datasets

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
	
# Performs LDA using scikit. Returns list with length equal to the number of
# unique classes plus 1 containing matrices with components [0,1...n_classes]
# removed
def scikit_lda(small_x, full_x, classes, n_comps):

    # Fits the LDA
    lda = LinearDiscriminantAnalysis(n_components=n_comps)
    lda.fit(small_x.transpose(), classes)

    # Removes up to n_comps components one by one
    mats = [full_x]
    full_x = full_x.transpose()
    for i in range(1, n_comps+1):
        temp_scalings = lda.scalings_[:, 0:i]
        norm_mat = full_x - np.matmul(full_x, np.matmul(temp_scalings, temp_scalings.transpose()))
        mats.append(norm_mat.transpose())
    #mats.append(np.matmul(full_x, np.matmul(lda.scalings_, lda.scalings_.transpose())))
    return mats
	
def inner_python_lda(small_x, full_x, classes, n_comps):

    # Gets boolean index variables
    classes = np.asarray(classes)
    X = np.transpose(small_x)
    a = np.sum(np.absolute(small_x), axis=0)
    b = np.sum(np.absolute(small_x), axis=1)
    X = X[np.ix_(a > 0, b > 0)]

    # Indexes classes and initializes scatter vectors
    classes = classes[a > 0]
    class_labels = np.unique(classes)
    dataMean = np.nanmean(X, 0)     # Column means
    Sw = np.zeros((X.shape[1], X.shape[1]))       # Empty n*n within-class scatter vector
    Sb = np.zeros((X.shape[1], X.shape[1]))       # Empty n*n between-class scatter vector

    # Computes scatter vectors for all classes
    for i in range(len(class_labels)):
        ind = np.where(classes == class_labels[i])[0]  # Finds all columns for a given class
        if ind.size == 0:
            continue
        # Just swapped to X[:, ind] from X[ind, :], since the latter threw an error
        classMean = np.nanmean(X[ind, :])   # Take mean of all columns in the class AND ALSO CHECK THIS OUT LATER!!!
        Sw = Sw + np.cov(X[ind, :], bias=1, rowvar=0)      # Find covariance of all columns in the class
        Sb = Sb + ind.size*np.transpose(classMean - dataMean)*(classMean - dataMean)    # CHECK THIS OUT LATER!!!
        # FOR REAL!!!

    # Gets matrix to decompose, and decomposes it
    eig_mat = np.linalg.pinv(Sw)*Sb
    U, D, V = np.linalg.svd(eig_mat)
    #a = np.diag(D)/max(np.diag(D))
    stopind = n_comps

    N = V[:, 0:stopind]
    Xnorm = np.transpose(full_x)
    a = np.sum(np.absolute(Xnorm), axis=0)
    b = np.sum(np.absolute(Xnorm), axis=1)

    Xnorm[np.ix_(b > 0, a > 0)] = Xnorm[np.ix_(b > 0, a > 0)] - np.matmul(Xnorm[np.ix_(b > 0, a > 0)], np.matmul(N, np.transpose(N)))
    Xnorm = Xnorm - np.matmul(Xnorm, np.matmul(N, np.transpose(N)))
    Xnorm = np.transpose(Xnorm)

    return Xnorm
	
def outer_python_lda(small_x, full_x, classes, n_comps):
    mats = [full_x]
    for i in range(1, n_comps+1):
        mats.append(inner_python_lda(small_x, full_x, classes, n_comps))
    return mats

def LDA_batch_normalization(dataset, sample_table, batch_col, output_folder, n_comps): # this is actually the batch normalization method
   
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

    # Map batch classes to integers
    batch_classes = get_batch_classes(dataset = [barcodes, filtered_conditions, filtered_matrix], sample_table = sample_table, batch_col = batch_col)
	
    # Checks number of classes and limits ncomps
    a = [x > 0 for x in np.sum(np.absolute(filtered_matrix), axis=0)]
    classes = np.asarray([batch_classes[i] for i in range(len(batch_classes)) if a[i]])
    n_samples = filtered_matrix.shape[0]
    n_classes = len(np.unique(classes))
    if n_samples == n_classes:
        print "ERROR: The number of samples is equal to the number of classes. Exiting"
    if n_classes <= n_comps:
        print "Fewer classes, " + str(n_classes) + ", than components. Setting components to " + str(n_classes-1)
        n_comps = n_classes-1

    # Runs LDA
    #Xnorm = scikit_lda(filtered_matrix, matrix, batch_classes, n_comps)
    Xnorm = outer_python_lda(filtered_matrix, matrix, batch_classes, n_comps)

    return [barcodes, conditions, Xnorm, n_comps]

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
    nonreplicating_conditions = get_nonreplicating_conditions(sample_table, batch_column, nondup_col_list)
    if verbosity >= 2:
        print nonreplicating_conditions
        print nonreplicating_conditions.shape

    # Filter the dataset to include only the nonreplicating conditions
    filtered_conditions, filtered_matrix = filter_dataset_by_conditions(conditions, matrix, nonreplicating_conditions)

    # Create a list to capture all of the LDA-batch-corrected matrices
    all_mats = []

	# Scikit LDA implementation
    final_dataset = [barcodes, filtered_conditions, filtered_matrix, conditions, matrix]
    corrected_dataset = LDA_batch_normalization(final_dataset, sample_table, batch_column, output_folder, max_comps_remove)
    all_mats = corrected_dataset[2]
    n_comps = corrected_dataset[3]

    # Turn the final dataset into a numpy array!
    components_removed = np.array(list(range(max_comps_remove + 1)))
    all_mats = [components_removed, barcodes, conditions, np.array(all_mats)]

    # Write out the 3-d array with different LDA components removed
    corrected_data_filename = get_corrected_data_filename(output_folder)
    dump_dataset(all_mats, corrected_data_filename)

    # Write info on the dumped stacked matrix to a file
    corrected_data_info_filename = get_corrected_data_info_filename(output_folder)
    write_corrected_data_info(all_mats, batch_column, nondup_col_list, corrected_data_info_filename, input_file)

    update_version_file(output_folder, VERSION)


# If this is run from the command line (which is the intention, then run the rest of this script!
if __name__ == '__main__':

    # Get arguments

    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_file', help = 'The dataset on which to perform batch correction.')
    parser.add_argument('sample_table', help = 'The sample table corresponding to the dataset.')
    parser.add_argument('batch_column', help = 'The column from the sample table that defines the batches.')
    parser.add_argument('nondup_columns', help = 'Comma delimited. The columns that contain condition identifiers that should not be duplicated within the same batch.')
    parser.add_argument('max_components', type = int, help = 'The maximum number of LDA components to remove.')
    parser.add_argument('output_folder', help = 'The folder to which results are exported.')
    parser.add_argument('-v', '--verbosity', type = int, default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')

    args = parser.parse_args()

    # Get the data ready to rumble!
    dataset_file = os.path.abspath(args.dataset_file)
    assert os.path.isfile(dataset_file), "Dataset file does not exist."
    dataset = load_dataset(args.dataset_file)

    nondup_col_list = args.nondup_columns.split(',')

    assert os.path.isfile(args.sample_table), "Sample table file does not exist."
    sample_table = read_sample_table(args.sample_table)

    output_folder = os.path.abspath(args.output_folder)
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    main(dataset, sample_table, args.batch_column, nondup_col_list, args.max_components, output_folder, dataset_file, args.verbosity)
