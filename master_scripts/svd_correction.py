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

def get_include_col_conditions(sample_table, include_col):

    # Allow a leading exclamation point to negate a column. Default is false
    negate = False
    if include_col.startswith('!'):
        negate = True
        include_col = include_col[1:]

    # Do some checks to make sure the columns are in the sample table
    assert include_col in sample_table.columns, "Specified include column '{}' not in the sample table".format(include_col)

    bool_dict = {'True': True, 'TRUE': True, 'False': False, 'FALSE': False}
    include_vals = [bool_dict[x] for x in sample_table[include_col]]
    if negate:
        include_vals = np.invert(include_vals)
    conditions = np.array(sample_table.loc[include_vals, ['screen_name', 'expt_id']])
    return conditions

def filter_dataset_by_conditions(conditions, matrix, conds_to_keep):
    
    inds_to_keep = np.array([i for i, val in enumerate(conditions) if a_is_row_in_b(val, conds_to_keep)])
    filtered_conditions = conditions[inds_to_keep]
    filtered_matrix = matrix[:, inds_to_keep]

    return [filtered_conditions, filtered_matrix]

def writeList(l, filename):
    with open(filename, 'wt') as f:
        for x in l:
            f.write(str(x) + '\n')

def get_corrected_data_filename(output_folder):
    return os.path.join(output_folder, 'svd_corrected_datasets.dump.gz')

def get_svd_components_filename(output_folder):
    return os.path.join(output_folder, 'svd_removed_components.dump.gz')

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

def svd_correction(matrix, filtered_matrix, n, verbosity):
    '''
    Computes the SVD on a subset of matrix (filtered_matrix),
    and successively removes SVD components from the full matrix.

    Returns a stack of matrices, each of which has from 0 to n
    SVD components removed
    '''

    # Get the SVD of the filtered matrix
    U, s, V = np.linalg.svd(filtered_matrix, full_matrices = False)

    # Loop through the individual components to remove and ensure
    # that the singular vectors have the right sign. (there is no
    # guarantee that the singular vectors have the correct sign)
    # the correct sign is ensured by creating a vector of "multipliers"
    # that are either 1 or -1.
    sign_multipliers = [1]
    # This loop goes from 0 to n-1, because that is how the components are
    # indexed in the numpy arrays
    for i in range(n):
        comp_mat_UsV = get_one_component_via_UsV(U, s, V, i, verbosity)
        comp_mat_proj = get_one_component_via_proj(U, filtered_matrix, i, verbosity)
        if np.allclose(comp_mat_UsV, comp_mat_proj):
            sign_multipliers.append(1)
        elif np.allclose(comp_mat_UsV, -comp_mat_proj):
            sign_multipliers.append(-1)
        else:
            assert False, "SVD not computed properly"

    # Set up the list that will contain the corrected matrices
    corrected_mats = []
    sign_multipliers = np.array(sign_multipliers)
    
    # Create a list of each component that is successively removed from the matrix
    removed_component_matrices = []

    # Successively remove component and put the results in
    # the corrected mats list
    components_removed = []
    for i in range(n + 1):
        if verbosity >= 2:
            print 'Removing {} components...'.format(i)
        components_removed.append(i)
        multipliers = sign_multipliers[0:i]
        corrected_mats.append(remove_n_components(matrix, U, i, multipliers, verbosity))
        removed_component_matrices.append(get_n_plus_oneth_component(matrix, U, i, sign_multipliers[i], verbosity))

    return [np.array(components_removed), np.array(corrected_mats), np.array(removed_component_matrices)]
       
def get_one_component_via_UsV(U, s, V, n, verbosity):

    s_copy = np.copy(s)
    s_copy[np.array(range(len(s))) != n] = 0
    S = np.zeros((U.shape[1], V.shape[0]))
    S[:len(s_copy), :len(s_copy)] = np.diag(s_copy)
    if verbosity >= 2:
        print "Singular values, vector form:"
        print s
        print s.shape
        print "Singular values, matrix form:"
        print S
        print S.shape
        print "Eigen-condition matrix:"
        print U
        print U.shape
        print "Eigen-strain matrix:"
        print V
        print V.shape

    return np.dot(U, np.dot(S, V))

def get_one_component_via_proj(U, mat, n, verbosity):

    component_n = U[:, n]
    contribution = np.dot(component_n, mat)
    component_n_col = np.reshape(component_n, (-1, 1))
    contribution_row = np.reshape(contribution, (1, -1))
    if verbosity >= 2:
        print "Eigenvector to remove:"
        print component_n_col
        print component_n_col.shape
        print "Contribution across all conditions:"
        print contribution_row
        print contribution_row.shape

    return np.multiply(component_n_col, contribution_row)

def remove_n_components(mat, U, n, mults, verbosity):

    if n == 0:
        return mat
    else:
        components_n = U[:, 0:n]
        components_n_cols = np.reshape(components_n, (-1, n))
        if verbosity >= 2:
            print components_n
            print components_n.shape
            print components_n_cols
            print components_n_cols.shape
        # correct the signs (should test if this actually matters...oh well it is good and safe)
        components_n_cols = mults * components_n_cols
        if verbosity >= 2:
            print components_n_cols
            print components_n_cols.shape
        
        contribution = np.dot(components_n_cols.T, mat)
        contribution_rows = np.reshape(contribution, (n, -1))
        if verbosity >= 2:
            print contribution
            print contribution.shape
            print contribution_rows
            print contribution_rows.shape
        
        mat_to_remove = np.dot(components_n_cols, contribution_rows)
        
        return mat - mat_to_remove

def get_n_plus_oneth_component(mat, U, n, mults, verbosity):

    components_n = U[:, n]
    components_n_cols = np.reshape(components_n, (-1, 1))
    if verbosity >= 2:
        print components_n
        print components_n.shape
        print components_n_cols
        print components_n_cols.shape
    # correct the signs (should test if this actually matters...oh well it is good and safe)
    components_n_cols = mults * components_n_cols
    if verbosity >= 2:
        print components_n_cols
        print components_n_cols.shape
    
    contribution = np.dot(components_n_cols.T, mat)
    contribution_rows = np.reshape(contribution, (1, -1))
    if verbosity >= 2:
        print contribution
        print contribution.shape
        print contribution_rows
        print contribution_rows.shape
    
    mat_to_remove = np.dot(components_n_cols, contribution_rows)
    
    return mat_to_remove

def main(dataset, sample_table, max_comps, include_column, output_folder, input_file, verbosity):

    barcodes, conditions, matrix = dataset
    # If an include column was specified:
    # include only conditions that were specified in the optional "include_column" argument
    if include_column is not None:
        include_cond = get_include_col_conditions(sample_table, include_column)
        if verbosity >= 2:
            print include_cond[0:10]
            print conditions[0:10]
        sample_table = filter_sample_table(sample_table, include_cond)
        filtered_conditions, filtered_matrix = filter_dataset_by_conditions(conditions, matrix, include_cond)
    else:
        filtered_matrix = np.copy(matrix)


    # Get the maximum number of LDA components to remove
    # Program will remove 0 through that number of components
    max_comps_remove = int(max_comps)

    # Load in the final zscore dataset
    barcodes, conditions, matrix = dataset

    # filter the sample table, so that no samples that have been eliminated
    # are used in further processing
    sample_table = filter_sample_table(sample_table, conditions)
    if verbosity >= 2:
        print sample_table
    
    # Do some checking of things:
    if verbosity >= 2:
        print 'number of conditions: {}'.format(len(conditions))
        print 'matrix shape: {}'.format(matrix.shape)
        print 'filtered_matrix shape: {}'.format(filtered_matrix.shape)

    # Perform the SVD correction
    # Remove 0 components up to and including the max number of comps specified by the user
    if verbosity >= 1:
        print "Performing batch effect correction and removing up to and including {} components".format(max_comps_remove)
    components_removed, all_svd_matrices, removed_component_matrices = svd_correction(matrix, filtered_matrix, max_comps, verbosity)

    # Put all the pieces of data together
    corrected_3d_dataset = [components_removed, barcodes, conditions, all_svd_matrices]
    components_3d_dataset = [components_removed, barcodes, conditions, removed_component_matrices]

    # Write out the 3-d array with different LDA components removed
    corrected_data_filename = get_corrected_data_filename(output_folder)
    svd_components_filename = get_svd_components_filename(output_folder)
    dump_dataset(corrected_3d_dataset, corrected_data_filename)
    dump_dataset(components_3d_dataset, svd_components_filename)

    # Don't think I need to write SVD info to file...there's not much to write except which samples were used,
    # which is in the sample info table and the script!
    # Write info on the dumped stacked matrix to a file
    # corrected_data_info_filename = get_corrected_data_info_filename(output_folder)
    # write_corrected_data_info(corrected_3d_dataset, corrected_data_info_filename, input_file)


# If this is run from the command line (which is the intention, then run the rest of this script!
if __name__ == '__main__':

    # Get arguments

    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_file', help = 'The dataset on which to perform batch correction.')
    parser.add_argument('sample_table', help = 'The sample table corresponding to the dataset.')
    parser.add_argument('max_components', help = 'The maximum number of SVD components to remove.')
    parser.add_argument('output_folder', help = 'The folder to which results are exported.')
    parser.add_argument('-incl', '--include_column', help = 'Column in the sample table that specifies True/False whether or not the condition in that row should be used when computing the SVD components to remove. If it starts with an exclamation point, then the values in the column are negated first.')
    parser.add_argument('-v', '--verbosity', help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')

    args = parser.parse_args()

    # Get the data ready to rumble!
    dataset_file = os.path.abspath(args.dataset_file)
    assert os.path.isfile(dataset_file), "Dataset file does not exist."
    dataset = load_dataset(args.dataset_file)

    assert os.path.isfile(args.sample_table), "Sample table file does not exist."
    sample_table = read_sample_table(args.sample_table)

    assert args.max_components.isdigit(), "The specified maximum number of components to remove is not an integer."
    max_components = int(args.max_components)

    output_folder = os.path.abspath(args.output_folder)
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    # Take care of the include column optional argument
    if args.include_column is not None:
        include_column = args.include_column
    else:
        include_column = None

    if args.verbosity.isdigit():
        verbosity = int(args.verbosity)
    else:
        verbosity = 1

    main(dataset, sample_table, max_components, include_column, output_folder, dataset_file, verbosity)
