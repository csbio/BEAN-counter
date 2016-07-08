#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

import scipy
from scipy.spatial.distance import pdist

import numpy as np
from numpy import random

import fastcluster

import Bio
from Bio.Cluster import Node,Tree

import StringIO
import sys
import os

def cluster(dataset, new_matrix = None):

    genes, conditions, matrix = dataset
    print 'shape of original matrix:'
    print matrix.shape

    # Get the indices of NaNs in the matrix
    nan_inds = np.isnan(matrix)
    num_nans = np.sum(nan_inds)

    # Replace the NaN values with extremely small noise values
    # (noise has standard deviation of 1 million times less than the data)
    np.random.seed(96857463)
    data_sd = np.nanstd(matrix)
    noise = np.random.randn(num_nans) * data_sd / 1e6
    matrix[nan_inds] = noise

    ## Remove rows and columns that do not have any values - they kill the clustering process!
    ## print 'row_nansum: {}'.format(np.nansum(matrix, axis = 1))
    ## print 'col_nansum: {}'.format(np.nansum(matrix, axis = 0))
    #good_rows = np.nansum(matrix, axis = 1).astype(np.bool)
    #good_cols = np.nansum(matrix, axis = 0).astype(np.bool)

    #print 'number of good rows: {}'.format(np.nansum(good_rows))
    #print 'number of good cols: {}'.format(np.nansum(good_cols))

    #genes = np.array(genes)[good_rows]
    #conditions = np.array(conditions)[good_cols]
    #matrix = matrix[np.ix_(good_rows, good_cols)]

    #print 'shape of good matrix:'
    #print matrix.shape

    num_genes = len(genes)
    num_conds = len(conditions)

    # Compute distance matrices
    cols_dist = pdist(matrix.transpose(), 'cosine')
    rows_dist = pdist(matrix, 'cosine')

    ## Get the names of rows and columns that have NaN dissimilarity values
    #rows_dist_nan_inds_1, rows_dist_nan_inds_2 = [x[np.isnan(rows_dist)] for x in np.triu_indices(matrix.shape[0], 1)]
    #cols_dist_nan_inds_1, cols_dist_nan_inds_2 = [x[np.isnan(cols_dist)] for x in np.triu_indices(matrix.shape[1], 1)]

    #row_names_nan_dist_1, row_names_nan_dist_2 = genes[rows_dist_nan_inds_1], genes[rows_dist_nan_inds_2]
    #col_names_nan_dist_1, col_names_nan_dist_2 = conditions[cols_dist_nan_inds_1], conditions[cols_dist_nan_inds_2]

    ## And print out the rows (strains) and columns (conditions) with NaN dissimilarity values
    #print "Strain pairs with NaN dissimilarity values:"
    #for i, row_name_1 in enumerate(row_names_nan_dist_1):
    #    print row_name_1, row_names_nan_dist_2[i]
    #print ""
    #
    #print "Condition pairs with NaN dissimilarity values:"
    #for i, col_name_1 in enumerate(col_names_nan_dist_1):
    #    print col_name_1, col_names_nan_dist_2[i]
    #print ""

    

    # Cluster the matrix using fastcluster!
    print 'clustering columns...'
    cols_clust_mat = fastcluster.average(cols_dist)
    print 'clustering rows...'
    rows_clust_mat = fastcluster.average(rows_dist)

    # Transform the values in the clustering matrices so they can be used with Bio.Cluster
    for i in range(num_genes - 1):
        if rows_clust_mat[i, 0] > (num_genes - 1):
            rows_clust_mat[i, 0] = -(rows_clust_mat[i, 0] - (num_genes - 1))
        if rows_clust_mat[i, 1] > (num_genes - 1):
            rows_clust_mat[i, 1] = -(rows_clust_mat[i, 1] - (num_genes - 1))


    for i in range(num_conds - 1):
        if cols_clust_mat[i, 0] > (num_conds - 1):
            cols_clust_mat[i, 0] = -(cols_clust_mat[i, 0] - (num_conds - 1))
        if cols_clust_mat[i, 1] > (num_conds - 1):
            cols_clust_mat[i, 1] = -(cols_clust_mat[i, 1] - (num_conds - 1))

    # Turn into lists of nodes
    cols_nodes_list = [Node(int(cols_clust_mat[i, 0]), int(cols_clust_mat[i, 1]), cols_clust_mat[i, 2]) for i in range(cols_clust_mat.shape[0])]
    rows_nodes_list = [Node(int(rows_clust_mat[i, 0]), int(rows_clust_mat[i, 1]), rows_clust_mat[i, 2]) for i in range(rows_clust_mat.shape[0])]

    # Create trees
    cols_tree = Tree(cols_nodes_list)
    rows_tree = Tree(rows_nodes_list)

    # Add the NaNs back into the matrix, so it can be visualized properly
    matrix[nan_inds] = np.nan

    # If a "new_matrix" was specified, that means that we wanted to use the original dataset to
    # get the clustering but then actually use a different matrix for the data. So, at this point
    # we set the variable "matrix" to be the values of "new_matrix"
    if new_matrix is not None:
        matrix = new_matrix

    # Create a giant text string so that the input data can be turned into a "record" object
    row1 = 'ORF\tNAME\tGWEIGHT\t' + '\t'.join(conditions)
    row2 = 'EWEIGHT\t\t\t' + '\t'.join(['1' for i in range(len(conditions))])
    rows_rest = [['' for i in range(len(conditions) + 3)] for j in range(len(genes))]
    for i in range(len(genes)):
        rows_rest[i][0:2] = [genes[i] for j in range(2)]
        rows_rest[i][2] = '1'
        for j in range(len(conditions)):
            rows_rest[i][j+3] = str(matrix[i, j])
    rows_rest_intermed = ['\t'.join(x) for x in rows_rest]
    rows_rest_final = '\n'.join(rows_rest_intermed)
    final_string = '%s\n%s\n%s' % (row1, row2, rows_rest_final)

    # Read in as a "record" object
    handle = StringIO.StringIO(final_string)
    record = Bio.Cluster.read(handle)

    return record, rows_tree, cols_tree

