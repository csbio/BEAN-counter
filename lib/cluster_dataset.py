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

def cluster(dataset):

    genes, conditions, matrix = dataset

    # Remove rows and columns that do not have any values - they kill the clustering process!
    good_rows = np.nansum(matrix, axis = 1).astype(np.bool)
    good_cols = np.nansum(matrix, axis = 0).astype(np.bool)

    genes = np.array(genes)[good_rows]
    conditions = np.array(conditions)[good_cols]
    matrix = matrix[np.ix_(good_rows, good_cols)]

    num_genes = len(genes)
    num_conds = len(conditions)

    # Compute distance matrices
    cols_dist = pdist(matrix.transpose(), 'cosine')
    rows_dist = pdist(matrix, 'cosine')

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

