#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This file defines some useful functions for computing average
# within- and between-group correlations from a symmetric
# correlation matrix and a vector of groups along the rows and columns
import numpy as np
import itertools as it


def get_group_correlations(cor_mat, group_names, within_group = True):

    # Group names must line up with rows/columns that are in those
    # respective groups
    group_names = np.array(group_names)
    group_names_unique = np.unique(group_names)

    if within_group:
        group_pairs = [(x, x) for x in group_names_unique]
    else:
        # Note that it.combinations, used in this manner, essentially
        # gives preference to the upper triangular of the matrix. An
        # example of this output is [('a', 'b'), ('a', 'c'), ('a', 'd'),
        # ('b', 'c'), ('b', 'd'), ('c', 'd')]
        group_pairs = list(it.combinations(group_names_unique, 2))

    # Get a dictionary that maps each batch to its row and column indices
    group_to_ind = {}
    for i, group in enumerate(group_names):
        if group in group_to_ind:
            group_to_ind[group] = np.append(group_to_ind[group], i)
        else:
            group_to_ind[group] = np.array([i], dtype = np.int)

    #print group_names.shape
    #print group_to_ind

    # Set the lower triangular of the correlation matrix to nan!
    cor_mat[np.tril_indices_from(cor_mat)] = np.nan

    ## Get upper triangular from the first diagonal and above.
    #row_triu_indices, col_triu_indices = np.triu_indices_from(cor_mat, 1)
    ## Convert upper triangular indices into relevant group labels
    #row_triu_group_labels = group_names[row_triu_indices]
    #col_triu_group_labels = group_names[col_triu_indices]

    # Get mean correlation for each within- or between-pair of groups
    # Also keep track of how many correlations contribute to the mean
    mean_cor_array = np.zeros(len(group_pairs))
    n_array = np.zeros(len(group_pairs))
    for i, pair in enumerate(group_pairs):

        ## Get boolean of which row/columns indices match the 
        ## current row/column pair
        #row_group_in_pair = row_triu_group_labels == pair[0]
        #col_group_in_pair = col_triu_group_labels == pair[1]
        #row_col_group_in_pair = row_group_in_pair & col_group_in_pair

        ## Use the boolean from above to get the row and columns
        ## indices for the current pair
        #row_inds = row_triu_indices[row_col_group_in_pair]
        #col_inds = col_triu_indices[row_col_group_in_pair]

        # Get the values
        corrs = cor_mat[np.ix_(group_to_ind[pair[0]], group_to_ind[pair[1]])]
        #print corrs
        #input()
        
        # Add the mean correlation to the list!
        mean_cor_array[i] = np.nanmean(corrs)

        # Add the number of non NaN correlations to the n_list
        n_array[i] = np.sum(np.invert(np.isnan(corrs)))

    return group_pairs, mean_cor_array, n_array





