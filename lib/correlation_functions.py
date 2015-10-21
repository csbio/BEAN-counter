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
        group_pairs = list(it.combinations(group_names_unique, 2))

    # Get upper triangular from the first diagonal and above.
    row_triu_indices, col_triu_indices = np.triu_indices_from(cor_mat, 1)
    # Convert upper triangular indices into relevant group labels
    row_triu_group_labels = group_names[row_triu_indices]
    col_triu_group_labels = group_names[col_triu_indices]

    # Get mean correlation for each within- or between-pair of groups
    # Also keep track of how many correlations contribute to the mean
    mean_cor_list = []
    n_list = []
    for pair in group_pairs:

        # Get boolean of which row/columns indices match the 
        # current row/column pair
        row_group_in_pair = row_triu_group_labels == pair[0]
        col_group_in_pair = col_triu_group_labels == pair[1]
        row_col_group_in_pair = row_group_in_pair & col_group_in_pair

        # Use the boolean from above to get the row and columns
        # indices for the current pair
        row_inds = row_triu_indices[row_col_group_in_pair]
        col_inds = col_triu_indices[row_col_group_in_pair]

        # Append the mean correlation to the list!
        mean_cor_list.append(np.nanmean(cor_mat[row_inds, col_inds]))

        # Append the number of non NaN correlations to the n_list
        n_list.append(np.sum(np.invert(np.isnan(cor_mat[row_inds, col_inds]))))

    return group_pairs, np.array(mean_cor_list), np.array(n_list)





