#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.5.0'

# This script takes in a dataset, a correlation cutoff, and a sample table with a couple of important
# columns to be specified by the user. These columns tell the user the groups that should be
# collapsed (or try to collapse them at least). The collapsing algorithm does the following:
# For each "collapse-by" group:
#   - Extract the profiles and compute a similarity matrix
#   - Binarize this matrix using the specified cutoff
#   - Convert to a graph and compute the connected components
#   - Choose the best connected component based on the average similarity
#   - Now that the best connected component has been chosen, grab the profiles and average them
#   - Also, grab the corresponding rows from the sample table and collapse the values in a
#       pre-specified manner (if the number is a float, take mean and sd, if it is a string
#       that is the same across all replicates, then collapse to one instance of that string,
#       for everything else, semicolon-delimit).
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

def get_groups_and_maps(sample_table, collapse_col, conditions, verbosity):

    # Generate a np array of the unique groups given by the "collapse_col"
    uniq_groups = np.unique(sample_table[collapse_col])

    # Generate a map from each condition to its column index in the matrix
    cond_to_index = {tuple(val):i for i, val in enumerate(conditions)}
    
    # Generate a map from each replicate group to a numpy array of condition names, and another one that maps to
    # matrix indices!
    group_to_cond_map = dict()
    group_to_ind_map = dict()
    tab_groups = sample_table.groupby(collapse_col)
    for key, tab in tab_groups:
        if verbosity >= 2:
            print key
            print tab
        group_to_cond_map[key] = np.array([(row[1]['screen_name'], row[1]['expt_id']) for row in tab.iterrows()])
        group_to_ind_map[key] = np.array([cond_to_index[tuple(x)] for x in group_to_cond_map[key]])
    
    return(uniq_groups, group_to_cond_map, group_to_ind_map)

def collapse_one_group(group_name, matrix, name_to_col_index, name_to_sample_id, sample_table, cor_cutoff, cols_to_keep, how_to_collapse, verbosity):

    # Get the matrix and conditions that correspond to each group
    group_mat = matrix[:, name_to_col_index[group_name]]
    group_conds = name_to_sample_id[group_name]
   
    if verbosity >= 2:
        print group_mat
        print group_conds

    if len(group_conds) < 2:
        return None, None

    # Compute the similarity matrix for the group of replicates
    group_sim_mat = np.corrcoef(group_mat, rowvar = 0)
    if verbosity >= 2:
        print group_sim_mat

    # Convert the group_sim_mat to a binary matrix so the replicates can be split out
    group_sim_mat_bin = group_sim_mat >= cor_cutoff

    # Get the connected components in this matrix
    rep_graph = nx.from_numpy_matrix(group_sim_mat_bin)
    conn_comps = np.array([list(x) for x in nx.connected_components(rep_graph)])
    conn_comp_sizes = [len(x) for x in conn_comps]
    
    # Iterate over the connected components and compute a mean
    # similarity for each. Exclude all components not larger
    # than one condition.
    comp_mean_sims = []
    for i, comp in enumerate(conn_comps):
        if conn_comp_sizes[i] < 2:
            continue
        if verbosity >= 2:
            print comp
        comp_sim_mat = group_sim_mat[np.ix_(comp, comp)]
        if verbosity >= 2:
            print comp_sim_mat
        comp_mean_sims.append(np.mean(comp_sim_mat[np.triu_indices_from(comp_sim_mat, k =1)]))

    # If no connected components of size two or greater are extracted,
    # then treat this as a failure to replicate and return None
    if len(comp_mean_sims) < 1:
        return None, None

    print 'test_1'
    comp_mean_sims = np.array(comp_mean_sims)

    # Sort the connected components and choose the one with the highest correlation
    comp_mean_sim_order = np.argsort(-comp_mean_sims)
    #if len(comp_mean_sim_order) < 2:
    #    best_comp = conn_comps[0]
    #else:
    #    best_comp = conn_comps[comp_mean_sim_order, :]
    best_comp = conn_comps[comp_mean_sim_order][0]
    best_comp_sim = comp_mean_sims[comp_mean_sim_order][0]
    best_comp_conds = group_conds[best_comp]
    if verbosity >= 2:
        print conn_comps
        print comp_mean_sim_order
        print best_comp
        print comp_mean_sims
        print best_comp_sim
        print group_conds
        print best_comp_conds
    best_comp_avg_profile = np.nanmean(group_mat[:, best_comp], axis = 1)

    print 'test_2'
    # Grab the desired rows from the sample table and collapse the columns as
    # specified by the user
    if verbosity >= 2:
        print best_comp_conds
        print sample_table
    best_comp_conds_tup = [tuple(x) for x in best_comp_conds]
    # Note: sample table has already been indexed outside of this function
    group_tab = sample_table.ix[best_comp_conds_tup]
    print 'test_3'
    
    # First, create the new table row and add the screen name and an empty string for
    # the new ID of the collapsed replicate group.
    new_screen_name_array = np.unique(best_comp_conds[:, 0])
    # Originally, I had this check below to prevent replicates from being collapsed across screens...
    # HOWEVER, this now seems silly and I am commenting this out. If the user wants to keep screens
    # separate, he/she can add the screen name to the column used for determining the replicate groups!
    # assert len(new_screen_name_array) < 2, 'Multiple screen names found in one replicate group;\nthis is not allowed'
    if len(new_screen_name_array) == 0:
        # Given the other checks I put in, this code should never be executed
        return None
    new_tab = pd.Series({'screen_name': str(new_screen_name_array[0])})
  
    # Start the column order array so that the series is ordered correctly when returned
    col_order = ['screen_name', 'expt_id', 'mean_correlation']
    
    # Add in the average replicate correlation, just for fun!
    new_tab['mean_correlation'] = best_comp_sim
    

    for i, col in enumerate(cols_to_keep):
        collapse = how_to_collapse[i]
        if verbosity >= 2:
            print col
            print collapse
        # if the option is 'uniq', then attempt to get one unique value. If not successful,
        # just concatenate as before.
        if collapse == 'uniq':
            val_array = np.unique(group_tab[col])
            if len(val_array) > 1:
                new_tab[col] = ';'.join(np.array(group_tab[col], dtype = np.str))
            else:
                new_tab[col] = str(val_array[0])
            col_order.append(col)
        # if the option is mean, then compute mean and standard deviation
        elif collapse == 'mean':
            vals = np.array(group_tab[col], dtype = np.float)
            mean = np.nanmean(vals)
            stdev = np.nanstd(vals)
            new_tab[col + '_mean'] = mean
            new_tab[col + '_stdev'] = stdev
            col_order.append(col + '_mean')
            col_order.append(col + '_stdev')
        # if the option given is 'concat' or any other not legitimate option,
        # then just join them together with semicolons
        elif collapse == 'concat':
            new_tab[col] = ';'.join(np.array(group_tab[col], dtype = np.str))
            col_order.append(col)
        else:
            new_tab[col] = ';'.join(np.array(group_tab[col], dtype = np.str))
            col_order.append(col)
    
    # Reorder the series!
    new_tab = new_tab.reindex(index = col_order)

    if verbosity >= 2:
        print new_tab
        print best_comp_avg_profile

    return new_tab, best_comp_avg_profile
    
    

def main(dataset, sample_table, collapse_col, cor_cutoff, output_folder, cols_to_keep, how_to_collapse, verbosity):

    barcodes, conditions, matrix = dataset

    # Set all NaNs in the matrix to zero. While great pains are undertaken to
    # make sure no data points are NaN, this can still happen, especially when
    # combining datasets.
    matrix[np.isnan(matrix)] = 0
    
    # Filter the sample table, so that no samples that have been eliminated
    # are used in further processing
    sample_table = filter_sample_table(sample_table, conditions)

    unique_groups, group_to_cond, group_to_ind = get_groups_and_maps(sample_table, collapse_col, conditions, verbosity)

    if verbosity >= 2:
        print unique_groups[0:10]
        print len(unique_groups)
        print group_to_cond.items()[0:10]
        print len(group_to_cond)
        print group_to_ind.items()[0:10]
        print len(group_to_ind)

    # Create containers for the new collapsed profiles and sample table rows
    sample_table['individual_rep_ids'] = ['{}_{}'.format(x[1]['screen_name'], x[1]['expt_id']) for x in sample_table.iterrows()]
    sample_table = sample_table.set_index(['screen_name', 'expt_id'])
    collapsed_profile_list = []
    sample_tab_row_list = []
    
    # Here I make sure that the column that is used to determine the groups of samples to
    # collapse is retained in the exported sample table.
    if cols_to_keep is None:
        cols_to_keep = []
        how_to_collapse = []
    cols_to_keep.append(collapse_col)
    how_to_collapse.append('uniq')
    # This allows automatic tracking of which replicates contributed to the collapsed profile
    cols_to_keep.append('individual_rep_ids')
    how_to_collapse.append('concat')
    
    for group in unique_groups:
        if verbosity >= 2:
            print group
        sample_tab_row, collapsed_profile = collapse_one_group(group, matrix, group_to_ind, group_to_cond, sample_table, cor_cutoff, cols_to_keep, how_to_collapse, verbosity)
        if collapsed_profile is None:
            continue
        collapsed_profile_list.append(collapsed_profile)
        sample_tab_row_list.append(sample_tab_row)
    assert len(collapsed_profile_list) == len(sample_tab_row_list), "The number of rows in the collapsed sample table and the number of columns in the collapsed matrix are not the same!"

    if verbosity >= 2:
        print collapsed_profile_list[0:3]
        print sample_tab_row_list[0:3]
        print np.unique([len(x) for x in sample_tab_row_list])
        #print pd.DataFrame

    # Combine profiles into matrix, rows into one table
    collapsed_matrix = np.vstack(collapsed_profile_list).T
    collapsed_sample_table = pd.DataFrame(sample_tab_row_list)
    
    if verbosity >= 2:
        print 'Collapsed matrix dimensions: {}'.format(collapsed_matrix.shape)

    # Create new condition ids (indexed from 1, like the sample tables)
    num_conds = collapsed_matrix.shape[1]
    cond_ints = list(range(1, num_conds + 1))
    num_digits = len(str(max(cond_ints)))
    string_conds = ['collapsed-' + str(x).zfill(num_digits + 1) for x in cond_ints]

    # Add the condition ids into the sample table
    collapsed_sample_table['expt_id'] = string_conds

    # Create an array of the new collapsed condition ids
    collapsed_conditions = np.array(collapsed_sample_table[['screen_name', 'expt_id']])

    if verbosity >= 2:
        print 'Collapsed sample table:'
        print collapsed_sample_table
        print 'Collapsed conditions:'
        print collapsed_conditions[0:10]

    # Write the new sample table out to file!
    table_filename = os.path.join(output_folder, 'collapsed_sample_table.txt')
    collapsed_sample_table.to_csv(table_filename, sep = '\t', header = True, index = False)

    # And write the collapsed dataset out to file!
    data_filename = os.path.join(output_folder, 'collapsed_dataset.dump.gz')
    dump_dataset(dataset = [barcodes, collapsed_conditions, collapsed_matrix], filename = data_filename)
    
    update_version_file(output_folder, VERSION)

    return None


# If this is run from the command line (which is the intention), then run the rest of this script!
if __name__ == '__main__':

    # Get arguments

    parser = argparse.ArgumentParser(epilog = 'Note that collapsing the control samples is completely customizable by defining which groups of DMSOs should be collapsed in the sample table.')
    parser.add_argument('dataset_file', help = 'The dataset on which to perform replicate collapsing.')
    parser.add_argument('sample_table', help = 'The sample table corresponding to the dataset.')
    parser.add_argument('collapse_by_col', help = 'The column from the sample table that defines the groups of replicates to collapse (could be a "names" columns, could be something more complex).')
    parser.add_argument('cor_cutoff', help = 'The minimum correlation two replicates must have in order to be collapsed')
    parser.add_argument('output_folder', help = 'The folder to which the resulting collapsed matrix and sample table are written')
    parser.add_argument('--cols_to_keep', nargs = '+', help = 'A space-delimited list of column names from the sample table to include in the new, collapsed sample table.')
    parser.add_argument('--how_to_collapse', nargs = '+', help = 'A space-delimited list, parallel to the list from "--cols_to_keep", giving the method by which each column should be collapsed. Current options are:\n"concat" to keep all values and semicolon delimit,\n"uniq" if the values for each condition in a replicate group are the same; and\n"mean" to compute a mean and standard deviation of the values (NaNs are removed first).')
    parser.add_argument('-v', '--verbosity', type = int, default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')

    args = parser.parse_args()

    # Get the data ready to rumble!
    dataset_file = os.path.abspath(args.dataset_file)
    assert os.path.isfile(dataset_file), "Dataset file does not exist."
    dataset = load_dataset(args.dataset_file)

    if args.verbosity >= 2:
        print dataset

    assert os.path.isfile(args.sample_table), "Sample table file does not exist."
    sample_table = read_sample_table(args.sample_table)

    collapse_col = args.collapse_by_col

    # Just making sure that we handle the columns-to-keep situation correctly.
    if args.cols_to_keep is None or args.how_to_collapse is None:
        cols_to_keep = None
        how_to_collapse = None
    else:
        cols_to_keep = args.cols_to_keep
        how_to_collapse = args.how_to_collapse
        assert len(cols_to_keep) == len(how_to_collapse), '"how_to_collapse" must contain the same number of entries as "cols_to_keep"'

    # Ensure that cor_cutoff is numeric
    try:
        cor_cutoff = float(args.cor_cutoff)
    except:
        assert False, 'correlation cutoff must be a floating-point number!'

    if args.verbosity >= 2:
        print cols_to_keep
        print how_to_collapse

    output_folder = args.output_folder
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    main(dataset, sample_table, collapse_col, cor_cutoff, output_folder, cols_to_keep, how_to_collapse, args.verbosity)
