#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.2.1'

# This script takes in a dataset and a corresponding sample table with
# a reduced number of rows. This script will reduce the dataset so that
# it contains only the samples corresponding to the rows in the
# dataset. It also exports a sample table that perfectly corresponds
# to the matrix (any extra rows in the sample table have been removed).
import argparse
import numpy as np, pandas as pd
import networkx as nx
import os, sys
import gzip
import cPickle

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'scripts'))
sys.path.append(os.path.join(barseq_path, 'lib'))

import cluster_dataset_wrappers as clus_wrap
from cg_common_functions import read_sample_table, read_barcode_table
from version_printing import update_version_file


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

def main(dataset, dataset_file, table, val_name, strain_table_f, strain_columns, sample_table_f, condition_columns, output_file, verbosity):

    strains, conditions, matrix = dataset

    full_fname = os.path.abspath(dataset_file)
    out_folder = full_fname.replace('.dump.gz', '')
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)

    ## Here is where the script determines what modifications to make to the row/column labels
    #if strain_table_f is not None:
    #    strain_table = read_sample_table(strain_table_f)
    #    custom_strains = clus_wrap.customize_strains(strains, strain_table, strain_columns)
    #
    #if sample_table_f is not None:
    #    sample_table = read_sample_table(sample_table_f)
    #    custom_conditions = clus_wrap.customize_conditions(conditions, sample_table, condition_columns)

    # If the table boolean is not True, write out three fiies: matrix, rownames, and colnames
    if not table:
        strain_fname = os.path.join(out_folder, 'strains.txt')
        cond_fname = os.path.join(out_folder, 'conditions.txt')
        mat_fname = os.path.join(out_folder, 'matrix.txt.gz')
        
        # Write out the strains!
        with open(strain_fname, 'wt') as strain_f:
            strain_f.write('Strain_ID\n')
            for strain in strains:
                strain_f.write(strain + '\n')

        # Write out the conditions!
        with open(cond_fname, 'wt') as cond_f:
            cond_f.write('screen_name\texpt_id\n')
            for cond in conditions:
                cond_f.write('\t'.join(cond) + '\n')

        # Write out the matrix!
        with gzip.open(mat_fname, 'wb') as mat_f:
            np.savetxt(mat_f, matrix, delimiter = '\t')
    
    # If the table boolean is True, then reshape into a table and write everything out as one table
    else:
        # Determine if custom information is to be added to the table and make sure column names check out
        if strain_table_f is not None and strain_columns is not None:
            strain_table = read_barcode_table(strain_table_f)
            strain_table_idx = {strain_table.loc[i, 'Strain_ID']:i for i in range(strain_table.shape[0])}
            strain_columns = strain_columns.split(',')
            assert all([x in strain_table.columns for x in strain_columns]), "Not all of the specified columns are in the provided strain table ({})".format(';'.join(list(set(strain_columns) - set(strain_table.columns))))
        else:
            strain_columns = []
            strain_table_idx = {}

        if sample_table_f is not None and strain_columns is not None:
            sample_table = read_sample_table(sample_table_f)
            sample_table_idx = {(sample_table.loc[i, 'screen_name'], sample_table.loc[i, 'expt_id']):i for i in range(sample_table.shape[0])}
            condition_columns = condition_columns.split(',')
            assert all([x in sample_table.columns for x in condition_columns]), "Not all of the specified columns are in the provided sample table ({})".format(';'.join(list(set(condition_columns) - set(sample_table.columns))))
        else:
            condition_columns = []
            sample_table_idx = {}

        # Determine identities of final table columns
        if val_name is None:
            val_name = 'value'
        default_cols = ['Strain_ID', 'screen_name', 'expt_id', val_name]
        extra_strain_cols = []
        extra_condition_cols = []
        for col in strain_columns:
            if col not in default_cols:
                extra_strain_cols.append(col)
        for col in condition_columns:
            if col not in default_cols:
                extra_condition_cols.append(col)

        # Create list of lists container for table-formatted data
        final_cols = default_cols + extra_strain_cols + extra_condition_cols
        final_data = [[] for x in final_cols]

        # Iterate over the rows and columns of the matrix, and fill in the lists!
        for i in range(matrix.shape[0]):
            strain_idx = strain_table_idx.get(strains[i])
            for j in range(matrix.shape[1]):
                k = 4
                condition_idx = sample_table_idx.get(tuple(conditions[j]))
                
                final_data[0].append(strains[i])
                final_data[1].append(conditions[j][0])
                final_data[2].append(conditions[j][1])
                final_data[3].append(matrix[i, j])
                for col in extra_strain_cols:
                    final_data[k].append(strain_table.loc[strain_idx, col])
                    k += 1
                for col in extra_condition_cols:
                    final_data[k].append(sample_table.loc[condition_idx, col])
                    k += 1

        # Create the data frame!
        table = pd.DataFrame(final_data).transpose()
        table.columns = final_cols

        # Write the table out to file!
        if output_file is None:
            tab_fname = os.path.join(out_folder, 'data_table.txt.gz')
        else:
            tab_fname = output_file
        with gzip.open(tab_fname, 'wb') as tab_f:
            table.to_csv(tab_f, sep = '\t', header = True, index = False)
    
    update_version_file(out_folder, VERSION)

# If this is run from the command line (which is the intention), then run the rest of this script!
if __name__ == '__main__':

    # Get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_file', help = 'The dataset to be reduced.')
    parser.add_argument('--table', action = 'store_true', help = 'Add this flag if the output should be in table form, instead of a matrix plus row and column label files. If --table is specified, the user can include additional information from the strain and sample tables to the output. For similar functionality with a matrix-formatted output, visualize_dataset.py provides output to CDT format which can be viewed in Java TreeView or a spreadsheet editor. Right now, the matrix prints out values with higher precision, but it should not make a difference for future computations.')
    parser.add_argument('--strain_table', help = 'The strain barcode table used to interpret the strain barcodes (is in the "data" folder by default). Only used if --table is specified')
    parser.add_argument('--sample_table', help = 'The sample table corresponding to the dataset. Only used if --table is specified.')
    parser.add_argument('--strain_columns', help = 'Comma-delimited list of columns from the strain table to add to the output, if --table is specified.')
    parser.add_argument('--condition_columns', help = 'Comma-delimited list of columns from the sample table to add to the output, if --table is specified.')
    parser.add_argument('--value_name', help = 'Optional argument to specify a name for the column that stores the matrix values')
    parser.add_argument('--output_file', help = 'Optional argument to specify a different output file. Default is "data_table.txt.gz", deposited inside a folder named after the input dataset.')
    parser.add_argument('-v', '--verbosity', type = int, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')

    args = parser.parse_args()

    # Handle the presence/absence of custom columns
    if args.strain_table is None or args.strain_columns is None:
        strain_table_f = None
        strain_columns = None
    else:
        strain_table_f = args.strain_table
        strain_columns = args.strain_columns

    if args.sample_table is None or args.condition_columns is None:
        sample_table_f = None
        condition_columns = None
    else:
        sample_table_f = args.sample_table
        condition_columns = args.condition_columns

    # Get the data ready to rumble!
    dataset_file = os.path.abspath(args.dataset_file)
    assert os.path.isfile(dataset_file), "Dataset file does not exist."
    dataset = load_dataset(args.dataset_file)

    if args.verbosity >= 2:
        print dataset

    table = args.table
    val_name = args.value_name

    main(dataset, dataset_file, table, val_name, strain_table_f, strain_columns, sample_table_f, condition_columns, args.output_file, args.verbosity)
