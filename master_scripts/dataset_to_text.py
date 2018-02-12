#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.2.0'

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
from cg_common_functions import read_sample_table
from version_printing import update_version_file

#def read_sample_table(tab_filename):
#
#    assert os.path.isfile(tab_filename), 'File {} does not exist'.format(tab_filename)
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

def main(dataset, dataset_file, table, val_name, strain_table_f, strain_columns, sample_table_f, condition_columns, verbosity):

    strains, conditions, matrix = dataset

    full_fname = os.path.abspath(dataset_file)
    out_folder = full_fname.replace('.dump.gz', '')
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)

    # Here is where the script determines what modifications to make to the row/columns labels
    if strain_table_f is not None:
        strain_table = read_sample_table(strain_table_f)
        custom_strains = clus_wrap.customize_strains(strains, strain_table, strain_columns)

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
        # Create containers for the rowname, colname, and value columns of the new table
        Strain_IDs = []
        screen_names = []
        expt_ids = []
        values = []
        
        # Iterate over the rows and columns of the matrix, and fill in the lists!
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                Strain_IDs.append(strains[i])
                screen_names.append(conditions[j][0])
                expt_ids.append(conditions[j][1])
                values.append(matrix[i, j])

        # Create the data frame!
        if val_name is None:
            val_name = 'value'
        table = pd.DataFrame({'Strain_ID': Strain_IDs, 'screen_name': screen_names, 'expt_id': expt_ids, val_name: values})

        # Reorder the data frame columns!
        table = table.reindex(columns = ['Strain_ID', 'screen_name', 'expt_id', val_name])

        # Write the table out to file!
        tab_fname = os.path.join(out_folder, 'data_table.txt.gz')
        with gzip.open(tab_fname, 'wb') as tab_f:
            table.to_csv(tab_f, sep = '\t', header = True, index = False)
    
    update_version_file(out_folder, VERSION)

# If this is run from the command line (which is the intention), then run the rest of this script!
if __name__ == '__main__':

    # Get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_file', help = 'The dataset to be reduced.')
    parser.add_argument('--table', action = 'store_true', help = 'Add this flag if the output should be in table form, instead of a matrix plus row and column label files. Right now, the matrix prints out values with higher precision, but it should not make a difference for future computations.')
    parser.add_argument('--strain_table', help = 'The strain barcode table used to interpret the strain barcodes (is in the "data" folder by default).')
    parser.add_argument('--sample_table', help = 'The sample table corresponding to the dataset.')
    parser.add_argument('--strain_columns', help = 'The columns from the barcode table to be included in the visualization.')
    parser.add_argument('--condition_columns', help = 'The columns from the sample table to be included in the visualization.')
    parser.add_argument('--value_name', help = 'Optional argument to specify a name for the column that stores the matrix values')
    parser.add_argument('-v', '--verbosity', type = int, default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')

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

    main(dataset, dataset_file, table, val_name, strain_table_f, strain_columns, sample_table_f, condition_columns, args.verbosity)
