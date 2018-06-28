#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.5.0'

# This script takes in a dataset in the original BEAN-counter
# format, with strains identified by pairs of (Barcode, Strain_ID),
# and converts it to the new format, with only unique Strain_IDs
# keeping track of the strains. This requires the user to specify
# old and new barcode tables, which will be matched based on the
# barcode (as nothing else may be consistent!).

import argparse

if __name__ == '__main__':

    # Get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset_file', help = 'The dataset to be reduced.')
    #parser.add_argument('original_barcode_table', help = 'The barcode table corresponding to the dataset in its current format')
    parser.add_argument('new_barcode_table', help = 'The barcode table corresponding to the dataset in its destination format')
    parser.add_argument('output_file', help = 'The file to which the strain-converted dataset will be written.')
    parser.add_argument('-v', '--verbosity', default = 1, help = 'The level of verbosity printed to stdout. Ranges from 0 to 3, 1 is default.')

    args = parser.parse_args()


import numpy as np, pandas as pd
import os, sys
import gzip
import cPickle

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'lib'))

from cg_common_functions import read_barcode_table
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

def main(dataset, new_barcode_tab, output_file):

    barcodes, conditions, matrix = dataset

    orig_barcode_labels = [x[1] for x in barcodes]
    new_barcode_to_strain_id = {row[1]['Barcode']:row[1]['Strain_ID'] for row in new_barcode_tab.iterrows()}

    new_barcodes = [new_barcode_to_strain_id[x] for x in orig_barcode_labels]

    of = gzip.open(output_file, 'wb')
    cPickle.dump([new_barcodes, conditions, matrix], of)
    of.close()

if __name__ == '__main__':

    # Get the data ready to rumble!
    dataset_file = os.path.abspath(args.dataset_file)
    assert os.path.isfile(dataset_file), "Dataset file does not exist."
    dataset = load_dataset(args.dataset_file)
    if args.verbosity >= 2:
        print dataset

    #assert os.path.isfile(args.orig_barcode_table), "Original barcode table file does not exist."
    #orig_barcode_table = read_barcode_table(args.orig_barcode_table)
    
    assert os.path.isfile(args.new_barcode_table), "New barcode table file does not exist."
    new_barcode_table = read_barcode_table(args.new_barcode_table)

    output_folder = os.path.dirname(args.output_file)
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    main(dataset, new_barcode_table, args.output_file)

