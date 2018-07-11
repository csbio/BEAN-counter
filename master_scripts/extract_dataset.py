#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.6.0'

# This brief script takes in a 3D dataset (stacked matrices) and exports,
# as a similar format, one of the stacked matrices, with its row and
# column labels.
import numpy as np
import gzip
import cPickle
import argparse
import os, sys

barseq_path = os.getenv('BARSEQ_PATH')
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'lib'))

from version_printing import update_version_file

def main(dataset, n, folder):

    assert len(dataset) == 4, "Dataset does not contain one or more of the following:\nNumber of components removed vector, strain_barcode_vector,\ncondition vector, or 3-D matrix"

    if not os.path.isdir(folder):
        os.makedirs(folder)

    extracted_dataset = np.array([dataset[1], dataset[2], dataset[3][n]])
    
    filename = os.path.join(folder, '{}_components_removed.dump.gz'.format(n))
    dump_dataset(extracted_dataset, filename)

    update_version_file(folder, VERSION)
   
def dump_dataset(dataset, filename):

    f = gzip.open(filename, 'wb')
    cPickle.dump(dataset, f)
    f.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stacked_dataset_file', help = 'The dataset of stacked matrices, from which one matrix will be extracted.')
    parser.add_argument('n', type = int, help = 'The index of the matrix to extract (Python indexes from zero)')
    args = parser.parse_args()

    filename = args.stacked_dataset_file
    f = gzip.open(filename)
    dataset = cPickle.load(f)
    assert filename.endswith('.dump.gz'), 'Not a valid dataset file - must have a ".dump.gz" extension'
    folder = filename.replace('.dump.gz', '')

    main(dataset, args.n, folder)
