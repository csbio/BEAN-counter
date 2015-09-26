# This script will read in a count matrix and, given a sample table that
# indicates which samples are controls and which should be excluded,
# generates chemical genetic interaction z-scores.
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

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append(os.path.join(barseq_path, 'lib'))

import config_file_parser as cfp
import compressed_file_opener as cfo
import cg_file_tools as cg_file

from mlabwrap import mlab


def get_sample_table(config_params):

    filename = config_params['sample_table_file']

    # Read everything in as a string, to prevent vexing
    # number interpretation problems! Methods further down
    # can coerce to different types.
    tab = pd.read_table(filename, dtype = 'S')
    return tab

def get_lane_data_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'intermediate', lane_id)

def get_lane_interactions_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'interactions', lane_id)

def get_detection_limits(config_params):

    sample_detection_limit = float(config_params['sample_detection_limit'])
    control_detection_limit = float(config_params['control_detection_limit'])
    
    return sample_detection_limit, control_detection_limit

def get_dumped_count_matrix_filename(config_params, lane_id):

    lane_folder = get_lane_data_path(config_params, lane_id)
    return os.path.join(lane_folder, '{}_barseq_matrix.dump.gz'.format(lane_id))

def load_dumped_count_matrix(config_params, lane_id):

    filename = get_dumped_count_matrix_filename(config_params, lane_id)
    f = gzip.open(filename, 'rb')
    barcodes, conditions, matrix = cPickle.load(f)
    dataset = [np.array(barcodes), np.array(conditions), matrix]
    f.close()

    return dataset
    
def filter_dataset_for_include(dataset, sample_table):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset

    bool_dict = {'True': True, 'False': False}
    include_bool_ind = np.array([bool_dict[x] for x in sample_table['include?']])
    include_table = sample_table[include_bool_ind]
    include_screen_names = include_table['screen_name']
    include_expt_ids = include_table['expt_id']
    include_condition_ids = ['{0}-{1}'.format(*x) for x in it.izip(include_screen_names, include_expt_ids)]
    
    include_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if cond_id in include_condition_ids])
    
    filtered_condition_ids = condition_ids[include_condition_indices]
    filtered_matrix = matrix[:, include_condition_indices]

    return [barcode_gene_ids, filtered_condition_ids, filtered_matrix]

def get_control_condition_ids(dataset, sample_table):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
   
    bool_dict = {'True': True, 'False': False}
    control_bool_ind = np.array([bool_dict[x] for x in sample_table['control?']])
    control_table = sample_table[control_bool_ind]
    control_screen_names = control_table['screen_name']
    control_expt_ids = control_table['expt_id']
    control_condition_ids = ['{0}-{1}'.format(*x) for x in it.izip(control_screen_names, control_expt_ids)]
    
    control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if cond_id in control_condition_ids])
    final_control_condition_ids = condition_ids[control_condition_indices]

    return final_control_condition_ids

def get_control_dataset(dataset, control_condition_ids):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
    
    control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if cond_id in control_condition_ids])
    
    control_condition_ids = condition_ids[control_condition_indices]
    control_matrix = matrix[:, control_condition_indices]

    return [barcode_gene_ids, control_condition_ids, control_matrix]


def wellbehaved(xx):
    return ( np.invert( np.isinf(xx) ) ) & ( np.invert( np.isnan(xx) ) )

def smooth(xx, yy):
    k = wellbehaved(xx) & wellbehaved(yy)
    yy_normalized = np.zeros(xx.shape) + np.nan
    if np.sum(yy[k]) == 0:
        return yy_normalized
    lowess_line = np.array( mlab.smooth(xx[k], yy[k], 0.3, 'rlowess')  ).transpose()[0]
    if np.sum(lowess_line) == 0:
        return yy_normalized
    yy_normalized[k] = xx[k] * yy[k] / lowess_line
    return yy_normalized

def normalizeUsingAllControlsAndSave(config_params, outfolder, dataset, control_condition_ids, lane_id):
    barcode_gene_ids, condition_ids, matrix = dataset
    
    # Make sure the counts are floats for all normalization procedures
    matrix = matrix.astype(np.float)

    # Get the detection limits
    sample_detection_limit, control_detection_limit = get_detection_limits(config_params)
    
    # Dump out the raw count matrix before processing
    raw_filename = os.path.join(outfolder, '{}_raw.dump.gz'.format(lane_id))
    raw_of = gzip.open(raw_filename, 'wb')
    cPickle.dump(dataset, raw_of)
    raw_of.close()

    control_matrix_gene_barcode_ids, control_matrix_condition_ids, control_matrix = get_control_dataset(dataset, control_condition_ids)

    # Convert the control matrix to floats so it can incorporate NaNs
    control_matrix = control_matrix.astype(np.float)

    # Set all strains in control conditions with counts below the control count detection limit to NaN
    control_matrix[control_matrix < control_detection_limit] = np.nan
    # Compute the mean control profile
    mean_control_profile = np.log( np.nanmean( control_matrix, axis=1 ))
    # print mean_control_profile
    # print mean_control_profile.mean()
    #print x.shape

    # Replace each profile in the matrix with a smoothed profile
    # In the process, set the counts for any strain under the sample count detection limit
    # to the sample count detection limit. This prevents logging of zero values
    for j in range(matrix.shape[1]):
        #print j
        y = matrix[:, j]
        y[y < sample_detection_limit] = sample_detection_limit
        y_log = np.log(y)
        #print y
        matrix[:, j] = smooth(mean_control_profile, y_log)
    
    # Dump out the lowess-normalized matrix
    lowess_dataset = [barcode_gene_ids, condition_ids, matrix]
    lowess_filename = os.path.join(outfolder, '{}_lowess_norm.dump.gz'.format(lane_id))
    lowess_of = gzip.open(lowess_filename, 'wb')
    cPickle.dump(lowess_dataset, lowess_of)
    lowess_of.close()

    return lowess_dataset, mean_control_profile

def deviations_globalmean(config_params, outfolder, lowess_dataset, mean_control_profile, lane_id):
    barcode_gene_ids, condition_ids, matrix = lowess_dataset
    
    # control_matrix_gene_barcode_ids, control_matrix_condition_ids, control_matrix = get_control_dataset(lowess_dataset, control_condition_ids)
    
    #print type(allcontrols)
    # Subtract the mean control profile from each profile in the lowess-normalized matrix
    # mean_control_profile = np.log( np.nanmean( control_matrix, axis=1 ))
    for j in range(matrix.shape[1]):
        matrix[:, j] -= mean_control_profile

    # Dump out the deviation matrix
    deviation_dataset = [barcode_gene_ids, condition_ids, matrix]
    deviation_filename = os.path.join(outfolder, '{}_deviation.dump.gz'.format(lane_id))
    deviation_of = gzip.open(deviation_filename, 'wb')
    cPickle.dump(deviation_dataset, deviation_of)
    deviation_of.close()

    return deviation_dataset

def getAsymmetricSigmaForScipyMatrix(raw_mean_control_profile, dev_control_matrix):
    dev_control_matrix_tall = np.zeros((0, 1))
    repeated_raw_mean_control_profile = np.zeros((0, 1))
    raw_mean_control_profile = np.matrix(raw_mean_control_profile).transpose()
    for i in range(dev_control_matrix.shape[1]):
        #print allcontrols.shape, diffcontrols[:, [i]].shape
        dev_control_matrix_tall = np.vstack((dev_control_matrix_tall, dev_control_matrix[:, [i]]))
        repeated_raw_mean_control_profile = np.vstack((repeated_raw_mean_control_profile, raw_mean_control_profile))
    dev_control_matrix_tall = np.array(dev_control_matrix_tall)
    pos = dev_control_matrix_tall >= 0
    neg = dev_control_matrix_tall < 0
    dev_control_matrix_tall_squared = dev_control_matrix_tall**2
    lowess = np.zeros(dev_control_matrix_tall.shape)
    #print xs[pos],  allcontrols_squared[pos]
    #print xs.shape, allcontrols_squared.shape, pos.shape, scipy.nonzero(pos)[0].shape[0]
    lowess[pos] = np.array( mlab.smooth(repeated_raw_mean_control_profile[pos],  dev_control_matrix_tall_squared[pos] , 0.3, 'lowess')  ).transpose()[0]
    lowess[neg] = np.array( mlab.smooth(repeated_raw_mean_control_profile[neg],  dev_control_matrix_tall_squared[neg] , 0.3, 'lowess')  ).transpose()[0]
    lowess_pos = np.zeros(raw_mean_control_profile.shape) + np.nan
    lowess_neg = np.zeros(raw_mean_control_profile.shape) + np.nan
    for i in range(raw_mean_control_profile.shape[0]):
        for j in range(i, dev_control_matrix_tall.shape[0], raw_mean_control_profile.shape[0]):
            if pos[j]:
                lowess_pos[i] = lowess[j]
            else:
                lowess_neg[i] = lowess[j]
    #print scipy.nonzero(wellbehaved(xs))[0].shape[0]
    #print scipy.nonzero(wellbehaved(allcontrols))[0].shape[0]
    lowess_symmetric = np.array( mlab.smooth(repeated_raw_mean_control_profile,  dev_control_matrix_tall, 'lowess')  ).transpose()[0]

    return np.sqrt( lowess_neg ).real , np.sqrt( lowess_pos ).real, np.sqrt(lowess_symmetric[range(raw_mean_control_profile.shape[0])]).real

def scaleInteractions(config_params, outfolder, deviation_dataset, raw_dataset, control_condition_ids, lane_id):
    barcode_gene_ids, condition_ids, matrix = deviation_dataset
    
    # Get the detection limits
    sample_detection_limit, control_detection_limit = get_detection_limits(config_params)
    
    scaled_dev_matrix = np.zeros(matrix.shape)
    control_matrix_gene_barcode_ids, control_matrix_condition_ids, control_matrix = get_control_dataset(deviation_dataset, control_condition_ids)

    # Get rid of controls that have an infinite or NaN as one of their strains
    did = sum(np.invert(np.add(np.isfinite(control_matrix), np.isnan(control_matrix)))) == 0
    control_matrix = control_matrix[:,did]

    # Get the control profiles from the raw dataset!
    control_raw_matrix_gene_barcode_ids, control_raw_matrix_condition_ids, control_raw_matrix = get_control_dataset(raw_dataset, control_condition_ids)
    # Set all strains in control conditions with counts below the control count detection limit to NaN
    control_raw_matrix = control_raw_matrix.astype(np.float)
    control_raw_matrix[control_raw_matrix < control_detection_limit] = np.nan
    # Compute the mean control profile on the raw control profiles (again)
    mean_control_profile = np.log( np.nanmean( control_raw_matrix, axis=1 ))

    lowess_neg, lowess_pos, lowess_symmetric = getAsymmetricSigmaForScipyMatrix(mean_control_profile, control_matrix)
    for i in range(control_matrix.shape[1]):
        deviation = np.array( control_matrix[:, [i]] )
        k = wellbehaved(deviation) # this remains same but nevertheless I have put it here
        I = np.argsort(abs(deviation[k]))
        I = I[range(int( I.shape[0]*0.75) )] # median 75%
        sigma = np.std(deviation[k][I])
        sigma = np.zeros(mean_control_profile.shape) + sigma
        total_sigma_neg = np.maximum(sigma, lowess_neg)
        total_sigma_pos = np.maximum(sigma, lowess_pos)
        zscores = np.zeros(deviation.shape) + np.nan
        zscores[deviation < 0] = deviation[deviation < 0] / total_sigma_neg[deviation < 0]
        zscores[deviation > 0] = deviation[deviation > 0] / total_sigma_pos[deviation > 0]
        zscores[deviation == 0] = 0
        scaled_dev_matrix[:, i] = zscores.transpose()

    # Dump out the iscaled deviation matrix
    scaled_dev_dataset = [barcode_gene_ids, condition_ids, scaled_dev_matrix]
    scaled_dev_filename = os.path.join(outfolder, '{}_scaled_dev.dump.gz'.format(lane_id))
    scaled_dev_of = gzip.open(scaled_dev_filename, 'wb')
    cPickle.dump(scaled_dev_dataset, scaled_dev_of)
    scaled_dev_of.close()

    return scaled_dev_dataset


def main(config_file, lane_id):
    
    # Read in the config params
    print 'parsing parameters...'
    config_params = cfp.parse(config_file)
    sample_table = get_sample_table(config_params)

    # Get the interactions output folder
    outfolder = get_lane_interactions_path(config_params, lane_id)
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    
    # Read in the count matrix from dumped file
    dataset = load_dumped_count_matrix(config_params, lane_id)

    # Filter out samples flagged as "do not include" (include? == True)
    filtered_dataset = filter_dataset_for_include(dataset, sample_table)

    # Get list of control samples (control? = True)
    control_condition_ids = get_control_condition_ids(dataset, sample_table)

    # Proceed with algorithm to obtain chemical genetic interaction zscores (scaled deviations)
    print "Normalizing ... ",
    normalized_dataset, mean_control_profile = normalizeUsingAllControlsAndSave(config_params, outfolder, filtered_dataset, control_condition_ids, lane_id)
    print "Done"
    print "Calculating deviations ... ",
    deviation_dataset = deviations_globalmean(config_params, outfolder, normalized_dataset, mean_control_profile, lane_id)
    print "Done"
    print "Scaling interactions ... ",
    scaled_dev_dataset = scaleInteractions(config_params, outfolder, deviation_dataset, filtered_dataset, control_condition_ids, lane_id)
    print "Done"
    
# call: python counts_to_zscores.py <config_file> <lane_id>
if '__name__' == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: python counts_to_zscores.py <config_file> <lane_id>'
    else:
        config_file = sys.argv[1]
        lane_id = sys.argv[2]
        main(config_file, lane_id)
