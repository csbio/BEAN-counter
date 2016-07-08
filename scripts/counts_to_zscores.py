#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script will read in a count matrix and, given a sample table that
# indicates which samples are controls and which should be excluded,
# generates chemical genetic interaction z-scores.
import pandas as pd
import numpy as np
import scipy
from scipy.stats import rankdata
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
from cg_common_functions import *
from cluster_dataset_wrappers import customize_strains, customize_conditions

sys.path.append(os.path.join(barseq_path, 'lib/python2.7/site-packages')) 
from mlabwrap import mlab


def get_sample_table(config_params):

    filename = config_params['sample_table_file']

    # Read everything in as a string, to prevent vexing
    # number interpretation problems! Methods further down
    # can coerce to different types.
    tab = pd.read_table(filename, dtype = 'S')
    return tab

def get_barcode_table(config_params):

    species_config_params = get_species_config_params(config_params)

    barseq_path = os.getenv('BARSEQ_PATH')
    filename = species_config_params['gene_barcode_file']
    full_path = os.path.join(barseq_path, 'data/barcodes', filename)

    tab = pd.read_table(full_path, dtype = 'S')
    return tab

def get_species_config_params(config_params):
    barseq_path = os.getenv('BARSEQ_PATH')
    species_config_file = os.path.join(barseq_path, 'data/species_config_file.txt')
    all_species_params = cfp.parse_species_config(species_config_file)
    species_id = config_params['species_ID']
    species_params = all_species_params[species_id]
    return species_params

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

def a_is_row_in_b(a, b):

    return np.any(np.all(a == b, axis = 1))
    
def filter_dataset_for_include(dataset, sample_table, config_params):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset

    bool_dict = {'True': True, 'TRUE': True, 'False': False, 'FALSE': False}
    include_bool_ind = np.array([bool_dict[x] for x in sample_table['include?']])
    include_table = sample_table[include_bool_ind]
    include_screen_names = include_table['screen_name']
    include_expt_ids = include_table['expt_id']
    include_condition_ids = np.array(list(it.izip(include_screen_names, include_expt_ids)))

    include_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, include_condition_ids)])
    
    if get_verbosity(config_params) >= 2:
        print condition_ids
        print include_condition_ids
        print include_condition_indices
    filtered_condition_ids = condition_ids[include_condition_indices]
    filtered_matrix = matrix[:, include_condition_indices]

    return [barcode_gene_ids, filtered_condition_ids, filtered_matrix]

def get_control_condition_ids(dataset, sample_table):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
   
    bool_dict = {'True': True, 'TRUE': True, 'False': False, 'FALSE': False}
    control_bool_ind = np.array([bool_dict[x] for x in sample_table['control?']])
    control_table = sample_table[control_bool_ind]
    control_screen_names = control_table['screen_name']
    control_expt_ids = control_table['expt_id']
    control_condition_ids = np.array(list(it.izip(control_screen_names, control_expt_ids)))
    
    control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, control_condition_ids)])
    final_control_condition_ids = condition_ids[control_condition_indices]

    return final_control_condition_ids

def get_control_dataset(dataset, control_condition_ids):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset
    
    control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, control_condition_ids)])
    
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
    if get_verbosity(config_params) >= 2:
        print 'mean control profile:'
        print mean_control_profile
        print 'mean control profile shape:'
        print mean_control_profile.shape
        print 'mean of mean control profile:'
        print np.nanmean(mean_control_profile)

    # Replace each profile in the matrix with a smoothed profile
    # In the process, set the counts for any strain under the sample count detection limit
    # to the sample count detection limit. This prevents logging of zero values
    for j in range(matrix.shape[1]):
        if get_verbosity(config_params) >= 3:
            print j
        y = matrix[:, j]
        y[y < sample_detection_limit] = sample_detection_limit
        y_log = np.log(y)
        if get_verbosity(config_params) >= 3:
            print y
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

def getAsymmetricSigmaForScipyMatrix(raw_mean_control_profile, dev_control_matrix, config_params):
    dev_control_matrix_tall = np.zeros((0, 1))
    repeated_raw_mean_control_profile = np.zeros((0, 1))
    raw_mean_control_profile = np.matrix(raw_mean_control_profile).transpose()
    for i in range(dev_control_matrix.shape[1]):
        dev_control_matrix_tall = np.vstack((dev_control_matrix_tall, dev_control_matrix[:, [i]]))
        repeated_raw_mean_control_profile = np.vstack((repeated_raw_mean_control_profile, raw_mean_control_profile))
    dev_control_matrix_tall = np.array(dev_control_matrix_tall)
    pos = dev_control_matrix_tall >= 0
    neg = dev_control_matrix_tall < 0
    dev_control_matrix_tall_squared = dev_control_matrix_tall**2
    lowess = np.zeros(dev_control_matrix_tall.shape)
    for i in range(dev_control_matrix.shape[1]):
        if get_verbosity(config_params) >= 3:
            print dev_control_matrix_tall[pos],  dev_control_matrix_tall_squared[pos]
            print dev_control_matrix_tall.shape, dev_control_matrix_tall_squared.shape, pos.shape, scipy.nonzero(pos)[0].shape[0]
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

    lowess_neg, lowess_pos, lowess_symmetric = getAsymmetricSigmaForScipyMatrix(mean_control_profile, control_matrix, config_params)
    for i in range(matrix.shape[1]):
        deviation = np.array( matrix[:, [i]] )
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

def generate_scatterplots(config_params, outfolder, mean_control_profile, raw_dataset, lowess_dataset, deviation_dataset, scaled_dev_dataset):

    scatter_outfolder = os.path.join(outfolder, 'scatterplots')
    if not os.path.isdir(scatter_outfolder):
        os.makedirs(scatter_outfolder)
        
    if config_params['scatter_label_scheme'] != '0':
        assert config_params['scatter_label_scheme'].isdigit(), '"scatter_label_scheme" must be an integer!'
        try:
            float(config_params['scatter_label_cutoff'])
        except ValueError:
            assert False, "'scatter_label_cutoff' must be a float!"

    # Load in sample table
    sample_table = get_sample_table(config_params)

    # Load in barcode table
    barcode_table = get_barcode_table(config_params)

    # Get formatted condition names
    condition_ids = raw_dataset[1]
    if config_params['scatter_condition_columns'] == '':
        condition_fmt_string = 'screen_name,expt_id'
    else:
        condition_fmt_string = config_params['scatter_condition_columns']
    condition_ids_final = customize_conditions(condition_ids, sample_table, condition_fmt_string)

    # If specified, get formatted strain names
    # If the strain labeling scheme is 0, then all the labels are blank!
    
    if config_params['scatter_label_scheme'] == '0':
        pass
        # strain_ids_final = ['' for i in range(raw_dataset[0].shape[0])]
    elif config_params['scatter_label_scheme'] in ['1', '2']:
        strain_ids = raw_dataset[0]
        strain_fmt_string = config_params['scatter_strain_columns']
        assert strain_fmt_string != '', 'scatter_strain_columns must be specified if scatter_label_scheme is not "0"'
        strain_ids_custom = customize_strains(strain_ids, barcode_table, strain_fmt_string)
    else:
        assert False, '"scatter_label_scheme" must be one of [0, 1, 2]!'
    

    # Loop over the columns of the datasets and make plots!
    for i in range(raw_dataset[2].shape[1]):

        # Extract the profiles to work with
        raw_prof = np.log(raw_dataset[2][:, i].squeeze())
        lowess_prof = lowess_dataset[2][:, i].squeeze()
        deviation_prof = deviation_dataset[2][:, i].squeeze()
        scaled_dev_prof = scaled_dev_dataset[2][:, i].squeeze()

        if get_verbosity(config_params) >= 2:
            print raw_prof
            print lowess_prof
            print deviation_prof
            print scaled_dev_prof

        # If labeling is to occur, get which strains to label from the final profile
        if config_params['scatter_label_scheme'] == '1':
            strains_to_label_inds = np.abs(scaled_dev_prof) > float(config_params['scatter_label_cutoff'])
            strain_ids_final = strain_ids_custom.copy()
            strain_ids_final[np.invert(strains_to_label_inds)] = ''
        elif config_params['scatter_label_scheme'] == '2':
            strains_to_label_inds = rankdata(-np.abs(scaled_dev_prof), method = 'min') <= float(config_params['scatter_label_cutoff'])
            strain_ids_final = strain_ids_custom.copy()
            strain_ids_final[np.invert(strains_to_label_inds)] = ''
        #strains_to_label = scaled_dev_dataset[0][strains_to_label_inds, :]

        # Set up and draw the 2x2 plot!
        f, axarr = plt.subplots(2, 2, sharex = True)
        
        # also create a big plot so x axis label is common
        #ax = f.add_subplot(111)

        # Turn off axis lines and ticks of the big subplot
        #ax.spines['top'].set_color('none')
        #ax.spines['bottom'].set_color('none')
        #ax.spines['left'].set_color('none')
        #ax.spines['right'].set_color('none')
        #ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

        # Draw scatterplots
        axarr[0, 0].scatter(mean_control_profile, raw_prof)
        axarr[0, 1].scatter(mean_control_profile, lowess_prof)
        axarr[1, 0].scatter(mean_control_profile, deviation_prof)
        axarr[1, 1].scatter(mean_control_profile, scaled_dev_prof)

        # Add point labels
        #if config_params['scatter_label_scheme'] == '0':
        #    for j, lab in enumerate(strain_ids_final):
        #        axarr[0, 0].annotate(lab, xy = (mean_control_profile[j], raw_prof[j]))
        #    for j, lab in enumerate(strain_ids_final):
        #        axarr[0, 1].annotate(lab, xy = (mean_control_profile[j], lowess_prof[j]))
        #    for j, lab in enumerate(strain_ids_final):
        #        axarr[1, 0].annotate(lab, xy = (mean_control_profile[j], deviation_prof[j]))
        #    for j, lab in enumerate(strain_ids_final):
        #        axarr[1, 1].annotate(lab, xy = (mean_control_profile[j], scaled_dev_prof[j]))

        # Add common x label
        #ax.set_xlabel('mean control profile (log counts)')

        # Add title!
        #ax.set_title(condition_ids_final[i])

        # Add individual y labels
        axarr[0, 0].set_ylabel('Read counts')
        axarr[0, 1].set_ylabel('Lowess-normalized read counts')
        axarr[1, 0].set_ylabel('Deviation from normalized counts')
        axarr[1, 1].set_ylabel('z-score')

        outfile = os.path.join(scatter_outfolder, '{}.pdf'.format(condition_ids_final[i]))

        plt.savefig(outfile, bbox_inches = 'tight')
        plt.close()



def main(config_file, lane_id):
    
    # Read in the config params
    config_params = cfp.parse(config_file)
    sample_table = get_sample_table(config_params)

    # Get the interactions output folder
    outfolder = get_lane_interactions_path(config_params, lane_id)
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    
    # Read in the count matrix from dumped file
    dataset = load_dumped_count_matrix(config_params, lane_id)

    # Filter out samples flagged as "do not include" (include? == True)
    filtered_dataset = filter_dataset_for_include(dataset, sample_table, config_params)

    # I think here is the best spot to split the dataset so that different
    # controls can be used for different samples.
    ### split dataset stuff!!!

    # Get list of control samples (control? = True)
    control_condition_ids = get_control_condition_ids(dataset, sample_table)

    # Proceed with algorithm to obtain chemical genetic interaction zscores (scaled deviations)
    if get_verbosity(config_params) >= 1:
        print "Normalizing ... "
    normalized_dataset, mean_control_profile = normalizeUsingAllControlsAndSave(config_params, outfolder, filtered_dataset, control_condition_ids, lane_id)
    if get_verbosity(config_params) >= 1:
        print "Column means: "
        print np.nanmean(normalized_dataset[2], axis = 0)
        print "Done"
        print "Calculating deviations ... "
    deviation_dataset = deviations_globalmean(config_params, outfolder, normalized_dataset, mean_control_profile, lane_id)
    if get_verbosity(config_params) >= 1:
        print "Column means: "
        print np.nanmean(deviation_dataset[2], axis = 0)
        print "Done"
        print "Scaling interactions ... "
    scaled_dev_dataset = scaleInteractions(config_params, outfolder, deviation_dataset, filtered_dataset, control_condition_ids, lane_id)
    if get_verbosity(config_params) >= 1:
        print "Column means: "
        print np.nanmean(scaled_dev_dataset[2], axis = 0)
        print "Done"
    if 'generate_scatterplots' in config_params:
	if config_params['generate_scatterplots'] == 'Y' and lane_id == 'all_lanes_filtered':
            if get_verbosity(config_params) >= 1:
                print "Generating scatterplots"
            generate_scatterplots(config_params, outfolder, mean_control_profile, filtered_dataset, normalized_dataset, deviation_dataset, scaled_dev_dataset)

# call: python counts_to_zscores.py <config_file> <lane_id>
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: python counts_to_zscores.py <config_file> <lane_id>'
    else:
        config_file = sys.argv[1]
        lane_id = sys.argv[2]
        main(config_file, lane_id)
