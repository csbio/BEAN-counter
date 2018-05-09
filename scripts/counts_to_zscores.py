#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

VERSION='2.3.0'

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
assert barseq_path is not None, "'BARSEQ_PATH' environment variable is not set. Please consult the instructions for setting up BEAN-counter."
sys.path.append(os.path.join(barseq_path, 'lib'))

import compressed_file_opener as cfo
import cg_file_tools as cg_file
from cg_common_functions import get_verbosity, get_barcode_table, get_sample_table, get_num_cores, parse_yaml, bool_dict
from version_printing import update_version_file
from cluster_dataset_wrappers import customize_strains, customize_conditions
from lowess import py_lowess
from contextlib import closing
from multiprocessing import Pool

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
   
    control_bool_ind = np.array([bool_dict[x] for x in sample_table['control?']])
    control_table = sample_table[control_bool_ind]
    control_screen_names = control_table['screen_name']
    control_expt_ids = control_table['expt_id']
    control_condition_ids = np.array(list(it.izip(control_screen_names, control_expt_ids)))
    
    control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, control_condition_ids)])
    final_control_condition_ids = condition_ids[control_condition_indices]

    return final_control_condition_ids

def get_control_dataset(dataset, control_condition_ids, control_detection_limit):
    
    [barcode_gene_ids, condition_ids, matrix] = dataset

    #control_condition_indices = np.array([i for i, cond_id in enumerate(condition_ids) if a_is_row_in_b(cond_id, control_condition_ids)])
    control_condition_indices = []
    # Get the ids of the conditions for which >= 75% of the profile is above
    # the control count threshold. This percentage could be changed to a
    # parameter in the future, but this is a relatively rare corner case that
    # only occurs in per-lane scoring when all controls in that lane are of low
    # counts/quality.
    for i, cond_id in enumerate(condition_ids):
        if a_is_row_in_b(cond_id, control_condition_ids):
            if np.nanmean(matrix[:, i] >= control_detection_limit) >= 0.75:
                control_condition_indices.append(i)
    control_condition_indices = np.array(control_condition_indices)

    # If there are less than two control profiles in the data, or if this many
    # are left after the above quality filtering, then return the entire
    # dataset as the control dataset.
    if control_condition_indices.size < 2:
        control_condition_ids = condition_ids
        control_matrix = matrix
    else:
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
    # The individual profiles will never be big enough to take advantage of
    # parallelization. If anything, this function should be called in parallel
    # instead. Setting num_cores to 1 just uses python's builtin `map` instead
    # of using multiprocessing with one core.
    lowess_line = py_lowess(xx[k], yy[k], f = 0.3, iter = 5, num_cores = 1)
    if np.sum(lowess_line) == 0:
        return yy_normalized
    yy_normalized[k] = xx[k] * yy[k] / lowess_line
    return yy_normalized

def normalize_one_profile(j):
    global matrix_
    global mean_control_profile_
    global sample_detection_limit_
    #print j
    y = matrix_[:, j]
    y[y < sample_detection_limit_] = sample_detection_limit_
    y_log = np.log(y)
    return smooth(mean_control_profile_, y_log)

def normalizeUsingAllControlsAndSave(config_params, outfolder, dataset, control_condition_ids, lane_id):
    global matrix_
    global mean_control_profile_
    global sample_detection_limit_

    barcode_gene_ids, condition_ids, matrix = dataset
    
    # Make sure the counts are floats for all normalization procedures
    matrix = matrix.astype(np.float).copy()

    # Get the detection limits
    sample_detection_limit, control_detection_limit = get_detection_limits(config_params)
    
    # Dump out the raw count matrix before processing
    raw_filename = os.path.join(outfolder, '{}_raw.dump.gz'.format(lane_id))
    raw_of = gzip.open(raw_filename, 'wb')
    cPickle.dump(dataset, raw_of)
    raw_of.close()

    control_matrix_gene_barcode_ids, control_matrix_condition_ids, control_matrix = get_control_dataset(dataset, control_condition_ids, control_detection_limit)

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
    matrix_ = matrix
    mean_control_profile_ = mean_control_profile
    sample_detection_limit_ = sample_detection_limit
    num_cores = get_num_cores(config_params)
    #print matrix.shape
    if num_cores > 1:
        with closing(Pool(processes = num_cores)) as pool:
            matrix = np.vstack(pool.map(normalize_one_profile, range(matrix.shape[1]))).T
    else:
        matrix = np.vstack(map(normalize_one_profile, range(matrix.shape[1]))).T
    #for j in range(matrix.shape[1]):
    #    if get_verbosity(config_params) >= 3:
    #        print j
    #    y = matrix[:, j]
    #    y[y < sample_detection_limit] = sample_detection_limit
    #    y_log = np.log(y)
    #    if get_verbosity(config_params) >= 3:
    #        print y
    #    matrix[:, j] = smooth(mean_control_profile, y_log)
    
    # Dump out the lowess-normalized matrix
    lowess_dataset = [barcode_gene_ids, condition_ids, matrix]
    lowess_filename = os.path.join(outfolder, '{}_lowess_norm.dump.gz'.format(lane_id))
    lowess_of = gzip.open(lowess_filename, 'wb')
    cPickle.dump(lowess_dataset, lowess_of)
    lowess_of.close()

    return lowess_dataset, mean_control_profile

def deviations_globalmean(config_params, outfolder, lowess_dataset, mean_control_profile, lane_id):
    barcode_gene_ids, condition_ids, matrix = lowess_dataset
    matrix = matrix.copy()
    
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
    lowess[pos] = py_lowess(np.asarray(repeated_raw_mean_control_profile[pos].ravel().tolist()[0]), dev_control_matrix_tall_squared[pos], 
                                   f=0.3, iter=1, num_cores = get_num_cores(config_params))
    lowess[neg] = py_lowess(np.asarray(repeated_raw_mean_control_profile[neg].ravel().tolist()[0]), dev_control_matrix_tall_squared[neg],
                                   f=0.3, iter=1, num_cores = get_num_cores(config_params))
    lowess_pos = np.zeros(raw_mean_control_profile.shape) + np.nan
    lowess_neg = np.zeros(raw_mean_control_profile.shape) + np.nan
    for i in range(raw_mean_control_profile.shape[0]):
        for j in range(i, dev_control_matrix_tall.shape[0], raw_mean_control_profile.shape[0]):
            if pos[j]:
                lowess_pos[i] = lowess[j]
            else:
                lowess_neg[i] = lowess[j]
   
    ## Symmetric lowess is not used
    ## Smooths symmetric matrix, removing NaNs and non-well-behaved input    
    #formatted_symmetric_x = np.asarray(repeated_raw_mean_control_profile.ravel().tolist()[0])
    #formatted_symmetric_y = np.reshape(dev_control_matrix_tall.flatten(), -1)
    #lowess_symmetric = np.zeros(formatted_symmetric_x.shape) + np.nan
    #k = wellbehaved(formatted_symmetric_x) & wellbehaved(formatted_symmetric_y)
    #lowess_result = py_lowess(formatted_symmetric_x[k], formatted_symmetric_y[k], f=0.1, iter=1)
    #for i in range(len(lowess_symmetric)):
    #    j = 0
    #    if k[i]: 
    #        lowess_symmetric[i] = lowess_result[j]
    #        j += 1
    
    #return np.sqrt( lowess_neg ).real , np.sqrt( lowess_pos ).real, np.sqrt(lowess_symmetric[range(raw_mean_control_profile.shape[0])]).real
    return np.sqrt( lowess_neg ).real , np.sqrt( lowess_pos ).real

def scaleInteractions(config_params, outfolder, deviation_dataset, raw_dataset, control_condition_ids, lane_id):
    barcode_gene_ids, condition_ids, matrix = deviation_dataset
    
    # Get the detection limits
    sample_detection_limit, control_detection_limit = get_detection_limits(config_params)
    
    scaled_dev_matrix = np.zeros(matrix.shape)
    control_matrix_gene_barcode_ids, control_matrix_condition_ids, control_matrix = get_control_dataset(deviation_dataset, control_condition_ids, control_detection_limit)

    # Get rid of controls that have an infinite or NaN as one of their strains
    did = sum(np.invert(np.add(np.isfinite(control_matrix), np.isnan(control_matrix)))) == 0
    control_matrix = control_matrix[:,did]

    # Get the control profiles from the raw dataset!
    control_raw_matrix_gene_barcode_ids, control_raw_matrix_condition_ids, control_raw_matrix = get_control_dataset(raw_dataset, control_condition_ids, control_detection_limit)
    # Set all strains in control conditions with counts below the control count detection limit to NaN
    control_raw_matrix = control_raw_matrix.astype(np.float)
    control_raw_matrix[control_raw_matrix < control_detection_limit] = np.nan
    # Compute the mean control profile on the raw control profiles (again)
    mean_control_profile = np.log( np.nanmean( control_raw_matrix, axis=1 ))

    #lowess_neg, lowess_pos, lowess_symmetric = getAsymmetricSigmaForScipyMatrix(mean_control_profile, control_matrix, config_params)
    lowess_neg, lowess_pos = getAsymmetricSigmaForScipyMatrix(mean_control_profile, control_matrix, config_params)
    for i in range(matrix.shape[1]):
        deviation = np.array( matrix[:, [i]] )
        #print 'deviation shape:', deviation.shape
        #print 'lowess_neg shape:', lowess_neg.shape
        #print 'lowess_pos shape:', lowess_pos.shape
        k = wellbehaved(deviation) # this remains same but nevertheless I have put it here
        I = np.argsort(abs(deviation[k]))
        I = I[range(int( I.shape[0]*0.75) )] # median 75%
        sigma = np.std(deviation[k][I])
        #sigma = np.zeros(mean_control_profile.shape) + sigma
        sigma = np.zeros(deviation.shape) + sigma
        total_sigma_neg = np.nanmax([sigma, lowess_neg], axis = 0)
        total_sigma_pos = np.nanmax([sigma, lowess_pos], axis = 0)
        zscores = np.zeros(deviation.shape) + np.nan
        #print 'zscores shape:', zscores.shape
        #print 'total_sigma_neg shape:', total_sigma_neg.shape
        #print 'total_sigma_pos shape:', total_sigma_pos.shape
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

def add_one_label(axes, label, x, y, final_y):
    
    if final_y < 0:
        axes.annotate(label, xy = (x, y), horizontalalignment = 'left', verticalalignment = 'top', size = 4)
    else:
        axes.annotate(label, xy = (x, y), horizontalalignment = 'right', verticalalignment = 'bottom', size = 4)

    return None

def make_individual_plots(mean_prof, prof_1, prof_2, prof_3, prof_4, strain_labels, condition_label, config_params, outfolder):

    # Set up and draw the 2x2 plot!
    fig1, fig2, fig3, fig4 = plt.figure(0), plt.figure(1), plt.figure(2), plt.figure(3)
    ax1, ax2, ax3, ax4 = [fig1.add_subplot(111),
                          fig2.add_subplot(111),
                          fig3.add_subplot(111),
                          fig4.add_subplot(111)
                          ]
    
    # Draw scatterplots
    ax1.scatter(mean_prof, prof_1, s = 1, marker = '.')
    ax2.scatter(mean_prof, prof_2, s = 1, marker = '.')
    ax3.scatter(mean_prof, prof_3, s = 1, marker = '.')
    ax4.scatter(mean_prof, prof_4, s = 1, marker = '.')


    # Add point labels
    if config_params['scatter_label_scheme'] != '0':
        for j, lab in enumerate(strain_labels):
            #print 'trying to annotate {}'.format(lab)
            if lab != '':
                print lab
                add_one_label(ax1, lab, mean_prof[j], prof_1[j], prof_4[j])
                add_one_label(ax2, lab, mean_prof[j], prof_2[j], prof_4[j])
                add_one_label(ax3, lab, mean_prof[j], prof_3[j], prof_4[j])
                add_one_label(ax4, lab, mean_prof[j], prof_4[j], prof_4[j])

    # Add common x label
    ax1.set_xlabel('Mean control profile (log counts)')
    ax2.set_xlabel('Mean control profile (log counts)')
    ax3.set_xlabel('Mean control profile (log counts)')
    ax4.set_xlabel('Mean control profile (log counts)')
    
    # Add title!
    ax1.set_title(condition_label)
    ax2.set_title(condition_label)
    ax3.set_title(condition_label)
    ax4.set_title(condition_label)

    # Add individual y labels
    ax1.set_ylabel('Read counts')
    ax2.set_ylabel('Lowess-normalized\nread counts')
    ax3.set_ylabel('Deviation from\nnormalized counts')
    ax4.set_ylabel('z-score')

    # Align y-axis labels
    #label_x = -0.2
    #ax1.yaxis.set_label_coords(label_x, 0.5)
    #ax2.yaxis.set_label_coords(label_x, 0.5)
    #ax3.yaxis.set_label_coords(label_x, 0.5)
    #ax4.yaxis.set_label_coords(label_x, 0.5)

    outfolder_1 = os.path.join(outfolder, '1_raw')
    outfolder_2 = os.path.join(outfolder, '2_lowess-normalized')
    outfolder_3 = os.path.join(outfolder, '3_deviation')
    outfolder_4 = os.path.join(outfolder, '4_scaled-deviation')
    
    for folder in [outfolder_1, outfolder_2, outfolder_3, outfolder_4]:
        if not os.path.isdir(folder):
            os.makedirs(folder)
    
    outfile_1 = os.path.join(outfolder_1, '{}_raw.pdf'.format(condition_label))
    outfile_2 = os.path.join(outfolder_2, '{}_lowess-normalized.pdf'.format(condition_label))
    outfile_3 = os.path.join(outfolder_3, '{}_deviation.pdf'.format(condition_label))
    outfile_4 = os.path.join(outfolder_4, '{}_scaled-deviation.pdf'.format(condition_label))

    # plt.savefig(outfile, bbox_inches = 'tight')
    #fig1.savefig(outfile_1, bbox_inches = 'tight')
    #fig2.savefig(outfile_2, bbox_inches = 'tight')
    #fig3.savefig(outfile_3, bbox_inches = 'tight')
    #fig4.savefig(outfile_4, bbox_inches = 'tight')
    fig1.savefig(outfile_1)
    fig2.savefig(outfile_2)
    fig3.savefig(outfile_3)
    fig4.savefig(outfile_4)

    plt.close('all')

def make_2x2_plot(mean_prof, prof_1, prof_2, prof_3, prof_4, strain_labels, condition_label, config_params, outfolder):

    # Set up and draw the 2x2 plot!
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax1, ax2, ax3, ax4 = [fig.add_subplot(221),
                          fig.add_subplot(222),
                          fig.add_subplot(223),
                          fig.add_subplot(224)
                          ]
    
    fig.subplots_adjust(left=0.2, wspace=0.6)

    # Turn off axis lines and ticks of the big subplot
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    # Draw scatterplots
    ax1.scatter(mean_prof, prof_1, s = 1, marker = '.')
    ax2.scatter(mean_prof, prof_2, s = 1, marker = '.')
    ax3.scatter(mean_prof, prof_3, s = 1, marker = '.')
    ax4.scatter(mean_prof, prof_4, s = 1, marker = '.')


    # Add point labels
    if config_params['scatter_label_scheme'] != '0':
        for j, lab in enumerate(strain_labels):
            #print 'trying to annotate {}'.format(lab)
            if lab != '':
                print lab
                add_one_label(ax1, lab, mean_prof[j], prof_1[j], prof_4[j])
                add_one_label(ax2, lab, mean_prof[j], prof_2[j], prof_4[j])
                add_one_label(ax3, lab, mean_prof[j], prof_3[j], prof_4[j])
                add_one_label(ax4, lab, mean_prof[j], prof_4[j], prof_4[j])

    # Add common x label
    ax.set_xlabel('Mean control profile (log counts)')

    # Add title!
    ax.set_title(condition_label, y = 1.08)

    # Add individual y labels
    ax1.set_ylabel('Read counts')
    ax2.set_ylabel('Lowess-normalized\nread counts')
    ax3.set_ylabel('Deviation from\nnormalized counts')
    ax4.set_ylabel('z-score')

    # Align y-axis labels
    label_x = -0.2
    ax1.yaxis.set_label_coords(label_x, 0.5)
    ax2.yaxis.set_label_coords(label_x, 0.5)
    ax3.yaxis.set_label_coords(label_x, 0.5)
    ax4.yaxis.set_label_coords(label_x, 0.5)

    combined_outfolder = os.path.join(outfolder, 'combined')
    if not os.path.isdir(combined_outfolder):
        os.makedirs(combined_outfolder)
    outfile = os.path.join(combined_outfolder, '{}_combined.pdf'.format(condition_label))

    # plt.savefig(outfile, bbox_inches = 'tight')
    plt.savefig(outfile)
    plt.close()

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
            if get_verbosity(config_params) >= 2:
                print np.sum(strains_to_label_inds)
            strain_ids_final = strain_ids_custom.copy()
            strain_ids_final[np.invert(strains_to_label_inds)] = ''
            if get_verbosity(config_params) >= 2:
                print strain_ids_final[strains_to_label_inds]
        elif config_params['scatter_label_scheme'] == '2':
            strains_to_label_inds = rankdata(-np.abs(scaled_dev_prof), method = 'min') <= float(config_params['scatter_label_cutoff'])
            if get_verbosity(config_params) >= 2:
                print np.sum(strains_to_label_inds)
            strain_ids_final = strain_ids_custom.copy()
            strain_ids_final[np.invert(strains_to_label_inds)] = ''
            if get_verbosity(config_params) >= 2:
                print strain_ids_final[strains_to_label_inds]
        #strains_to_label = scaled_dev_dataset[0][strains_to_label_inds, :]

        # Functions to plot both individual versions of the data processing
        # plots and a combined version!
        make_individual_plots(mean_control_profile, raw_prof, lowess_prof, deviation_prof, scaled_dev_prof,
                              strain_ids_final, condition_ids_final[i], config_params, scatter_outfolder)
        
        make_2x2_plot(mean_control_profile, raw_prof, lowess_prof, deviation_prof, scaled_dev_prof,
                      strain_ids_final, condition_ids_final[i], config_params, scatter_outfolder)


def get_batch_inds(sample_table, column, dataset):
    if column == 'None' or column.replace(' ', '') == '':    
        return [''], [np.arange(len(dataset[1]))]
    assert column in sample_table, '\n"sub_screen_column" parameter "{}" is not a column in the sample information table,\n' \
            'found here: {}'.format(column, config_params['sample_table_file'])
    sample_table = sample_table.set_index(['screen_name', 'expt_id'], drop = False)
    cond_tuples = [tuple(x) for x in dataset[1]]
    batches = sample_table.loc[cond_tuples, column].values
    batches_uniq = list(np.unique(batches))
    batch_index_dict = {x:i for i,x in enumerate(batches_uniq)}
    batch_index_list = [[] for i in range(len(batches_uniq))]
    for i, cond in enumerate(cond_tuples):
        batch = sample_table.loc[cond, column]
        batch_index_list[batch_index_dict[batch]].append(i)

    return batches_uniq, [np.array(x) for x in batch_index_list]

def main(config_file, lane_id):
    
    # Read in the config params
    config_params = parse_yaml(config_file)
    sample_table = get_sample_table(config_params)

    # Get the interactions output folder
    outfolder = get_lane_interactions_path(config_params, lane_id)
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    
    # Read in the count matrix from dumped file
    dataset = load_dumped_count_matrix(config_params, lane_id)

    # I think here is the best spot to split the dataset so that different
    # controls can be used for different samples.
    ### split dataset stuff!!!
    # Function returns an empty string if there is only one batch
    batch_col = str(config_params.get('sub_screen_column'))
    batches, batch_inds = get_batch_inds(sample_table, batch_col, dataset)

    batch_outfolders = [os.path.join(outfolder, 'subscreen-{}'.format(x)) if x != '' else outfolder for x in batches]
    for folder in batch_outfolders:
        if not os.path.isdir(folder):
            os.makedirs(folder)

    # Keep track of individual batch datasets so they can be joined together again.
    strains = dataset[0]
    batch_condition_list = []
    normalized_dataset_list = []
    mean_control_profile_list = []
    deviation_dataset_list = []
    scaled_dev_dataset_list = []

    for i, batch in enumerate(batches):
        batch_ind = batch_inds[i]
        # If there are no indices for this batch, then skip to the next batch
        # (unlikely but possible, and it definitely would break things)
        if len(batch_ind) == 0:
            continue
        batch_outfolder = batch_outfolders[i]
        batch_conditions = dataset[1][batch_ind]
        batch_dataset = [dataset[0], batch_conditions, dataset[2][:, batch_ind]]
    
        # Filter out samples flagged as "do not include" (include? == False)
        filtered_dataset = filter_dataset_for_include(batch_dataset, sample_table, config_params)
        batch_condition_list.append(filtered_dataset[1])

        # Get list of control samples (control? = True)
        control_condition_ids = get_control_condition_ids(batch_dataset, sample_table)

        # Proceed with algorithm to obtain chemical genetic interaction zscores (scaled deviations)
        if get_verbosity(config_params) >= 1:
            print "Normalizing ... "
        normalized_dataset, mean_control_profile = normalizeUsingAllControlsAndSave(config_params, batch_outfolder, filtered_dataset, control_condition_ids, lane_id)
        if get_verbosity(config_params) >= 1:
            print "Column means: "
            print np.nanmean(normalized_dataset[2], axis = 0)
            print "Done"
            print "Calculating deviations ... "
        normalized_dataset_list.append(normalized_dataset)
        mean_control_profile_list.append(mean_control_profile)
        deviation_dataset = deviations_globalmean(config_params, batch_outfolder, normalized_dataset, mean_control_profile, lane_id)
        if get_verbosity(config_params) >= 1:
            print "Column means: "
            print np.nanmean(deviation_dataset[2], axis = 0)
            print "Done"
            print "Scaling interactions ... "
        deviation_dataset_list.append(deviation_dataset)
        scaled_dev_dataset = scaleInteractions(config_params, batch_outfolder, deviation_dataset, filtered_dataset, control_condition_ids, lane_id)
        if get_verbosity(config_params) >= 1:
            print "Column means: "
            print np.nanmean(scaled_dev_dataset[2], axis = 0)
            print "Done"
        scaled_dev_dataset_list.append(scaled_dev_dataset)
        # For another time
        #if 'generate_scatterplots' in config_params:
        #    if config_params['generate_scatterplots'] == 'Y' and lane_id == 'all_lanes_filtered':
        #        if get_verbosity(config_params) >= 1:
        #            print "Generating scatterplots"
        #        generate_scatterplots(config_params, outfolder, mean_control_profile, filtered_dataset, normalized_dataset, deviation_dataset, scaled_dev_dataset)
    # If the number of batches was greater than 1, then the data from the
    # different stages of interaction scoring have been exported into
    # batch-specific folders. Here we reconstruct the full datasets and export
    # them.
    # If one or more batch names were specified, then still dump out - for per-lane scoring
    if (len(batches) > 1) or (batches != ['']):
        combined_conditions = np.vstack(batch_condition_list)
        combined_norm_dataset = [strains, combined_conditions, np.hstack([x[2] for x in normalized_dataset_list])]
        combined_dev_dataset = [strains, combined_conditions, np.hstack([x[2] for x in deviation_dataset_list])]
        combined_scaled_dev_dataset = [strains, combined_conditions, np.hstack([x[2] for x in scaled_dev_dataset_list])]
        with gzip.open(os.path.join(outfolder, '{}_lowess_norm.dump.gz'.format(lane_id)), 'wb') as f_norm:
            cPickle.dump(combined_norm_dataset, f_norm)
        with gzip.open(os.path.join(outfolder, '{}_deviation.dump.gz'.format(lane_id)), 'wb') as f_dev:
            cPickle.dump(combined_dev_dataset, f_dev)
        with gzip.open(os.path.join(outfolder, '{}_scaled_dev.dump.gz'.format(lane_id)), 'wb') as f_scaleddev:
            cPickle.dump(combined_scaled_dev_dataset, f_scaleddev)
        # Should I dump out mean control profile here?

    update_version_file(outfolder, VERSION)
    update_version_file(config_params['output_folder'], VERSION)

# call: python counts_to_zscores.py <config_file> <lane_id>
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: python counts_to_zscores.py <config_file> <lane_id>'
    else:
        config_file = sys.argv[1]
        lane_id = sys.argv[2]
        main(config_file, lane_id)
