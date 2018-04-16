#!/usr/bin/env python

#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# This script takes in the barseq_counter main configuration and species configuration files
# and a sequencing lane identifier. It exports 1) reports on the total number
# of reads, the number of reads matching the common primer, and the number of reads that
# match index tags and barcodes; 2) a **temporary** file of "index_tag\tbarcode" for all
# sequences matching the common primer; and 3) a matrix of counts, rows being the barcodes
# (genes) and columns being the conditions.

import pandas as pd
import numpy as np
import sys, os, gzip
import jellyfish as jf
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import cPickle
import itertools as it
from Bio.trie import trie

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append(os.path.join(barseq_path, 'lib'))

import compressed_file_opener as cfo
import cg_file_tools as cg_file
from cg_common_functions import get_verbosity, get_sample_table, get_amplicon_struct_params, get_barcode_table, parse_yaml

import pdb

def get_lane_location_table(config_params):

    filename = config_params['lane_location_file']
    tab = pd.read_table(filename, dtype = 'S')
    return tab

def get_lane_folder(lane_id, lane_location_tab):

    lane_location_tab = lane_location_tab.set_index('lane')
    lane_folder = lane_location_tab.loc[lane_id, 'location']
    return lane_folder

def get_lane_data_paths(config_params, lane_id):

    if 'tmp_output_dir' in config_params:
        if os.path.isdir(config_params['tmp_output_dir']):
            return [os.path.join(config_params['tmp_output_dir'], 'intermediate', lane_id),
                    os.path.join(config_params['output_folder'], 'intermediate', lane_id)
                    ]
    
    return [os.path.join(config_params['output_folder'], 'intermediate', lane_id),
            os.path.join(config_params['output_folder'], 'intermediate', lane_id)
            ]

def get_lane_reports_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'reports', lane_id)

def get_barseq_filename(config_params, lane_id):

    path = get_lane_data_paths(config_params, lane_id)[0]
    return os.path.join(path, '{0}_{1}'.format(lane_id, 'barseq.txt'))

def remove_barseq_file(config_params, lane_id):
    fname = get_barseq_filename(config_params, lane_id)
    os.remove(fname)
    return None

#def get_barcode_to_gene(config_params):
#    # Need to update for paired-end configuration
#    #barseq_path = os.getenv('BARSEQ_PATH')
#
#    amplicon_struct_params = get_amplicon_struct_params(config_params)
#    
#    filename = config_params['gene_barcode_file']
#    #full_path = os.path.join(barseq_path, 'data', 'barcodes', filename)
#    barcode_tab = read_barcode_table(filename)
#
#    # Get all possible common primer/index tag/barcode parameters, then determine
#    # how to proceed.
#    read_1_params = get_seq_params(amplicon_struct_params, 'read_1')
#    read_2_params = get_seq_params(amplicon_struct_params, 'read_2')
#    read_type_dict = determine_read_type(read_1_params, read_2_params)
#
#    # Which barcode columns are required?
#    barcode_cols = [amplicon_struct_params[x]['genetic_barcode']['barcode_file_column'] for x in read_type_dict['barcode']]
#    assert all([x in barcode_tab.columns for x in barcode_cols]), 'Not all specified barcode columns are in the barcode table.\nRequired columns: {}\nBarcode table: {}'.format(barcode_cols, filename)
#
#    # Check for duplicates in the barcodes
#    num_dup_barcodes = sum(barcode_tab[barcode_cols].duplicated())
#    assert num_dup_barcodes == 0, '{} duplicate barcodes detected. Please review the barcode to strain mapping, found here: {}'.format(num_dup_barcodes, filename)
#
#    # If there is only one read containing the barcode, add a dummy
#    # column of empty strings for the other read.
#    if len(read_type_dict['barcode']) == 1:
#        barcode_tab['dummy_barcode'] = ['' for i in range(barcode_tab.shape[0])]
#        if read_type_dict['barcode'] == ['read_1']:
#            barcode_cols.append('dummy_barcode')
#        else:
#            barcode_cols.insert(0, 'dummy_barcode')
#
#    # Convert the table to a dictionary from barcode(s) to strain_id
#    # Note: the key is always a length-2 tuple. If there is only one
#    # barcode, the entry for the other barcode is just empty strings.
#    barcode_to_strain = barcode_tab.set_index(barcode_cols)['Strain_ID'].to_dict()
#
#    return barcode_to_strain
#
#def get_index_tag_to_condition(config_params, lane_id):
#
#    sample_table = get_sample_table(config_params)
#    sample_table = sample_table.set_index('lane')
#    sample_table_lane = sample_table.loc[lane_id]
#
#    # Get all possible common primer/index tag/barcode parameters, then determine
#    # how to proceed.
#    amplicon_struct_params = get_amplicon_struct_params(config_params)
#    read_1_params = get_seq_params(amplicon_struct_params, 'read_1')
#    read_2_params = get_seq_params(amplicon_struct_params, 'read_2')
#    read_type_dict = determine_read_type(read_1_params, read_2_params)
#
#    # Which index tag columns are required?
#    index_tag_cols = [amplicon_struct_params[x]['index_tag']['sample_table_column'] for x in ['read_1', 'read_2'] if x in amplicon_struct_params]
#    index_tag_reads = [x for x in ['read_1', 'read_2'] if x in amplicon_struct_params]
#    assert all([x in sample_table_lane.columns for x in index_tag_cols]), 'Not all specified index_tag columns are in the sample table.\nRequired columns: {}\nSample table: {}'.format(index_tag_cols, config_params['sample_table_file'])
#
#    # Check for duplicates in the index tags
#    num_dup_index_tags = sum(sample_table_lane[index_tag_cols].duplicated())
#    assert num_dup_index_tags == 0, '{} duplicate index tags within lane "{}" detected. Please review the sample table, found here: {}'.format(num_dup_index_tags, lane_id, config_params['sample_table_file'])
#
#    # If this is a single read experiment, add a dummy column of
#    # empty strings for the other read's empty index tag.
#    if len(index_tag_cols) == 1:
#        sample_table_lane['dummy_index_tag'] = ['' for i in range(sample_table_lane.shape[0])]
#        if index_tag_reads == ['read_1']:
#            index_tag_cols.append('dummy_index_tag')
#        else:
#            index_tag_cols.insert(0, 'dummy_index_tag')
#    
#    # Convert the table to a dictionary from index tag(s) to screen_name
#    # and expt_id.
#    # Note: the key is always a length-2 tuple. If there is only one
#    # index tag, the entry for the other index tag is just empty strings.
#    sample_table_lane['screen_name_and_expt_id'] = [tuple(x[1]) for x in sample_table_lane[['screen_name', 'expt_id']].iterrows()]
#    index_tag_to_condition = sample_table_lane.set_index(index_tag_cols)['screen_name_and_expt_id'].to_dict()
#
#    return index_tag_to_condition
#    
#def fastq_to_barseq(config_params, lane_id):
#
#    lane_location_tab = get_lane_location_table(config_params)
#    
#    # Get the raw data folder 
#    raw_folder = get_lane_folder(lane_id, lane_location_tab)
#
#    # Set up the output folders for the temporary barseq file, the output count matrix,
#    # and the read statistics
#    tmp_data_path = get_lane_data_paths(config_params, lane_id)[0]
#    data_path = get_lane_data_paths(config_params, lane_id)[1]
#    reports_path = get_lane_reports_path(config_params, lane_id)
#    cg_file.create_output_dir(tmp_data_path)
#    cg_file.create_output_dir(data_path)
#    cg_file.create_output_dir(reports_path)
#
#    total_counts, common_primer_counts, barcodes_in_data, index_tags_in_data, read_type_dict, parsed_coords = read_fastq(config_params, raw_folder, data_path, lane_id)
#
#    return total_counts, common_primer_counts, barcodes_in_data, index_tags_in_data, read_type_dict, parsed_coords


def is_fastq_filename(filename):
    '''
    Checks to see if the last extension is "fastq" or if the
    second-to-last extension is "fastq" with an appropriate extension.
    Designed to accommodate compressed files and exclude md5 checksums.
    If files are doubly-compressed for some reason, this will not read
    them.
    '''

    ext_split_f = filename.rsplit('.', 2)
    split_len = len(ext_split_f)
    # If there is no extension, ignore the file
    if ext_split_f[0] == filename:
        return False
    # Looks for .fastq.gz and .fastq.bz
    elif ext_split_f[split_len - 2].lower() == 'fastq':
        if ext_split_f[split_len - 1].lower() in ['gz', 'bz']:
            return True
    # Looks for only .fastq
    elif ext_split_f[split_len - 1].lower() == 'fastq':
        return True
    else:
        return False

def get_seq_params(amplicon_struct_params, read):
    try:
        common_primer_start = int(amplicon_struct_params[read]['common_primer']['start'])
        common_primer_seq = amplicon_struct_params[read]['common_primer']['sequence']
        common_primer_end = common_primer_start + len(common_primer_seq)
    except (KeyError, TypeError):
        common_primer_start = -1
        common_primer_seq = ''
        common_primer_end = -1
    try:
        index_tag_start = int(amplicon_struct_params[read]['index_tag']['start'])
        index_tag_end = index_tag_start + int(amplicon_struct_params[read]['index_tag']['length'])
    except (KeyError, TypeError):
        index_tag_start = -1
        index_tag_end = -1
    try:
        barcode_start = int(amplicon_struct_params[read]['genetic_barcode']['start'])
        barcode_end = barcode_start + int(amplicon_struct_params[read]['genetic_barcode']['length'])
    except (KeyError, TypeError):
        barcode_start = -1
        barcode_end = -1

    #return [common_primer_start, common_primer_seq, common_primer_end, index_tag_start, index_tag_end, barcode_start, barcode_end]
    return {'common_primer' : { 'start' : common_primer_start , 'end' : common_primer_end, 'seq' : common_primer_seq },
            'index_tag' : { 'start' : index_tag_start, 'end' : index_tag_end },
            'barcode' : {'start' : barcode_start, 'end' : barcode_end }
            }

def determine_read_type(read_1_params, read_2_params):
    #base_case = [-1, '', -1, -1, -1, -1, -1]
    base_case =  {'common_primer' : { 'start' : -1 , 'end' : -1, 'seq' : '' },
                  'index_tag' : { 'start' : -1, 'end' : -1 },
                  'barcode' : {'start' : -1, 'end' : -1 }
                  }
    #required_read_1_present = all([x is not base_case[i] for i,x in enumerate(read_1_params) if i in range(5)])
    #required_read_2_present = all([x is not base_case[i] for i,x in enumerate(read_2_params) if i in range(5)])
    required_read_1_list = []
    required_read_2_list = []
    for x in ['common_primer', 'index_tag']:
        for y in base_case[x].keys():
            required_read_1_list.append(read_1_params[x][y] is not base_case[x][y])
            required_read_2_list.append(read_2_params[x][y] is not base_case[x][y])
    #required_read_1_present = all([read_1_params[x][y] is not base_case[x][y] for y in [base_case[x].keys() for x in ['common_primer', 'index_tag']]])
    #required_read_2_present = all([read_2_params[x][y] is not base_case[x][y] for y in [base_case[x].keys() for x in ['common_primer', 'index_tag']]])
    #barcode_read_1_present = all([x is not base_case[i] for i,x in enumerate(read_1_params) if i in [5, 6]])
    #barcode_read_2_present = all([x is not base_case[i] for i,x in enumerate(read_2_params) if i in [5, 6]])
    barcode_read_1_list = []
    barcode_read_2_list = []
    for y in base_case['barcode'].keys():
        barcode_read_1_list.append(read_1_params['barcode'][y] is not base_case['barcode'][y])
        barcode_read_2_list.append(read_2_params['barcode'][y] is not base_case['barcode'][y])
    #barcode_read_1_present = all([read_1_params[x][y] is not base_case[x][y] for y in base_case[x].keys() for x in ['barcode']])
    #barcode_read_2_present = all([read_2_params[x][y] is not base_case[x][y] for y in base_case[x].keys() for x in ['barcode']])
    #print required_read_1_present, barcode_read_1_present, required_read_2_present, barcode_read_2_present
    required_read_1_present = all(required_read_1_list)
    required_read_2_present = all(required_read_2_list)
    barcode_read_1_present = all(barcode_read_1_list)
    barcode_read_2_present = all(barcode_read_2_list)

    if required_read_1_present and barcode_read_1_present and not required_read_2_present and not barcode_read_2_present:
        return {'type': 'single', 'barcode': ['read_1']}, [read_1_params]
    elif not required_read_1_present and not barcode_read_1_present and required_read_2_present and barcode_read_2_present:
        return {'type': 'single', 'barcode': ['read_2']}, [read_2_params]
    elif required_read_1_present and barcode_read_1_present and required_read_2_present and not barcode_read_2_present:
        return {'type': 'paired', 'barcode': ['read_1']}, [read_1_params, read_2_params]
    elif required_read_1_present and not barcode_read_1_present and required_read_2_present and barcode_read_2_present:
        return {'type': 'paired', 'barcode': ['read_2']}, [read_1_params, read_2_params]
    elif required_read_1_present and barcode_read_1_present and required_read_2_present and barcode_read_2_present:
        return {'type': 'paired', 'barcode': ['read_1', 'read_2']}, [read_1_params, read_2_params]
        pass # paired end, both barcode
    assert False, 'Required "read_1" and/or "read_2" parameters were not supplied. Please check the "amplicon_struct_file"'

def get_fastq_filename_list(folder, read_type, barcode_reads):
    if read_type == 'single' and barcode_reads == ['read_1']:
        filenames = [[os.path.join(folder, x), None] for x in os.listdir(folder) if is_fastq_filename(x)]
        filenames.sort()
    elif read_type == 'single' and barcode_reads == ['read_2']:
        filenames = [[None, os.path.join(folder, x)] for x in os.listdir(folder) if is_fastq_filename(x)]
        filenames.sort()
    else:
        base_filenames = [x for x in os.listdir(folder) if is_fastq_filename(x)]
        # Sorting ensures the files are in "R1", "R2" order
        base_filenames.sort()
        # While sorting the filenames would be enough to pair them up, it is
        # not particularly robust to missing files and such. So we perform a
        # more explicit approach here.
        no_r1_r2_filenames = [x.replace('_R1', '_R0').replace('_R2', '_R0') for x in base_filenames]
        no_r1_r2_filenames_dict = {}
        for i, no_r1_r2_filename in enumerate(no_r1_r2_filenames):
            if no_r1_r2_filename in no_r1_r2_filenames_dict:
                no_r1_r2_filenames_dict[no_r1_r2_filename].append(base_filenames[i])
            else:
                no_r1_r2_filenames_dict[no_r1_r2_filename] = [base_filenames[i]]
        assert len(base_filenames) == 2 * len(no_r1_r2_filenames_dict), 'Did not detect an equal number of "R1" and "R2" files in this raw data directory:\n{}'.format(folder)
        filenames = [[os.path.join(folder, y) for y in x] for x in no_r1_r2_filenames_dict.values()]

    assert len(filenames) > 0, '\nDid not detect any fastq files for lane: {}.\n' \
            'Please move/copy/symlink your raw data to this folder:\n{}\n'.format(os.path.basename(folder), folder)

    return filenames


#def read_fastq(config_params, folder, out_path, lane_id):
#
#    amplicon_struct_params = get_amplicon_struct_params(config_params)
#
#    # Get all possible common primer/index tag/barcode parameters, then determine
#    # how to proceed.
#    read_1_params = get_seq_params(amplicon_struct_params, 'read_1')
#    read_2_params = get_seq_params(amplicon_struct_params, 'read_2')
#
#    read_type_dict = determine_read_type(read_1_params, read_2_params)
#
#    common_primer_tolerance = int(config_params['common_primer_tolerance'])
#    
#    out_filename = get_barseq_filename(config_params, lane_id)
#    # print out_filename
#    # print common_primer_tolerance
#
#    of = open(out_filename, 'wt')
#
#    # Need a function to return pairs of filenames, for which a None object
#    # is one of the filenames if this is a single read run
#    #fastq_filenames = [os.path.join(folder, x) for x in os.listdir(folder) if is_fastq_filename(x)]
#    fastq_filenames = get_fastq_filename_list(folder, read_type_dict['type'], read_type_dict['barcode'])
#
#    common_primer_count = 0
#    index_tags = set()
#    barcodes = set()
#    for filename_pair in fastq_filenames:
#        handles = [cfo.get_compressed_file_handle(x) if x is not None else it.cycle('\n') for x in filename_pair]
#        line_count = 0
#        while True:
#        #for line_count, line in enumerate(f):
#            try:
#                line_1 = handles[0].next()
#                line_2 = handles[1].next()
#            except StopIteration:
#                break
#            if get_verbosity(config_params) >= 3:
#                print line_count, line_1, line_2
#            if line_count % 4 == 1:
#                #string = line.strip()
#                common_primer_1 = line_1[read_1_params[0]:read_1_params[2]]
#                common_primer_dist_1 = jf.hamming_distance(unicode(common_primer_1), unicode(read_1_params[1]))
#                common_primer_2 = line_2[read_2_params[0]:read_2_params[2]]
#                common_primer_dist_2 = jf.hamming_distance(unicode(common_primer_2), unicode(read_2_params[1]))
#                if get_verbosity(config_params) >= 3:
#                    print common_primer_1, read_1_params[1], common_primer_dist_1
#                    print common_primer_2, read_2_params[1], common_primer_dist_2
#                if common_primer_dist_1 <= common_primer_tolerance and common_primer_dist_2 <= common_primer_tolerance:
#                    common_primer_count += 1
#                    index_tag_1 = line_1[read_1_params[3]:read_1_params[4]]
#                    barcode_1 = line_1[read_1_params[5]:read_1_params[6]]
#                    index_tag_2 = line_2[read_2_params[3]:read_2_params[4]]
#                    barcode_2 = line_2[read_2_params[5]:read_2_params[6]]
#                    index_tags.add((index_tag_1, index_tag_2))
#                    barcodes.add((barcode_1, barcode_2))
#                    # print "index_tag, barcode : {}, {}".format(index_tag, barcode)
#                    of.write('{0}\t{1}\t{2}\t{3}\n'.format(index_tag_1, barcode_1, index_tag_2, barcode_2))
#            line_count += 1
#                
#        for h in handles:
#            if hasattr(h, 'read'):
#                h.close()
#
#    total_counts = (line_count + 1)/ 4
#
#    of.close()
#
#    parsed_lengths = [read_1_params[4] - read_1_params[3], read_1_params[6] - read_1_params[5],
#            read_2_params[4] - read_2_params[3], read_2_params[6] - read_2_params[5]]
#
#    shifts = range(4)
#    parsed_coords = []
#    for i, length in enumerate(parsed_lengths):
#        if i == 0:
#            parsed_coords.append(0)
#        else:
#            parsed_coords.append(parsed_coords[-1] + 1)
#        parsed_coords.append(parsed_coords[-1] + parsed_lengths[i])
#
#    return total_counts, common_primer_count, barcodes, index_tags, read_type_dict, parsed_coords

#def correct_barcode_map(config_params, barcodes_in_data, barcode_to_gene):
#
#    barcode_tolerance = int(config_params['barcode_tolerance'])
#
#    barcode_correcting_maps = [None for i in (0, 1)]
#    for i in (0, 1):
#        orig_barcodes = set([x[i] for x in barcode_to_gene.keys()])
#        observed_barcodes = set([x[i] for x in barcodes_in_data])
#        full_map = {x:x for x in orig_barcodes}
#        correcting_map = {}
#
#        unmatched_barcodes = observed_barcodes.difference(orig_barcodes)
#        for orig_barcode in orig_barcodes:
#            for unmatched_barcode in unmatched_barcodes:
#                barcode_dist = jf.hamming_distance(unicode(orig_barcode), unicode(unmatched_barcode))
#                if barcode_dist <= barcode_tolerance:
#                    if get_verbosity(config_params) >= 3:
#                        print 'bad : corrected --> {0} : {1}'.format(unmatched_barcode, orig_barcode)
#                    if correcting_map.has_key(unmatched_barcode):
#                        correcting_map[unmatched_barcode].append(orig_barcode)
#                    else:
#                        correcting_map[unmatched_barcode] = [orig_barcode]
#
#        # Now, filter out any unmatched barcodes that map to multiple original barcodes
#        for key in correcting_map.keys():
#            if len(correcting_map[key]) > 1:
#                correcting_map[key].pop()
#
#        # The corrected barcodes are still lists - turn them back into strings!
#        corrected_barcodes = correcting_map.keys()
#        for barcode in corrected_barcodes:
#            correcting_map[barcode] = correcting_map[barcode][0]
#
#        # Update the mapping of original barcodes to themselves with the mapping of
#        # unmatched barcodes to original barcodes
#        full_map.update(correcting_map)
#        barcode_correcting_maps[i] = full_map
#    
#    return barcode_correcting_maps


#def get_barseq_matrix(config_params, lane_id, parsed_coords, barcode_to_gene, barcode_correcting_map, index_tag_to_condition):
#
#    # Iterates over barseq text file to fill the matrix
#
#    # Set up the matrix
#    barcodes = barcode_to_gene.keys()
#    bc_ind = {j:i for i,j in enumerate(barcodes)}
#
#    index_tags = index_tag_to_condition.keys()
#    tag_ind = {j:i for i,j in enumerate(index_tags)}
#    
#    matrix = np.zeros([len(barcodes), len(index_tags)], dtype = np.int)
#
#    # Open the barseq file
#    barseq_filename = get_barseq_filename(config_params, lane_id)
#    f = open(barseq_filename, 'rt')
#
#    for line in f:
#        index_tag = (line[parsed_coords[0]:parsed_coords[1]], line[parsed_coords[4]:parsed_coords[5]])
#        barcode = (line[parsed_coords[2]:parsed_coords[3]], line[parsed_coords[6]:parsed_coords[7]])
#        #print index_tag
#        #print barcode
#        
#        try:
#            # Since each barcode is now a 2-tuple, must correct each one
#            # individually (might be a bit slower, but much less complex
#            # than building a dictionary of every possible 2-tuple variation.
#            barcode_fixed = (barcode_correcting_map[0][barcode[0]], barcode_correcting_map[1][barcode[1]])
#            #print barcode_fixed
#            row = bc_ind[barcode_fixed]
#            col = tag_ind[index_tag]
#        except KeyError:
#            continue
#        
#        matrix[row, col] += 1
#
#    # Here is where I implement the "dumb filter," as Justin called it, to
#    # automatically remove any strains or conditions that have zero counts.
#    # Raamesh's version originally did this with a threshold of >0, and I
#    # removed it for some reason. However, this filter was not necessarily
#    # stringent enough, as I accidentally told someone to test the code using
#    # the full genome barcode file instead of the minipool file, leading to
#    # ~183 strains that averaged less then 0.4 counts per condition and
#    # exposing a bug in the new python lowess implementation. I will set this
#    # so that each strain must average >= 1 count in each condition, for added
#    # stringency.
#    mean_thresh = 1
#    sufficient_count_barcodes = np.mean(matrix, axis = 1) > mean_thresh
#    sufficient_count_tags = np.mean(matrix, axis = 0) > mean_thresh
#
#    if get_verbosity(config_params) >= 1:
#        print "number of barcodes with mean counts per condition < {}: {}".format(mean_thresh, sum(np.invert(sufficient_count_barcodes)))
#        print "number of index tags with mean counts per strain < {}: {}".format(mean_thresh, sum(np.invert(sufficient_count_tags)))
#    
#    filtered_matrix = matrix[np.ix_(sufficient_count_barcodes, sufficient_count_tags)]
#    filtered_barcodes = [x for i,x in enumerate(barcodes) if sufficient_count_barcodes[i]]
#    filtered_index_tags = [x for i,x in enumerate(index_tags) if sufficient_count_tags[i]]
#
#    assert len(barcodes) == matrix.shape[0], "number of barcodes does not match number of rows in matrix after sufficient-count filtering."
#    assert len(index_tags) == matrix.shape[1], "number of index tags does not match number of columns in matrix after sufficient-count filtering."
#
#    f.close()
#
#    return filtered_barcodes, filtered_index_tags, filtered_matrix, matrix

def get_sorted_counts(label, counts):

    label = np.array(label)
    
    # Take "no_match" and "multi_match" into account
    no_or_multi_label = np.array(['multi_match', 'no_match'])
    no_or_multi_counts = np.array([counts[label == x] for x in no_or_multi_label])
    
    no_or_multi_bool = (label == 'no_match') | (label == 'multi_match')
    label = label[~no_or_multi_bool]
    counts = counts[~no_or_multi_bool]

    # Get in descending order - the [::-1] reverses the array
    sort_indices = counts.argsort()[::-1]
    labels_sorted = label[sort_indices]
    counts_sorted = counts[sort_indices]

    labels_sorted = np.append(labels_sorted, no_or_multi_label)
    counts_sorted = np.append(counts_sorted, no_or_multi_counts)

    return labels_sorted, counts_sorted

def write_distribution(labels, counts, output_folder, lane_id, file_string, table_label, xlabel):

    plot_filename = os.path.join(output_folder, '{0}_{1}_distribution.png'.format(lane_id, file_string))
    x = np.array(range(len(counts)))

    plt.figure()
    plt.fill_between(x, y1 = 0, y2 = counts, step = 'post')
    plt.ylabel('Number of occurrences')
    plt.xlabel(xlabel)
    plt.savefig(plot_filename)

    text_filename = os.path.join(output_folder, '{0}_{1}_distribution.txt'.format(lane_id, file_string))
    of = open(text_filename, 'wt')
    of.write('{0}\tNumber of occurrences\n'.format(table_label))
    for i in range(len(labels)):
        of.write('{0}\t{1}\n'.format(labels[i], counts[i]))
    of.close()

def write_summary(total_counts, common_primer_counts, total_index_barcode_counts, matrix_shape, orig_matrix_shape, output_folder, lane_id):

    filename = os.path.join(output_folder, '{0}_summary.txt'.format(lane_id))

    common_primer_vs_total_percent = float(common_primer_counts) / float(total_counts) * 100
    tag_barcode_match_vs_total_percent = float(total_index_barcode_counts) / float(total_counts) * 100
    tag_barcode_match_vs_common_primer_percent = float(total_index_barcode_counts) / float(common_primer_counts) * 100

    of = open(filename, 'wt')
    of.write('Total number of reads: {0}\n'.format(total_counts))
    of.write('Number of reads with common primer: {0} ({1:.2f} %)\n'.format(common_primer_counts, common_primer_vs_total_percent))
    of.write('Number of reads that match index tags and genetic barcodes: {0} ({1:.2f} % of total counts, {2:.2f} % of common primer counts)\n'.format(total_index_barcode_counts, tag_barcode_match_vs_total_percent, tag_barcode_match_vs_common_primer_percent))

    # Write original and final matrix sizes
    of.write('Original size of the data: {0} barcodes (gene mutants) x {1} index tags (conditions)\n'.format(*orig_matrix_shape))
    of.write('Final size of the data (zero-count barcodes/tags removed): {0} barcodes (gene mutants) x {1} index tags (conditions)\n'.format(*matrix_shape))

    of.close()

def generate_reports(config_params, lane_id, gen_barcodes, index_tags, matrix, orig_barcode_map, orig_index_tag_map, orig_matrix, total_counts, common_primer_counts):

    output_folder = get_lane_reports_path(config_params, lane_id)

    gen_barcode_counts = orig_matrix.sum(axis = 1)
    index_tag_counts = orig_matrix.sum(axis = 0)
    total_index_barcode_counts = orig_matrix.sum()
    orig_matrix_shape = (len(orig_barcode_map), len(orig_index_tag_map))

    index_tags_sorted, index_tag_counts_sorted = get_sorted_counts(orig_index_tag_map.keys(), index_tag_counts)
    gen_barcodes_sorted, gen_barcode_counts_sorted = get_sorted_counts(orig_barcode_map.keys(), gen_barcode_counts)

    write_distribution(index_tags_sorted, index_tag_counts_sorted, output_folder, lane_id, 'index_tag', 'Index tags (sorted)')
    write_distribution(gen_barcodes_sorted, gen_barcode_counts_sorted, output_folder, lane_id, 'barcode', 'Genetic barcodes (sorted)')

    write_summary(total_counts, common_primer_counts, total_index_barcode_counts, matrix.shape, orig_matrix_shape, output_folder, lane_id)

def generate_reports(config_params, lane_id, count_array_dataset, count_matrix_dataset, filtered_count_matrix_dataset, seq_info_dict, total_counts, common_primer_counts):

    # seq_info_dict contains: match_dicts, read_ids, seq_types, and column_names

    # Get and make output folder
    output_folder = get_lane_reports_path(config_params, lane_id)
    if not os.path.isdir(output_folder):
        os.makedirs(output_folder)

    # Get total counts for each index tag/barcode
    count_distributions = [None] * len(count_array_dataset[1].shape)
    for i in range(len(count_distributions)):
        sum_axes = tuple(j for j in range(len(count_distributions)) if j != i)
        count_distributions[i] = count_array_dataset[1].sum(axis = sum_axes)

    # Get total counts for each strain/condition
    strain_counts = count_matrix_dataset[2].sum(axis = 1)
    condition_counts = count_matrix_dataset[2].sum(axis = 0)

    # Get total counts for read/read pairs that successfully mapped from all barcodes/index tags
    mapped_read_counts = count_array_dataset[1][tuple(slice(-2) for x in count_array_dataset[1].shape)].sum()

    # Plot total counts for individual index tags/barcodes (minus the "no_match" and "multi_match" instances - those only go into the tables)
    label_dict = {'index_tag': 'Index tags (sorted)', 'barcode': 'Genetic barcodes (sorted)'} 
    for i, read_id in enumerate(seq_info_dict['read_ids']):
        seq_type = seq_info_dict['seq_types'][i]
        seq_table_label = seq_info_dict['column_names'][i]
        seqs_sorted, counts_sorted = get_sorted_counts(count_array_dataset[0][i], count_distributions[i])
        write_distribution(seqs_sorted, counts_sorted, output_folder, lane_id, '{}_{}'.format(seq_type, read_id), seq_table_label, label_dict[seq_type]) 

    # Plot total counts for strains and conditions
    strains_sorted, strain_counts_sorted = get_sorted_counts(count_matrix_dataset[0], strain_counts)
    conditions_sorted, condition_counts_sorted = get_sorted_counts(count_matrix_dataset[1], condition_counts)
    # write_distribution needs to be able to take in extra info columns!

    pdb.set_trace()

    gen_barcode_counts = orig_matrix.sum(axis = 1)
    index_tag_counts = orig_matrix.sum(axis = 0)
    total_index_barcode_counts = orig_matrix.sum()
    orig_matrix_shape = (len(orig_barcode_map), len(orig_index_tag_map))

    index_tags_sorted, index_tag_counts_sorted = get_sorted_counts(orig_index_tag_map.keys(), index_tag_counts)
    gen_barcodes_sorted, gen_barcode_counts_sorted = get_sorted_counts(orig_barcode_map.keys(), gen_barcode_counts)

    write_distribution(index_tags_sorted, index_tag_counts_sorted, output_folder, lane_id, 'index_tag', 'Index tags (sorted)')
    write_distribution(gen_barcodes_sorted, gen_barcode_counts_sorted, output_folder, lane_id, 'barcode', 'Genetic barcodes (sorted)')

    write_summary(total_counts, common_primer_counts, total_index_barcode_counts, matrix.shape, orig_matrix_shape, output_folder, lane_id)

def dump_count_matrix(config_params, lane_id, barcodes, conditions, matrix):

    out_path = get_lane_data_paths(config_params, lane_id)[1]
    out_filename = os.path.join(out_path, '{0}_{1}'.format(lane_id, 'barseq_matrix.dump.gz'))

    of = gzip.open(out_filename, 'wb')

    dataset = [barcodes, conditions, matrix]
    cPickle.dump(dataset, of)

    of.close()

def get_fastq_filename_list(folder, read_type):
    # Change from first implementation: just return a list of filenames if single read,
    # and pairs of filenames if paired end read. The software can figure it out form there.
    if read_type == 'single':
        filenames = [os.path.join(folder, x) for x in os.listdir(folder) if is_fastq_filename(x)]
    else:
        base_filenames = [x for x in os.listdir(folder) if is_fastq_filename(x)]
        # Sorting ensures the files are in "R1", "R2" order
        base_filenames.sort()
        # While sorting the filenames would be enough to pair them up, it is
        # not particularly robust to missing files and such. So we perform a
        # more explicit approach here.
        no_r1_r2_filenames = [x.replace('_R1', '_R0').replace('_R2', '_R0') for x in base_filenames]
        no_r1_r2_filenames_dict = {}
        for i, no_r1_r2_filename in enumerate(no_r1_r2_filenames):
            if no_r1_r2_filename in no_r1_r2_filenames_dict:
                no_r1_r2_filenames_dict[no_r1_r2_filename].append(base_filenames[i])
            else:
                no_r1_r2_filenames_dict[no_r1_r2_filename] = [base_filenames[i]]
        assert len(base_filenames) == 2 * len(no_r1_r2_filenames_dict), 'Did not detect an equal number of "R1" and "R2" files in this raw data directory:\n{}'.format(folder)
        filenames = [[os.path.join(folder, y) for y in x] for x in no_r1_r2_filenames_dict.values()]

    assert len(filenames) > 0, '\nDid not detect any fastq files for lane: {}.\n' \
            'Please move/copy/symlink your raw data to this folder:\n{}\n'.format(os.path.basename(folder), folder)

    return filenames

def line_gen(folder, read_type):

    filenames = get_fastq_filename_list(folder, read_type)

    if read_type == 'single':
        for fname in filenames:
            with cfo.get_compressed_file_handle(fname) as f:
                for i, line in enumerate(f):
                    if i % 4 == 1:
                        yield [line.rstrip()]
    else:
        for fname_pair in filenames:
            with cfo.get_compressed_file_handle(fname_pair[0]) as f1, cfo.get_compressed_file_handle(fname_pair[1]) as f2:
                for i, (line1, line2) in enumerate(it.izip(f1, f2)):
                    if i % 4 == 1:
                        yield [line1.rstrip(), line2.rstrip()]

def gen_seq_tries(seqs):
    lengths = list(set([len(seq) for seq in seqs]))
    lengths.sort(key = lambda x: -x)
    trie_list = [trie() for i in range(len(lengths))]
    len_dict = {x:i for i,x in enumerate(lengths)}
    for seq in seqs:
        trie_list[len_dict[len(seq)]][seq] = 0

    return trie_list, lengths

def initialize_dicts_arrays(read_type_dict, amplicon_struct_params, config_params, sample_tab, barcode_tab):
    # Need to add in return of sequence tries, valid lengths, and mismatches
    '''
    read_inds, seq_types, match_dicts, array_ind_dicts, seq_trie_lists, seq_lengths, tols, column_names, array = initialize_dicts_arrays(...)
    '''
    if read_type_dict['type'] == 'single':
        index_tag_col = amplicon_struct_params[read_type_dict['barcode'][0]]['index_tag']['sample_table_column']
        assert index_tag_col in sample_tab.columns, 'sample_table_column "{}" specified in amplicon_struct_file is not present in the sample table.'.format(index_tag_col)
        index_tag_tries, index_tag_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col])
        index_tags = {x:i for i, x in enumerate(sample_tab.loc[:, index_tag_col])}
        index_tags['multi_match'] = len(index_tags)
        index_tags['no_match'] = len(index_tags)
        barcode_col = amplicon_struct_params[read_type_dict['barcode'][0]]['genetic_barcode']['barcode_file_column']
        assert barcode_col in barcode_tab.columns, 'barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col)
        barcode_tries, barcode_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col])
        barcodes = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col])}
        barcodes['multi_match'] = len(barcodes)
        barcodes['no_match'] = len(barcodes)
        index_tag_tol = config_params.get('index_tag_tolerance', 0)
        barcode_tol = config_params.get('barcode_tolerance', 0)
        count_array = np.zeros((len(index_tags), len(barcodes)), dtype = np.int)
        # While the scheme is the same whether or not the single read is read_1
        # or read_2 (would it actually ever be read_2?!), I will return the id
        # of the read for use later on
        if read_type_dict['barcode'] == ['read_1']:
            read_id = 'read_1'
        elif read_type_dict['barcode'] == ['read_2']:
            read_id = 'read_2'
        else:
            assert False, "{} is not a valid read id".format(read_type_dict['barcode'])
        return [0, 0], ['index_tag', 'barcode'], [read_id, read_id], [{}, {}], [index_tags, barcodes], [index_tag_tries, barcode_tries], [index_tag_lengths, barcode_lengths], [index_tag_tol, barcode_tol], [index_tag_col, barcode_col], count_array
    elif read_type_dict['type'] == 'paired':
        # Since there are always 2 index tags in a paired read scheme, I can write this code once.
        index_tag_col_1 = amplicon_struct_params['read_1']['index_tag']['sample_table_column']
        assert index_tag_col_1 in sample_tab.columns, 'read_1 sample_table_column "{}" specified in amplicon_struct_file is not present in the sample table.'.format(index_tag_col)
        index_tag_col_2 = amplicon_struct_params['read_2']['index_tag']['sample_table_column']
        assert index_tag_col_2 in sample_tab.columns, 'read_2 sample_table_column "{}" specified in amplicon_struct_file is not present in the sample table.'.format(index_tag_col)
        if read_type_dict['barcode'] == ['read_1']:
            index_tag_1_tries, index_tag_1_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_1])
            index_tags_1 = {x:i for i, x in enumerate(sample_tab.loc[:, index_tag_col_1])}
            index_tags_1['multi_match'] = len(index_tags_1)
            index_tags_1['no_match'] = len(index_tags_1)
            barcode_col_1 = amplicon_struct_params['read_1']['genetic_barcode']['barcode_file_column']
            assert barcode_col_1 in barcode_tab.columns, 'read_1 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col_1)
            barcode_1_tries, barcode_1_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col_1])
            barcodes_1 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_1])}
            barcodes_1['multi_match'] = len(barcodes_1)
            barcodes_1['no_match'] = len(barcodes_1)
            index_tag_2_tries, index_tag_2_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_2])
            index_tags_2 = {x:i for i, x in enumerate(sample_tab.loc[:, index_tag_col_2])}
            index_tags_2['multi_match'] = len(index_tags_2)
            index_tags_2['no_match'] = len(index_tags_2)
            index_tag_tol = config_params.get('index_tag_tolerance', 0)
            barcode_tol = config_params.get('barcode_tolerance', 0)
            count_array = np.zeros((len(index_tags_1), len(barcodes_1), len(index_tags_2)), dtype = np.int)
            return [0, 0, 1], ['index_tag', 'barcode', 'index_tag'], ['read_1', 'read_1', 'read_2'], [{}, {}, {}], [index_tags_1, barcodes_1, index_tags_2], [index_tag_1_tries, barcode_1_tries, index_tag_2_tries], [index_tag_1_lengths, barcode_1_lengths, index_tag_2_lengths], [index_tag_tol, barcode_tol, index_tag_tol], [index_tag_col_1, barcode_col_1, index_tag_col_2], count_array
        elif read_type_dict['barcode'] == ['read_2']:
            index_tag_1_tries, index_tag_1_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_1])
            index_tags_1 = {x:i for i, x in enumerate(sample_tab.loc[:, index_tag_col_1])}
            index_tags_1['multi_match'] = len(index_tags_1)
            index_tags_1['no_match'] = len(index_tags_1)
            index_tag_2_tries, index_tag_2_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_2])
            index_tags_2 = {x:i for i, x in enumerate(sample_tab.loc[:, index_tag_col_2])}
            index_tags_2['multi_match'] = len(index_tags_2)
            index_tags_2['no_match'] = len(index_tags_2)
            barcode_col_2 = amplicon_struct_params['read_2']['genetic_barcode']['barcode_file_column']
            assert barcode_col_2 in barcode_tab.columns, 'read_2 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col_2)
            barcode_2_tries, barcode_2_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col_2])
            barcodes_2 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_2])}
            barcodes_2['multi_match'] = len(barcodes_2)
            barcodes_2['no_match'] = len(barcodes_2)
            index_tag_tol = config_params.get('index_tag_tolerance', 0)
            barcode_tol = config_params.get('barcode_tolerance', 0)
            count_array = np.zeros((len(index_tags_1), len(index_tags_2), len(barcodes_2)), dtype = np.int)
            return [0, 1, 1], ['index_tag', 'index_tag', 'barcode'], ['read_1', 'read_2', 'read_2'], [{}, {}, {}], [index_tags_1, index_tags_2, barcodes_2], [index_tag_1_tries, index_tag_2_tries, barcode_2_tries], [index_tag_1_lengths, index_tag_2_lengths, barcode_2_lengths], [index_tag_tol, index_tag_tol, barcode_tol], [index_tag_col_1, index_tag_col_2, barcode_col_1], count_array
        elif read_type_dict['barcode'] == ['read_1', 'read_2']:
            index_tag_1_tries, index_tag_1_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_1])
            index_tags_1 = {x:i for i, x in enumerate(sample_tab.loc[:, index_tag_col_1])}
            index_tags_1['multi_match'] = len(index_tags_1)
            index_tags_1['no_match'] = len(index_tags_1)
            barcode_col_1 = amplicon_struct_params['read_1']['genetic_barcode']['barcode_file_column']
            assert barcode_col_1 in barcode_tab.columns, 'read_1 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col_1)
            barcode_1_tries, barcode_1_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col_1])
            barcodes_1 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_1])}
            barcodes_1['multi_match'] = len(barcodes_1)
            barcodes_1['no_match'] = len(barcodes_1)
            index_tag_2_tries, index_tag_2_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_2])
            index_tags_2 = {x:i for i, x in enumerate(sample_tab.loc[:, index_tag_col_2])}
            index_tags_2['multi_match'] = len(index_tags_2)
            index_tags_2['no_match'] = len(index_tags_2)
            barcode_col_2 = amplicon_struct_params['read_2']['genetic_barcode']['barcode_file_column']
            assert barcode_col_2 in barcode_tab.columns, 'read_2 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene barcode table.'.format(barcode_col_2)
            barcode_2_tries, barcode_2_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col_2])
            barcodes_2 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_2])}
            barcodes_2['multi_match'] = len(barcodes_2)
            barcodes_2['no_match'] = len(barcodes_2)
            index_tag_tol = config_params.get('index_tag_tolerance', 0)
            barcode_tol = config_params.get('barcode_tolerance', 0)
            count_array = np.zeros((len(index_tags_1), len(barcodes_1), len(index_tags_2), len(barcodes_2)), dtype = np.int)
            return [0, 0, 1, 1], ['index_tag', 'barcode', 'index_tag', 'barcode'], ['read_1', 'read_1', 'read_2', 'read_2'], [{}, {}, {}, {}], [index_tags_1, barcodes_1, index_tags_2, barcodes_2], [index_tag_1_tries, barcode_1_tries, index_tag_2_tries, barcode_2_tries], [index_tag_1_lengths, barcode_1_lengths, index_tag_2_lengths, barcode_2_lengths], [index_tag_tol, barcode_tol, index_tag_tol, barcode_tol], [index_tag_col_1, barcode_col_1, index_tag_col_2, barcode_col_2], count_array
    

def match_seq(seq, seq_trie_length_list, n_mismatch, lengths):
    '''
    Returns first sequence that can be found in the given trie
    that matches the query sequences within the number of given mismatches.
    Comparisons are only performed between same-length sequences, aka the query
    string is sliced to the the same length as each set of reference sequences
    and the comparisons proceed starting with the longest slices first.

    Returns "multi_match" if the query sequence matches two reference sequences
    with the same distance. Returns "no_match" if it goes through everything
    and no match is found.
    '''

    for i, l in enumerate(lengths):
        res = seq_trie_length_list[i].get_approximate(seq[0:l], n_mismatch)
        if len(res) == 0:
            continue

        # If res meets the mismatch criterion, return before examining any others
        # But must account for multi-matches!
        if len(res) > 1:
            res = list(set(res))
            if len(res) > 1:
                res.sort(key = lambda x: x[2])
                if res[0][2] == res[1][2]:
                    return 'multi_match'
                else:
                    return res[0][0]
            else:
                return res[0][0]

        else:
            return res[0][0]

    # If nothing matched at all, return "no_match"!
    return 'no_match'

def initialize_cp_matchers(read_type_dict, amplicon_struct_params, config_params):


    if read_type_dict['type'] == 'single':
        cp_read_id = [amplicon_struct_params[read_type_dict['barcode'][0]]]
        if cp_read_id == 'read_1':
            cp_read_inds = [0]
        else:
            cp_read_inds = [1]
        cp_seq_list = [amplicon_struct_params[read_type_dict['barcode'][0]]['common_primer']['sequence']]
    else:
        cp_read_ids = ['read_1', 'read_2']
        cp_read_inds = [0, 1]
        cp_seq_list = [amplicon_struct_params[read_id]['common_primer']['sequence'] for read_id in ['read_1', 'read_2']]
        
    cp_tol = config_params.get('common_primer_tolerance', 0)
    bases = ['A', 'C', 'G', 'T']
    cp_dicts = [{}]
    for i in range(len(cp_seq_list)):
        cp_seq = cp_seq_list[i]
        for mismatch_locs in it.combinations(range(len(cp_seq)), cp_tol):
            this_seq = list(cp_seq)
            for ind in mismatch_locs:
                this_seq[ind] = bases
            for mutated_seq_list in it.product(*this_seq):
                cp_dicts[i][''.join(mutated_seq_list)] = 0

    return cp_read_inds, cp_dicts


def parse_seqs(lane_id, config_params):

    amplicon_struct_params = get_amplicon_struct_params(config_params)
    sample_tab = get_sample_table(config_params)
    barcode_tab = get_barcode_table(config_params)

    # Filter the sample table and return nothing if no samples exist with the given lane_id
    sample_tab = sample_tab[sample_tab.lane == lane_id]

    # Get all possible common primer/index tag/barcode parameters, then determine
    # how to proceed.
    read_params = [get_seq_params(amplicon_struct_params, 'read_1'), get_seq_params(amplicon_struct_params, 'read_2')]

    read_type_dict, read_params_clean = determine_read_type(read_params[0], read_params[1])

    cp_read_inds, cp_dicts = initialize_cp_matchers(read_type_dict, amplicon_struct_params, config_params)

    read_inds, seq_types, read_ids, match_dicts, array_ind_dicts, seq_trie_lists, seq_lengths, tols, column_names, array = initialize_dicts_arrays(read_type_dict, amplicon_struct_params, config_params, sample_tab, barcode_tab)

    lane_location_tab = get_lane_location_table(config_params)
    folder = get_lane_folder(lane_id, lane_location_tab)

    #pdb.set_trace()

    n = len(read_inds)
    idxs = [None] * n
    for counter, line_list in enumerate(line_gen(folder, read_type_dict['type'])):
    #for line_list in line_gen(folder, read_type_dict['type']):
        
        if counter % 1000000 == 0:
            print counter

        # For testing
        if counter / 1000000 == 2:
            break
        
        # Check first for common primer!
        cp_match = [False] * len(cp_dicts)
        for i in range(len(cp_dicts)):
            start_coord = read_params_clean[i]['common_primer']['start']
            end_coord = start_coord + len(read_params_clean[i]['common_primer']['seq'])
            seq = line_list[i][start_coord:end_coord]
            try:
                cp_dicts[i][seq] += 1
                cp_match[i] = True
            except KeyError as e:
                pass

        if not all(cp_match):
            continue

        # Match on index tags/barcodes
        for i in range(n):
            start_coord = read_params_clean[read_inds[i]][seq_types[i]]['start']
            end_coord = read_params_clean[read_inds[i]][seq_types[i]]['end']
            seq = line_list[read_inds[i]][start_coord:end_coord]

            try:
                corrected_seq = match_dicts[i][seq][0]
                match_dicts[i][seq][1] += 1
            except KeyError as e:
                # match_seq returns the corrected sequence
                match_dicts[i][seq] = [match_seq(seq, seq_trie_lists[i], tols[i], seq_lengths[i]), 1]
                corrected_seq = match_dicts[i][seq][0]

            idxs[i] = array_ind_dicts[i][corrected_seq]

        array[tuple(idxs)] += 1

    return array, array_ind_dicts, match_dicts, cp_dicts, counter, seq_types, column_names, read_ids

def map_counts_to_strains_conditions(array, array_ind_dicts, seq_types, column_names, config_params):

    seq_types = np.array(seq_types)
    column_names = np.array(column_names)

    sample_tab = get_sample_table(config_params)
    barcode_tab = get_barcode_table(config_params)
    
    # Filter the sample table and return nothing if no samples exist with the given lane_id
    sample_tab = sample_tab[sample_tab.lane == lane_id]

    # Make mappings from each strain and condition to their respective barcode and index tag combinations
    strain_to_barcode = {x[1]['Strain_ID']:tuple(x[1][column_names[seq_types == 'barcode']]) for x in barcode_tab.iterrows()}
    condition_to_index_tag = {tuple(x[1][['screen_name', 'expt_id']]):tuple(x[1][column_names[seq_types == 'index_tag']]) for x in sample_tab.iterrows()}
    
    # Make an empty final matrix
    strains = np.array(strain_to_barcode.keys())
    conditions = np.array(condition_to_index_tag.keys())
    matrix = np.zeros((len(strains), len(conditions)), dtype = np.int)

    # Fill in the matrix!
    ndim = len(array.shape)
    for i, strain in enumerate(strains):
        for j, condition in enumerate(conditions):
            idx = [None] * ndim
            index_tag_counter = 0
            barcode_counter = 0
            for k in range(ndim):
                if seq_types[k] == 'index_tag':
                    idx[k] = array_ind_dicts[k][condition_to_index_tag[tuple(condition)][index_tag_counter]]
                    index_tag_counter += 1
                else:
                    idx[k] = array_ind_dicts[k][strain_to_barcode[strain][barcode_counter]]
                    barcode_counter += 1
            matrix[i, j] = array[tuple(idx)]

    return strains, conditions, matrix

def filter_count_matrix(matrix, strains, conditions, config_params):
    # Here is where I implement the "dumb filter," as Justin called it, to
    # automatically remove any strains or conditions that have zero counts.
    # Raamesh's version originally did this with a threshold of >0, and I
    # removed it for some reason. However, this filter was not necessarily
    # stringent enough, as I accidentally told someone to test the code using
    # the full genome barcode file instead of the minipool file, leading to
    # ~183 strains that averaged less then 0.4 counts per condition and
    # exposing a bug in the new python lowess implementation. I will set this
    # so that each strain must average >= 1 count in each condition, for added
    # stringency.
    mean_thresh = 1
    sufficient_count_strains = np.mean(matrix, axis = 1) > mean_thresh
    sufficient_count_conditions = np.mean(matrix, axis = 0) > mean_thresh

    if get_verbosity(config_params) >= 1:
        print "number of strains with mean counts per condition < {}: {}".format(mean_thresh, sum(np.invert(sufficient_count_strains)))
        print "number of conditions with mean counts per strain < {}: {}".format(mean_thresh, sum(np.invert(sufficient_count_conditions)))
    
    filtered_matrix = matrix[np.ix_(sufficient_count_strains, sufficient_count_conditions)]
    filtered_strains = [x for i,x in enumerate(strains) if sufficient_count_strains[i]]
    filtered_conditions = [x for i,x in enumerate(conditions) if sufficient_count_conditions[i]]

    assert len(strains) == matrix.shape[0], "number of strains does not match number of rows in matrix after sufficient-count filtering."
    assert len(conditions) == matrix.shape[1], "number of index tags does not match number of columns in matrix after sufficient-count filtering."

    return filtered_strains, filtered_conditions, filtered_matrix

def main(config_file, lane_id):

    config_params = parse_yaml(config_file)
    amplicon_struct_params = get_amplicon_struct_params(config_params)

    # Check for common primer and put all reads into a master array
    count_array, array_ind_dicts, match_dicts, cp_dicts, total_reads, seq_types, column_names, read_ids = parse_seqs(lane_id, config_params)
    
    #pdb.set_trace()
    
    # Map the counts to actual conditions and strains
    strains, conditions, count_matrix = map_counts_to_strains_conditions(count_array, array_ind_dicts, seq_types, column_names, config_params)
   
    #pdb.set_trace()
    
    # Filter count matrix to exclude strains/conditions
    filtered_strains, filtered_conditions, filtered_count_matrix = filter_count_matrix(count_matrix, strains, conditions, config_params)

    # Generate reports for index tags and barcodes
    if get_verbosity(config_params) >= 1:
        print 'generating reports...'
    # Get sequences that go along each axis of the count_array
    seqs = [None] * len(array_ind_dicts)
    for i, d in enumerate(array_ind_dicts):
        seq_list = list(array_ind_dicts[i].iteritems())
        seq_list.sort(key = lambda x: x[1])
        seqs[i] = [x[0] for x in seq_list]

    pdb.set_trace()

#def generate_reports(config_params, lane_id, count_array_dataset, count_matrix_dataset, filtered_count_matrix_dataset, match_dicts, seq_info_dict, total_counts, common_primer_counts):
    generate_reports(config_params, lane_id,
            [seqs, count_array],
            [strains, conditions, count_matrix],
            [filtered_strains, filtered_conditions, filtered_count_matrix],
            {'match_dicts': match_dicts,
                'read_ids': read_ids,
                'seq_types': seq_types,
                'column_names': column_names},
            total_reads, count_array.sum())
    
    ## Get maps of barcode to barcode_gene (keeps the strains unique/traceable), and index tag to condition
    #if get_verbosity(config_params) >= 1:
    #    print 'creating mappings from barcodes and index tags...'
    #barcode_to_gene = get_barcode_to_gene(config_params)
    #index_tag_to_condition = get_index_tag_to_condition(config_params, lane_id)

    ## Loop over the raw fastq files, write out the "index_tag\tbarcode" file,
    ## and return all encountered index tags and barcodes
    #if get_verbosity(config_params) >= 1:
    #    print 'parsing fastq file(s)...'
    #total_counts, common_primer_counts, barcodes_in_data, index_tags_in_data, read_type_dict, parsed_coords = fastq_to_barseq(config_params, lane_id)

    #if get_verbosity(config_params) >= 1:
    #    print 'barcode to gene map: {}'.format(barcode_to_gene.items()[0:5])
    #    print 'index tag to condition map: {}'.format(index_tag_to_condition.items()[0:5])

    ## Correct the barcodes within the specified error tolerance
    ## (There is no function to correct the index tags - this could easily be
    ## written in later, although we see no need for it)
    #if get_verbosity(config_params) >= 1:
    #    print 'correcting barcodes...'
    #barcode_correcting_map = correct_barcode_map(config_params, barcodes_in_data, barcode_to_gene)
    #if get_verbosity(config_params) >= 1:
    #    print 'number of read_1 barcodes that will be counted: {}'.format(len(barcode_correcting_map[0]))
    #    print 'number of read_2 barcodes that will be counted: {}'.format(len(barcode_correcting_map[1]))

    ## Loop over the barseq file (index_tag_1\tbarcode_1\tindex_tag_2\tbarcode_2\n)
    ## and assemble the matrix of read counts
    #if get_verbosity(config_params) >= 1:
    #    print 'generating barseq matrix...'
    #corrected_barcodes, index_tags, matrix, unfiltered_matrix = get_barseq_matrix(config_params, lane_id, parsed_coords, barcode_to_gene, barcode_correcting_map, index_tag_to_condition)

    #if get_verbosity(config_params) >= 1:
    #    print 'number of barcodes: {}'.format(len(corrected_barcodes))
    #    print 'number of index tags: {}'.format(len(index_tags))
    #    print 'matrix shape: {0} rows x {1} columns'.format(*matrix.shape)

    # Generate reports for index tags and barcodes
    if get_verbosity(config_params) >= 1:
        print 'generating reports...'
    generate_reports(config_params, lane_id, corrected_barcodes, index_tags, matrix, barcode_to_gene, index_tag_to_condition, unfiltered_matrix, total_counts, common_primer_counts)

    # Convert the barcodes to their condition names and unique gene/barcode names
    barcode_gene_ids = np.array([barcode_to_gene[bc] for bc in corrected_barcodes])
    condition_ids = np.array([index_tag_to_condition[tag] for tag in index_tags])

    # Dump out the final count matrix to file - other scripts will read it and turn it into a readable matrix/CDT
    if get_verbosity(config_params) >= 1:
        print 'dumping count matrix...'
    dump_count_matrix(config_params, lane_id, barcode_gene_ids, condition_ids, matrix)

    # Remove the temporary barseq file
    remove_barseq_file(config_params, lane_id)
 

# call: python fastq_to_count_matrix.py <config_file> <lane_id>
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: python fastq_to_count_matrix.py <config_file> <lane_id>'
    else:
        config_file = sys.argv[1]
        lane_id = sys.argv[2]
        main(config_file, lane_id)
