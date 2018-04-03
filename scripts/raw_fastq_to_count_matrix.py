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

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append(os.path.join(barseq_path, 'lib'))

import compressed_file_opener as cfo
import cg_file_tools as cg_file
from cg_common_functions import get_verbosity, get_sample_table, get_amplicon_struct_params, read_barcode_table, parse_yaml

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

def get_barcode_to_gene(config_params):
    # Need to update for paired-end configuration
    #barseq_path = os.getenv('BARSEQ_PATH')

    amplicon_struct_params = get_amplicon_struct_params(config_params)
    
    filename = config_params['gene_barcode_file']
    #full_path = os.path.join(barseq_path, 'data', 'barcodes', filename)
    barcode_tab = read_barcode_table(filename)

    # Get all possible common primer/index tag/barcode parameters, then determine
    # how to proceed.
    read_1_params = get_seq_params(amplicon_struct_params, 'read_1')
    read_2_params = get_seq_params(amplicon_struct_params, 'read_2')
    read_type_dict = determine_read_type(read_1_params, read_2_params)

    # Which barcode columns are required?
    barcode_cols = [amplicon_struct_params[x]['genetic_barcode']['barcode_file_column'] for x in read_type_dict['barcode']]
    assert all([x in barcode_tab.columns for x in barcode_cols]), 'Not all specified barcode columns are in the barcode table.\nRequired columns: {}\nBarcode table: {}'.format(barcode_cols, filename)

    # Check for duplicates in the barcodes
    num_dup_barcodes = sum(barcode_tab[barcode_cols].duplicated())
    assert num_dup_barcodes == 0, '{} duplicate barcodes detected. Please review the barcode to strain mapping, found here: {}'.format(num_dup_barcodes, filename)

    # If there is only one read containing the barcode, add a dummy
    # column of empty strings for the other read.
    if len(read_type_dict['barcode']) == 1:
        barcode_tab['dummy_barcode'] = ['' for i in range(barcode_tab.shape[0])]
        if read_type_dict['barcode'] == ['read_1']:
            barcode_cols.append('dummy_barcode')
        else:
            barcode_cols.insert(0, 'dummy_barcode')

    # Convert the table to a dictionary from barcode(s) to strain_id
    # Note: the key is always a length-2 tuple. If there is only one
    # barcode, the entry for the other barcode is just empty strings.
    barcode_to_strain = barcode_tab.set_index(barcode_cols)['Strain_ID'].to_dict()

    return barcode_to_strain

def get_index_tag_to_condition(config_params, lane_id):

    sample_table = get_sample_table(config_params)
    sample_table = sample_table.set_index('lane')
    sample_table_lane = sample_table.loc[lane_id]

    # Get all possible common primer/index tag/barcode parameters, then determine
    # how to proceed.
    amplicon_struct_params = get_amplicon_struct_params(config_params)
    read_1_params = get_seq_params(amplicon_struct_params, 'read_1')
    read_2_params = get_seq_params(amplicon_struct_params, 'read_2')
    read_type_dict = determine_read_type(read_1_params, read_2_params)

    # Which index tag columns are required?
    index_tag_cols = [amplicon_struct_params[x]['index_tag']['sample_table_column'] for x in ['read_1', 'read_2'] if x in amplicon_struct_params]
    index_tag_reads = [x for x in ['read_1', 'read_2'] if x in amplicon_struct_params]
    assert all([x in sample_table_lane.columns for x in index_tag_cols]), 'Not all specified index_tag columns are in the sample table.\nRequired columns: {}\nSample table: {}'.format(index_tag_cols, config_params['sample_table_file'])

    # Check for duplicates in the index tags
    num_dup_index_tags = sum(sample_table_lane[index_tag_cols].duplicated())
    assert num_dup_index_tags == 0, '{} duplicate index tags within lane "{}" detected. Please review the sample table, found here: {}'.format(num_dup_index_tags, lane_id, config_params['sample_table_file'])

    # If this is a single read experiment, add a dummy column of
    # empty strings for the other read's empty index tag.
    if len(index_tag_cols) == 1:
        sample_table_lane['dummy_index_tag'] = ['' for i in range(sample_table_lane.shape[0])]
        if index_tag_reads == ['read_1']:
            index_tag_cols.append('dummy_index_tag')
        else:
            index_tag_cols.insert(0, 'dummy_index_tag')
    
    # Convert the table to a dictionary from index tag(s) to screen_name
    # and expt_id.
    # Note: the key is always a length-2 tuple. If there is only one
    # index tag, the entry for the other index tag is just empty strings.
    sample_table_lane['screen_name_and_expt_id'] = [tuple(x[1]) for x in sample_table_lane[['screen_name', 'expt_id']].iterrows()]
    index_tag_to_condition = sample_table_lane.set_index(index_tag_cols)['screen_name_and_expt_id'].to_dict()

    return index_tag_to_condition
    
def fastq_to_barseq(config_params, lane_id):

    lane_location_tab = get_lane_location_table(config_params)
    
    # Get the raw data folder 
    raw_folder = get_lane_folder(lane_id, lane_location_tab)

    # Set up the output folders for the temporary barseq file, the output count matrix,
    # and the read statistics
    tmp_data_path = get_lane_data_paths(config_params, lane_id)[0]
    data_path = get_lane_data_paths(config_params, lane_id)[1]
    reports_path = get_lane_reports_path(config_params, lane_id)
    cg_file.create_output_dir(tmp_data_path)
    cg_file.create_output_dir(data_path)
    cg_file.create_output_dir(reports_path)

    total_counts, common_primer_counts, barcodes_in_data, index_tags_in_data, read_type_dict, parsed_coords = read_fastq(config_params, raw_folder, data_path, lane_id)

    return total_counts, common_primer_counts, barcodes_in_data, index_tags_in_data, read_type_dict, parsed_coords


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
    base_case = [-1, '', -1, -1, -1, -1, -1]
    required_read_1_present = all([x is not base_case[i] for i,x in enumerate(read_1_params) if i in range(5)])
    required_read_2_present = all([x is not base_case[i] for i,x in enumerate(read_2_params) if i in range(5)])
    barcode_read_1_present = all([x is not base_case[i] for i,x in enumerate(read_1_params) if i in [5, 6]])
    barcode_read_2_present = all([x is not base_case[i] for i,x in enumerate(read_2_params) if i in [5, 6]])
    #print required_read_1_present, barcode_read_1_present, required_read_2_present, barcode_read_2_present
    if required_read_1_present and barcode_read_1_present and not required_read_2_present and not barcode_read_2_present:
        return {'type': 'single', 'barcode': ['read_1']}
    elif not required_read_1_present and not barcode_read_1_present and required_read_2_present and barcode_read_2_present:
        return {'type': 'single', 'barcode': ['read_2']}
    elif required_read_1_present and barcode_read_1_present and required_read_2_present and not barcode_read_2_present:
        return {'type': 'paired', 'barcode': ['read_1']}
    elif required_read_1_present and not barcode_read_1_present and required_read_2_present and barcode_read_2_present:
        return {'type': 'paired', 'barcode': ['read_2']}
    elif required_read_1_present and barcode_read_1_present and required_read_2_present and barcode_read_2_present:
        return {'type': 'paired', 'barcode': ['read_1', 'read_2']}
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


def read_fastq(config_params, folder, out_path, lane_id):

    amplicon_struct_params = get_amplicon_struct_params(config_params)

    # Get all possible common primer/index tag/barcode parameters, then determine
    # how to proceed.
    read_1_params = get_seq_params(amplicon_struct_params, 'read_1')
    read_2_params = get_seq_params(amplicon_struct_params, 'read_2')

    read_type_dict = determine_read_type(read_1_params, read_2_params)

    common_primer_tolerance = int(config_params['common_primer_tolerance'])
    
    out_filename = get_barseq_filename(config_params, lane_id)
    # print out_filename
    # print common_primer_tolerance

    of = open(out_filename, 'wt')

    # Need a function to return pairs of filenames, for which a None object
    # is one of the filenames if this is a single read run
    #fastq_filenames = [os.path.join(folder, x) for x in os.listdir(folder) if is_fastq_filename(x)]
    fastq_filenames = get_fastq_filename_list(folder, read_type_dict['type'], read_type_dict['barcode'])

    common_primer_count = 0
    index_tags = set()
    barcodes = set()
    for filename_pair in fastq_filenames:
        handles = [cfo.get_compressed_file_handle(x) if x is not None else it.cycle('\n') for x in filename_pair]
        line_count = 0
        while True:
        #for line_count, line in enumerate(f):
            try:
                line_1 = handles[0].next()
                line_2 = handles[1].next()
            except StopIteration:
                break
            if get_verbosity(config_params) >= 3:
                print line_count, line_1, line_2
            if line_count % 4 == 1:
                #string = line.strip()
                common_primer_1 = line_1[read_1_params[0]:read_1_params[2]]
                common_primer_dist_1 = jf.hamming_distance(unicode(common_primer_1), unicode(read_1_params[1]))
                common_primer_2 = line_2[read_2_params[0]:read_2_params[2]]
                common_primer_dist_2 = jf.hamming_distance(unicode(common_primer_2), unicode(read_2_params[1]))
                if get_verbosity(config_params) >= 3:
                    print common_primer_1, read_1_params[1], common_primer_dist_1
                    print common_primer_2, read_2_params[1], common_primer_dist_2
                if common_primer_dist_1 <= common_primer_tolerance and common_primer_dist_2 <= common_primer_tolerance:
                    common_primer_count += 1
                    index_tag_1 = line_1[read_1_params[3]:read_1_params[4]]
                    barcode_1 = line_1[read_1_params[5]:read_1_params[6]]
                    index_tag_2 = line_2[read_2_params[3]:read_2_params[4]]
                    barcode_2 = line_2[read_2_params[5]:read_2_params[6]]
                    index_tags.add((index_tag_1, index_tag_2))
                    barcodes.add((barcode_1, barcode_2))
                    # print "index_tag, barcode : {}, {}".format(index_tag, barcode)
                    of.write('{0}\t{1}\t{2}\t{3}\n'.format(index_tag_1, barcode_1, index_tag_2, barcode_2))
            line_count += 1
                
        for h in handles:
            if hasattr(h, 'read'):
                h.close()

    total_counts = (line_count + 1)/ 4

    of.close()

    parsed_lengths = [read_1_params[4] - read_1_params[3], read_1_params[6] - read_1_params[5],
            read_2_params[4] - read_2_params[3], read_2_params[6] - read_2_params[5]]

    shifts = range(4)
    parsed_coords = []
    for i, length in enumerate(parsed_lengths):
        if i == 0:
            parsed_coords.append(0)
        else:
            parsed_coords.append(parsed_coords[-1] + 1)
        parsed_coords.append(parsed_coords[-1] + parsed_lengths[i])

    return total_counts, common_primer_count, barcodes, index_tags, read_type_dict, parsed_coords

def correct_barcode_map(config_params, barcodes_in_data, barcode_to_gene):

    barcode_tolerance = int(config_params['barcode_tolerance'])

    barcode_correcting_maps = [None for i in (0, 1)]
    for i in (0, 1):
        orig_barcodes = set([x[i] for x in barcode_to_gene.keys()])
        observed_barcodes = set([x[i] for x in barcodes_in_data])
        full_map = {x:x for x in orig_barcodes}
        correcting_map = {}

        unmatched_barcodes = observed_barcodes.difference(orig_barcodes)
        for orig_barcode in orig_barcodes:
            for unmatched_barcode in unmatched_barcodes:
                barcode_dist = jf.hamming_distance(unicode(orig_barcode), unicode(unmatched_barcode))
                if barcode_dist <= barcode_tolerance:
                    if get_verbosity(config_params) >= 3:
                        print 'bad : corrected --> {0} : {1}'.format(unmatched_barcode, orig_barcode)
                    if correcting_map.has_key(unmatched_barcode):
                        correcting_map[unmatched_barcode].append(orig_barcode)
                    else:
                        correcting_map[unmatched_barcode] = [orig_barcode]

        # Now, filter out any unmatched barcodes that map to multiple original barcodes
        for key in correcting_map.keys():
            if len(correcting_map[key]) > 1:
                correcting_map[key].pop()

        # The corrected barcodes are still lists - turn them back into strings!
        corrected_barcodes = correcting_map.keys()
        for barcode in corrected_barcodes:
            correcting_map[barcode] = correcting_map[barcode][0]

        # Update the mapping of original barcodes to themselves with the mapping of
        # unmatched barcodes to original barcodes
        full_map.update(correcting_map)
        barcode_correcting_maps[i] = full_map
    
    return barcode_correcting_maps


def get_barseq_matrix(config_params, lane_id, parsed_coords, barcode_to_gene, barcode_correcting_map, index_tag_to_condition):

    # Iterates over barseq text file to fill the matrix

    # Set up the matrix
    barcodes = barcode_to_gene.keys()
    bc_ind = {j:i for i,j in enumerate(barcodes)}

    index_tags = index_tag_to_condition.keys()
    tag_ind = {j:i for i,j in enumerate(index_tags)}
    
    matrix = np.zeros([len(barcodes), len(index_tags)], dtype = np.int)

    # Open the barseq file
    barseq_filename = get_barseq_filename(config_params, lane_id)
    f = open(barseq_filename, 'rt')

    for line in f:
        index_tag = (line[parsed_coords[0]:parsed_coords[1]], line[parsed_coords[4]:parsed_coords[5]])
        barcode = (line[parsed_coords[2]:parsed_coords[3]], line[parsed_coords[6]:parsed_coords[7]])
        #print index_tag
        #print barcode
        
        try:
            # Since each barcode is now a 2-tuple, must correct each one
            # individually (might be a bit slower, but much less complex
            # than building a dictionary of every possible 2-tuple variation.
            barcode_fixed = (barcode_correcting_map[0][barcode[0]], barcode_correcting_map[1][barcode[1]])
            #print barcode_fixed
            row = bc_ind[barcode_fixed]
            col = tag_ind[index_tag]
        except KeyError:
            continue
        
        matrix[row, col] += 1

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
    sufficient_count_barcodes = np.mean(matrix, axis = 1) > mean_thresh
    sufficient_count_tags = np.mean(matrix, axis = 0) > mean_thresh

    if get_verbosity(config_params) >= 1:
        print "number of barcodes with mean counts per condition < {}: {}".format(mean_thresh, sum(np.invert(sufficient_count_barcodes)))
        print "number of index tags with mean counts per strain < {}: {}".format(mean_thresh, sum(np.invert(sufficient_count_tags)))
    
    filtered_matrix = matrix[np.ix_(sufficient_count_barcodes, sufficient_count_tags)]
    filtered_barcodes = [x for i,x in enumerate(barcodes) if sufficient_count_barcodes[i]]
    filtered_index_tags = [x for i,x in enumerate(index_tags) if sufficient_count_tags[i]]

    assert len(barcodes) == matrix.shape[0], "number of barcodes does not match number of rows in matrix after sufficient-count filtering."
    assert len(index_tags) == matrix.shape[1], "number of index tags does not match number of columns in matrix after sufficient-count filtering."

    f.close()

    return filtered_barcodes, filtered_index_tags, filtered_matrix, matrix

def get_sorted_counts(label, counts):

    label = np.array(label)
    # Get in descending order - the [::-1] reverses the array
    sort_indices = counts.argsort()[::-1]
    labels_sorted = label[sort_indices]
    counts_sorted = counts[sort_indices]

    return labels_sorted, counts_sorted

def write_distribution(labels, counts, output_folder, lane_id, file_string, xlabel):

    plot_filename = os.path.join(output_folder, '{0}_{1}_distribution.png'.format(lane_id, file_string))
    x = np.array(range(len(counts)))

    plt.figure()
    plt.fill_between(x, y1 = 0, y2 = counts, step = 'post')
    plt.ylabel('Number of occurrences')
    plt.xlabel(xlabel)
    plt.savefig(plot_filename)

    text_filename = os.path.join(output_folder, '{0}_{1}_distribution.txt'.format(lane_id, file_string))
    of = open(text_filename, 'wt')
    of.write('{0}\tNumber of occurrences\n'.format(file_string))
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

def dump_count_matrix(config_params, lane_id, barcodes, conditions, matrix):

    out_path = get_lane_data_paths(config_params, lane_id)[1]
    out_filename = os.path.join(out_path, '{0}_{1}'.format(lane_id, 'barseq_matrix.dump.gz'))

    of = gzip.open(out_filename, 'wb')

    dataset = [barcodes, conditions, matrix]
    cPickle.dump(dataset, of)

    of.close()


def main(config_file, lane_id):

    config_params = parse_yaml(config_file)
    amplicon_struct_params = get_amplicon_struct_params(config_params)

    # Get maps of barcode to barcode_gene (keeps the strains unique/traceable), and index tag to condition
    if get_verbosity(config_params) >= 1:
        print 'creating mappings from barcodes and index tags...'
    barcode_to_gene = get_barcode_to_gene(config_params)
    index_tag_to_condition = get_index_tag_to_condition(config_params, lane_id)

    # Loop over the raw fastq files, write out the "index_tag\tbarcode" file,
    # and return all encountered index tags and barcodes
    if get_verbosity(config_params) >= 1:
        print 'parsing fastq file(s)...'
    total_counts, common_primer_counts, barcodes_in_data, index_tags_in_data, read_type_dict, parsed_coords = fastq_to_barseq(config_params, lane_id)

    if get_verbosity(config_params) >= 1:
        print 'barcode to gene map: {}'.format(barcode_to_gene.items()[0:5])
        print 'index tag to condition map: {}'.format(index_tag_to_condition.items()[0:5])

    # Correct the barcodes within the specified error tolerance
    # (There is no function to correct the index tags - this could easily be
    # written in later, although we see no need for it)
    if get_verbosity(config_params) >= 1:
        print 'correcting barcodes...'
    barcode_correcting_map = correct_barcode_map(config_params, barcodes_in_data, barcode_to_gene)
    if get_verbosity(config_params) >= 1:
        print 'number of read_1 barcodes that will be counted: {}'.format(len(barcode_correcting_map[0]))
        print 'number of read_2 barcodes that will be counted: {}'.format(len(barcode_correcting_map[1]))

    # Loop over the barseq file (index_tag_1\tbarcode_1\tindex_tag_2\tbarcode_2\n)
    # and assemble the matrix of read counts
    if get_verbosity(config_params) >= 1:
        print 'generating barseq matrix...'
    corrected_barcodes, index_tags, matrix, unfiltered_matrix = get_barseq_matrix(config_params, lane_id, parsed_coords, barcode_to_gene, barcode_correcting_map, index_tag_to_condition)

    if get_verbosity(config_params) >= 1:
        print 'number of barcodes: {}'.format(len(corrected_barcodes))
        print 'number of index tags: {}'.format(len(index_tags))
        print 'matrix shape: {0} rows x {1} columns'.format(*matrix.shape)

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
