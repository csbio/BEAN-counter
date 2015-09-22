# This script takes in the barseq_counter main configuration and species configuration files
# and a sequencing lane identifier. It exports 1) reports on the total number
# of reads, the number of reads matching the common primer, and the number of reads that
# match index tags and barcodes; 2) a file of "index_tag\tbarcode" for all sequences matching
# the common primer; and 3) a matrix of counts, rows being the barcodes (genes) and columns
# being the conditions.

import pandas as pd
import numpy as np
import sys, os, gzip
import jellyfish as jf
import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
import matplotlib.pyplot as plt
import cPickle

barseq_path = os.getenv('BARSEQ_PATH')
sys.path.append(os.path.join(barseq_path, 'lib'))

import config_file_parser as cfp
import compressed_file_opener as cfo
import cg_file_tools as cg_file

def get_species_config_params(config_params):
    species_config_file = './data/species_config_file.txt'
    all_species_params = cfp.parse_species_config(species_config_file)
    species_id = config_params['species_ID']
    species_params = all_species_params[species_id]
    return species_params

def get_sample_table(config_params):

    filename = config_params['sample_table_file']
    
    # Read everything in as a string, to prevent vexing
    # number interpretation problems! Methods further down
    # can coerce to different types.
    tab = pd.read_table(filename, dtype = 'S')
    return tab

# def get_barcode_table(config_params):
# 
#     filename = config_params['barcode_table_file']
#     tab = pd.read_table(filename)
#     return tab

def get_lane_location_table(config_params):

    filename = config_params['lane_location_file']
    tab = pd.read_table(filename, dtype = 'S')
    return tab

def get_lane_folder(lane_id, lane_location_tab):

    lane_location_tab = lane_location_tab.set_index('lane_id')
    lane_folder = lane_location_tab.loc[lane_id, 'location']
    return lane_folder

def get_lane_data_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'intermediate', lane_id)

def get_lane_reports_path(config_params, lane_id):

    output_folder = config_params['output_folder']
    return os.path.join(output_folder, 'reports', lane_id)

def get_barseq_filename(config_params, lane_id):

    path = get_lane_data_path(config_params, lane_id)
    return os.path.join(path, '{0}_{1}'.format(lane_id, 'barseq.txt.gz'))

def get_barcode_to_gene(species_config_params):

    filename = species_config_params['gene_barcode_file']
    full_path = os.path.join('./data/barcodes/{0}'.format(filename))

    f = open(full_path)
    barcode_2_gene = {}
    for line in f:
        barcode, gene = line.rstrip().split('\t')
        if barcode_2_gene.has_key(barcode):
            sys.exit('duplicate barcodes detected')
        else:
            barcode_2_gene[barcode] = '{0}_{1}'.format(gene, barcode)

    f.close()
    return barcode_2_gene

def get_index_tag_to_condition(config_params, lane_id):

    sample_table = get_sample_table(config_params)
    sample_table = sample_table.set_index('lane')

    sample_table_lane = sample_table.loc[lane_id]

    index_tags = list(sample_table_lane['index_tag'])
    screen_names = list(sample_table_lane['screen_name'])
    expt_ids = list(sample_table_lane['expt_id'])

    index_tag_to_condition = {}

    for i in range(len(index_tags)):
        index_tag = index_tags[i]
        screen_name = screen_names[i]
        expt_id = expt_ids[i]
        
        index_tag_to_condition[index_tag] = '{0}-{1}'.format(screen_name, expt_id)

    return index_tag_to_condition

def fastq_to_barseq(config_params, species_config_params, lane_id):

    lane_location_tab = get_lane_location_table(config_params)
    
    # Get the raw data folder 
    raw_folder = get_lane_folder(lane_id, lane_location_tab)

    # Set up the output folders for data and read stats
    data_path = get_lane_data_path(config_params, lane_id)
    reports_path = get_lane_reports_path(config_params, lane_id)
    cg_file.create_output_dir(data_path)
    cg_file.create_output_dir(reports_path)

    total_counts, common_primer_counts, barcodes_in_data, index_tags_in_data = read_fastq(config_params, species_config_params, raw_folder, data_path, lane_id)

    return total_counts, common_primer_counts, barcodes_in_data, index_tags_in_data


def read_fastq(config_params, species_config_params, folder, out_path, lane_id):

    common_primer_start = int(species_config_params['common_primer_start'])
    common_primer_end = common_primer_start + int(species_config_params['common_primer_length'])
    common_primer_seq = species_config_params['common_primer_sequence']
    index_tag_start = int(species_config_params['index_tag_start'])
    index_tag_end = index_tag_start + int(species_config_params['index_tag_length'])
    barcode_start = int(species_config_params['genetic_barcode_start'])
    barcode_end = barcode_start + int(species_config_params['genetic_barcode_length'])

    out_filename = get_barseq_filename(config_params, lane_id)

    of = gzip.open(out_filename, 'wb')

    fastq_filenames = [os.path.join(folder, x) for x in os.listdir(folder) if x.count('fastq') > 0]

    line_count = 0
    common_primer_count = 0
    barcodes = set()
    index_tags = set()
    for filename in fastq_filenames:
        compressed_file_opener = cfo.get_compressed_file(filename)
        f = compressed_file_opener.open()
        for line in f:
            line_count += 1
            if line_count % 4 == 2:
                string = line.strip()
                common_primer = string[common_primer_start:common_primer_end]
                common_primer_dist = jf.hamming_distance(common_primer, common_primer_seq)
                if common_primer_dist <= 2:
                    common_primer_count += 1
                    index_tag = string[index_tag_start:index_tag_end]
                    barcode = string[barcode_start:barcode_end]
                    index_tags.update(set([index_tag]))
                    barcodes.update(set([barcode]))
                    of.write('{0}\t{1}\n'.format(index_tag, barcode))
        f.close()

    total_counts = line_count / 4

    of.close()

    return total_counts, common_primer_count, barcodes, index_tags

def correct_barcode_map(config_params, barcodes_in_data, barcode_to_gene):

    barcode_error_cutoff = int(config_params['barcode_error_cutoff'])

    orig_barcodes = set(barcode_to_gene.keys())
    full_map = {x:x for x in orig_barcodes}
    correcting_map = {}

    unmatched_barcodes = set(barcodes_in_data).difference(orig_barcodes)
    for orig_barcode in barcode_to_gene.keys():
        for unmatched_barcode in unmatched_barcodes:
            barcode_dist = jf.hamming_distance(orig_barcode, unmatched_barcode)
            if barcode_dist <= barcode_error_cutoff:
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
    
    return full_map

def get_barseq_matrix(config_params, lane_id, barcode_to_gene, barcode_correcting_map, index_tag_to_condition):

    # Iterates over barseq text file to fill the matrix

    # Set up the matrix
    barcodes = barcode_to_gene.keys()
    bc_ind = {j:i for i,j in enumerate(barcodes)}

    index_tags = index_tag_to_condition.keys()
    tag_ind = {j:i for i,j in enumerate(index_tags)}
    
    matrix = np.zeros([len(barcodes), len(index_tags)], dtype = np.int)

    # Open the barseq file
    barseq_filename = get_barseq_filename(config_params, lane_id)
    f = gzip.open(barseq_filename)

    for line in f:
        index_tag, barcode = line.rstrip().split('\t')
        
        try:
            barcode_fixed = barcode_correcting_map[barcode]
            row = bc_ind[barcode_fixed]
            col = tag_ind[index_tag]
        except KeyError:
            continue
        
        matrix[row, col] += 1

    f.close()

    return barcodes, index_tags, matrix

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
    plt.bar(x, counts)
    plt.ylabel('Number of occurrences')
    plt.xlabel(xlabel)
    plt.savefig(plot_filename)

    text_filename = os.path.join(output_folder, '{0}_{1}_distribution.txt'.format(lane_id, file_string))
    of = open(text_filename, 'wt')
    of.write('{0}\tNumber of occurrences\n'.format(file_string))
    for i in range(len(labels)):
        of.write('{0}\t{1}\n'.format(labels[i], counts[i]))
    of.close()

def write_summary(total_counts, common_primer_counts, total_index_barcode_counts, matrix_shape, output_folder, lane_id):

    filename = os.path.join(output_folder, '{0}_summary.txt'.format(lane_id))

    common_primer_vs_total_percent = float(common_primer_counts) / float(total_counts) * 100
    tag_barcode_match_vs_total_percent = float(total_index_barcode_counts) / float(total_counts) * 100
    tag_barcode_match_vs_common_primer_percent = float(total_index_barcode_counts) / float(common_primer_counts) * 100

    of = open(filename, 'wt')
    of.write('Total number of reads: {0}\n'.format(total_counts))
    of.write('Number of reads with common primer: {0} ({1:.2f} %)\n'.format(common_primer_counts, common_primer_vs_total_percent))
    of.write('Number of reads that match index tags and genetic barcodes: {0} ({1:.2f} % of total counts, {2:.2f} % of common primer counts)\n'.format(total_index_barcode_counts, tag_barcode_match_vs_total_percent, tag_barcode_match_vs_common_primer_percent))

    # Write matrix size
    of.write('Size of the data: {0} barcodes (gene mutants) x {1} index tags (conditions)'.format(*matrix_shape))

    of.close()

def generate_reports(config_params, lane_id, gen_barcodes, index_tags, matrix, total_counts, common_primer_counts):

    output_folder = get_lane_reports_path(config_params, lane_id)

    gen_barcode_counts = matrix.sum(axis = 1)
    index_tag_counts = matrix.sum(axis = 0)
    total_index_barcode_counts = matrix.sum()

    index_tags_sorted, index_tag_counts_sorted = get_sorted_counts(index_tags, index_tag_counts)
    gen_barcodes_sorted, gen_barcode_counts_sorted = get_sorted_counts(gen_barcodes, gen_barcode_counts)

    write_distribution(index_tags_sorted, index_tag_counts_sorted, output_folder, lane_id, 'index_tag', 'Index tags (sorted)')
    write_distribution(gen_barcodes_sorted, gen_barcode_counts_sorted, output_folder, lane_id, 'barcode', 'Genetic barcodes (sorted)')

    write_summary(total_counts, common_primer_counts, total_index_barcode_counts, matrix.shape, output_folder, lane_id)

def dump_count_matrix(config_params, lane_id, barcodes, conditions, matrix):

    out_path = get_lane_data_path(config_params, lane_id)
    out_filename = os.path.join(out_path, '{0}_{1}'.format(lane_id, 'barseq_matrix.dump.gz'))

    of = gzip.open(out_filename, 'wb')

    dataset = [barcodes, conditions, matrix]
    cPickle.dump(dataset, of)

    of.close()


def main(config_file, lane_id):

    print 'parsing parameters...'
    config_params = cfp.parse(config_file)
    species_config_params = get_species_config_params(config_params)

    # Loop over the raw fastq files, write out the "index_tag\tbarcode" file,
    # and return all encountered index tags and barcodes
    print 'parsing fastq file(s)...'
    total_counts, common_primer_counts, barcodes_in_data, index_tags_in_data = fastq_to_barseq(config_params, species_config_params, lane_id)

    # Get maps of barcode to barcode_gene (keeps the strains unique/traceable), and index tag to condition
    print 'creating mappings from barcodes and index tags...'
    barcode_to_gene = get_barcode_to_gene(species_config_params)
    index_tag_to_condition = get_index_tag_to_condition(config_params, lane_id)

    print 'barcode to gene map: {}'.format(barcode_to_gene.items()[0:5])
    print 'index tag to condition map: {}'.format(index_tag_to_condition.items()[0:5])

    # Correct the barcodes within the specified error tolerance
    # (There is no function to correct the index tags - this could easily be
    # written in later, although we see no need for it)
    print 'correcting barcodes...'
    barcode_correcting_map = correct_barcode_map(config_params, barcodes_in_data, barcode_to_gene)
    print 'number of barcodes that will be counted: {}'.format(np.array([len(x) for x in barcode_correcting_map]).sum())

    # Loop over the barseq file (index_tag\tbarcode) and assemble the matrix of read counts
    print 'generating barseq matrix...'
    corrected_barcodes, index_tags, matrix = get_barseq_matrix(config_params, lane_id, barcode_to_gene, barcode_correcting_map, index_tag_to_condition)

    print 'number of barcodes: {}'.format(len(corrected_barcodes))
    print 'number of index tags: {}'.format(len(index_tags))
    print 'matrix shape: {0} rows x {1} columns'.format(*matrix.shape)

    # Generate reports for index tags and barcodes
    print 'generating reports...'
    generate_reports(config_params, lane_id, corrected_barcodes, index_tags, matrix, total_counts, common_primer_counts)

    # Convert the barcodes to their condition names and unique gene/barcode names
    barcode_gene_ids = np.array([barcode_to_gene[bc] for bc in corrected_barcodes])
    condition_ids = np.array([index_tag_to_condition[tag] for tag in index_tags])

    # Dump out the final count matrix to file - other scripts will read it and turn it into a readable matrix/CDT
    print 'dumping count matrix...'
    dump_count_matrix(config_params, lane_id, barcode_gene_ids, condition_ids, matrix)
 

# call: python fastq_to_count_matrix.py <config_file> <lane_id>
if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: python fastq_to_count_matrix.py <config_file> <lane_id>'
    else:
        config_file = sys.argv[1]
        lane_id = sys.argv[2]
        main(config_file, lane_id)