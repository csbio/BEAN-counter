
from Bio.trie import trie

# This is where I'm putting all potential new code before I attempt to
# integrate it into the old code.

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
                for line in f:
                    yield [line.rstrip()]
    else:
        for fname_pair in filenames:
            with cfo.get_compressed_file_handle(fname_pair[0]) as f1, cfo.get_compressed_file_handle(fname_pair[1]) as f2:
                for line1, line2 in it.izip(f1, f2):
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
    read_inds, titles, match_dicts, array_ind_dicts, seq_trie_list, array = initialize_dicts_arrays(...)
    '''
    if read_type_dict['type'] == 'single':
        index_tag_col = amplicon_struct_params[read_type_dict['barcode'][0]]['index_tag']['sample_table_column']
        assert index_tag_col in sample_tab.columns, 'sample_table_column "{}" specified in amplicon_struct_file is not present in the sample table.'.format(index_tag_col)
        index_tag_tries, index_tag_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col])
        index_tags = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col])}
        index_tags['multi_match'] = len(index_tags)
        index_tags['no_match'] = len(index_tags) + 1
        barcode_col = amplicon_struct_params[read_type_dict['barcode'][0]]['genetic_barcode']['barcode_file_column']
        assert barcode_col in barcode_tab.columns, 'barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col)
        barcode_tries, barcode_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col])
        barcodes = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col])}
        barcodes['multi_match'] = len(barcodes)
        barcodes['no_match'] = len(barcodes) + 1
        index_tag_tol = 0
        barcode_tol = config_params['barcode_tolerance']
        count_array = np.zeros((len(index_tags) + 2, len(barcodes) + 2), dtype = np.int)
        return [0, 0], ['index_tag', 'barcode'], [{}, {}], [index_tags, barcodes], [index_tag_tries, barcode_tries], [index_tag_lengthss, barcode_lengths], [index_tag_tol, barcode_tol], count_array
    elif read_type_dict['type'] == 'paired':
        # Since there are always 2 index tags in a paired read scheme, I can write this code once.
        index_tag_col_1 = amplicon_struct_params['read_1']['index_tag']['sample_table_column']
        assert index_tag_col_1 in sample_tab.columns, 'read_1 sample_table_column "{}" specified in amplicon_struct_file is not present in the sample table.'.format(index_tag_col)
        index_tag_col_2 = amplicon_struct_params['read_2']['index_tag']['sample_table_column']
        assert index_tag_col_2 in sample_tab.columns, 'read_2 sample_table_column "{}" specified in amplicon_struct_file is not present in the sample table.'.format(index_tag_col)
        if read_type_dict['barcode'] == ['read_1']:
            index_tag_1_tries, index_tag_1_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_1])
            index_tags_1 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_1])}
            index_tags_1['multi_match'] = len(index_tags_1)
            index_tags_1['no_match'] = len(index_tags_1) + 1
            barcode_col_1 = amplicon_struct_params['read_1']['genetic_barcode']['barcode_file_column']
            assert barcode_col_1 in barcode_tab.columns, 'read_1 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col_1)
            barcode_1_tries, barcode_1_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col_1])
            barcodes_1 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_1])}
            barcodes_1['multi_match'] = len(barcodes_1)
            barcodes_1['no_match'] = len(barcodes_1) + 1
            index_tag_2_tries, index_tag_2_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_2])
            index_tags_2 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_2])}
            index_tags_2['multi_match'] = len(index_tags_2)
            index_tags_2['no_match'] = len(index_tags_2) + 1
            index_tag_tol = 0
            barcode_tol = config_params['barcode_tolerance']
            count_array = np.zeros((len(index_tags_1) + 2, len(barcodes_1) + 2, len(index_tags_2) + 2), dtype = np.int)
            return [0, 0, 1], ['index_tag', 'barcode', 'index_tag'], [{}, {}, {}], [index_tags_1, barcodes_1, index_tags_2], [index_tag_1_tries, barcode_1_tries, index_tag_2_tries], [index_tag_1_lengths, barcode_1_lengths, index_tag_2_lengths], [index_tag_tol, barcode_tol, index_tag_tol], count_array
        elif read_type_dict['barcode'] == ['read_2']:
            index_tag_1_tries, index_tag_1_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_1])
            index_tags_1 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_1])}
            index_tags_1['multi_match'] = len(index_tags_1)
            index_tags_1['no_match'] = len(index_tags_1) + 1
            index_tag_2_tries, index_tag_2_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_2])
            index_tags_2 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_2])}
            index_tags_2['multi_match'] = len(index_tags_2)
            index_tags_2['no_match'] = len(index_tags_2) + 1
            barcode_col_2 = amplicon_struct_params['read_2']['genetic_barcode']['barcode_file_column']
            assert barcode_col_2 in barcode_tab.columns, 'read_2 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col_2)
            barcode_2_tries, barcode_2_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col_2])
            barcodes_2 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_2])}
            barcodes_2['multi_match'] = len(barcodes_2)
            barcodes_2['no_match'] = len(barcodes_2) + 1
            index_tag_tol = 0
            barcode_tol = config_params['barcode_tolerance']
            count_array = np.zeros((len(index_tags_1) + 2, len(index_tags_2) + 2, len(barcodes_2) + 2), dtype = np.int)
            return [0, 1, 1], ['index_tag', 'index_tag', 'barcode'], [{}, {}, {}], [index_tags_1, index_tags_2, barcodes_2], [index_tag_1_tries, index_tag_2_tries, barcode_2_tries], [index_tag_1_lengths, index_tag_2_lengths, barcode_2_lengths], [index_tag_tol, index_tag_tol, barcode_tol], count_array
        elif read_type_dict['barcode'] == ['read_1', 'read_2']:
            index_tag_1_tries, index_tag_1_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_1])
            index_tags_1 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_1])}
            index_tags_1['multi_match'] = len(index_tags_1)
            index_tags_1['no_match'] = len(index_tags_1) + 1
            barcode_col_1 = amplicon_struct_params['read_1']['genetic_barcode']['barcode_file_column']
            assert barcode_col_1 in barcode_tab.columns, 'read_1 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col_1)
            barcode_1_tries, barcode_1_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col_1])
            barcodes_1 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_1])}
            barcodes_1['multi_match'] = len(barcodes_1)
            barcodes_1['no_match'] = len(barcodes_1) + 1
            index_tag_2_tries, index_tag_2_lengths = gen_seq_tries(sample_tab.loc[:, index_tag_col_2])
            index_tags_2 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_2])}
            index_tags_2['multi_match'] = len(index_tags_2)
            index_tags_2['no_match'] = len(index_tags_2) + 1
            barcode_col_2 = amplicon_struct_params['read_2']['genetic_barcode']['barcode_file_column']
            assert barcode_col_2 in barcode_tab.columns, 'read_2 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene barcode table.'.format(barcode_col_2)
            barcode_2_tries, barcode_2_lengths = gen_seq_tries(barcode_tab.loc[:, barcode_col_2])
            barcodes_2 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_2])}
            barcodes_2['multi_match'] = len(barcodes_2)
            barcodes_2['no_match'] = len(barcodes_2) + 1
            index_tag_tol = 0
            barcode_tol = config_params['barcode_tolerance']
            count_array = np.zeros((len(index_tags_1) + 2, len(barcodes_1) + 2, len(index_tags_2) + 2, len(barcodes_2) + 2), dtype = np.int)
            return [0, 0, 1, 1], ['index_tag', 'barcode', 'index_tag', 'barcode'], [{}, {}, {}, {}], [index_tags_1, barcodes_1, index_tags_2, barcodes_2], [index_tag_1_tries, barcode_1_tries, index_tag_2_tries, barcode_2_tries], [index_tag_1_lengths, barcode_1_lengths, index_tag_2_lengths, barcode_2_lengths], [index_tag_tol, index_tag_tol, barcode_tol, barcode_tol], count_array
    

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

    for l in lengths:
        res = seq_trie_length_list[i].get_approximate(seq[0:l], n)
        if len(res) == 0:
            continue

        # If res meets the mismatch criterion, return before examining any others
        # But must account for multi-matches!
        if len(res) > 1:
            res = list(set(res))
            res.sort(key = lambda x: x[2])
            if res[0][2] == res[1][2]:
                return 'multi_match'
            else:
                return res[0]
        else:
            return res[0]

    # If nothing matched at all, return "no_match"!
    return ['no_match']


def parse_seqs(lane_id, config_params):

    amplicon_struct_params = get_amplicon_struct_params(config_params)

    # Get all possible common primer/index tag/barcode parameters, then determine
    # how to proceed.
    read_params = [get_seq_params(amplicon_struct_params, 'read_1'), get_seq_params(amplicon_struct_params, 'read_2')]

    read_type_dict = determine_read_type(read_1_params, read_2_params)

    read_inds, seq_types, match_dicts, array_ind_dicts, seq_trie_lists, seq_lengths, tols, array = initialize_dicts_arrays(read_type_dict, amplicon_struct_params, config_params, sample_tab, barcode_tab)

    lane_location_tab = get_lane_location_table(config_params))
    folder = get_lane_folder(lane_id, lane_location_tab)

    n = len(read_inds)
    idxs = [None] * len(read_inds)
    for line_list in line_gen(folder, read_type_dict['type']):
        # Check for common primer...
        
        # Match on index tags/barcodes
        for i in n:
            start_coord = read_params[read_inds[i]][seq_types[i]]['start']
            end_coord = read_params[read_inds[i]][seq_types[i]]['end']
            seq = line_list[start_coord:end_coord]

            try:
                corrected_seq = match_dicts[i][seq][0]
                match_dicts[i][seq][1] += 1
            except KeyError as e:
                # need to write match_seq function
                # Returns a list of [corrected_seq, 1] where the "1" represents the first occurrence
                match_dicts[i][seq] = match_seq(seq, seq_trie_lists[i], tols[i], seq_lengths)
                corrected_seq = match_dicts[i][seq][0]

            idxs[i] = array_ind_dicts[i][corrected_seq]

        array[tuple(idxs)] += 1
