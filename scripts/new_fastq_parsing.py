

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

def initialize_dicts_arrays(read_type_dict, amplicon_struct_params, sample_tab, barcode_tab):
    '''
    read_inds, titles, match_dicts, array_ind_dicts, array = initialize_dicts_arrays(...)
    '''
    if read_type_dict['type'] == 'single':
        index_tag_col = amplicon_struct_params[read_type_dict['barcode'][0]]['index_tag']['sample_table_column']
        assert index_tag_col in sample_tab.columns, 'sample_table_column "{}" specified in amplicon_struct_file is not present in the sample table.'.format(index_tag_col)
        index_tags = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col])}
        index_tags['multi_match'] = len(index_tags)
        index_tags['no_match'] = len(index_tags)
        barcode_col = amplicon_struct_params[read_type_dict['barcode'][0]]['genetic_barcode']['barcode_file_column']
        assert barcode_col in barcode_tab.columns, 'barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col)
        barcodes = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col])}
        barcodes['multi_match'] = len(barcodes)
        barcodes['no_match'] = len(barcodes)
        count_array = np.zeros((len(index_tags), len(barcodes)), dtype = np.int)
        return [0, 0], ['index_tag', 'barcode'], [{}, {}], [index_tags, barcodes], count_array
    elif read_type_dict['type'] == 'paired':
        # Since there are always 2 index tags in a paired read scheme, I can write this code once.
        index_tag_col_1 = amplicon_struct_params['read_1']['index_tag']['sample_table_column']
        assert index_tag_col_1 in sample_tab.columns, 'read_1 sample_table_column "{}" specified in amplicon_struct_file is not present in the sample table.'.format(index_tag_col)
        index_tag_col_2 = amplicon_struct_params['read_2']['index_tag']['sample_table_column']
        assert index_tag_col_2 in sample_tab.columns, 'read_2 sample_table_column "{}" specified in amplicon_struct_file is not present in the sample table.'.format(index_tag_col)
        if read_type_dict['barcode'] == ['read_1']:
            index_tags_1 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_1])}
            index_tags_1['multi_match'] = len(index_tags_1)
            index_tags_1['no_match'] = len(index_tags_1)
            barcode_col_1 = amplicon_struct_params['read_1']['genetic_barcode']['barcode_file_column']
            assert barcode_col_1 in barcode_tab.columns, 'read_1 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col_1)
            barcodes_1 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_1])}
            barcodes_1['multi_match'] = len(barcodes_1)
            barcodes_1['no_match'] = len(barcodes_1)
            index_tags_2 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_2])}
            index_tags_2['multi_match'] = len(index_tags_2)
            index_tags_2['no_match'] = len(index_tags_2)
            count_array = np.zeros((len(index_tags_1), len(barcodes_1), len(index_tags_2)), dtype = np.int)
            return [0, 0, 1], ['index_tag_1', 'barcode_1', 'index_tag_2'], [{}, {}, {}], [index_tags_1, barcodes_1, index_tags_2], count_array
        elif read_type_dict['barcode'] == ['read_2']:
            index_tags_1 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_1])}
            index_tags_1['multi_match'] = len(index_tags_1)
            index_tags_1['no_match'] = len(index_tags_1)
            index_tags_2 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_2])}
            index_tags_2['multi_match'] = len(index_tags_2)
            index_tags_2['no_match'] = len(index_tags_2)
            barcode_col_2 = amplicon_struct_params['read_2']['genetic_barcode']['barcode_file_column']
            assert barcode_col_2 in barcode_tab.columns, 'read_2 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col_2)
            barcodes_2 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_2])}
            barcodes_2['multi_match'] = len(barcodes_2)
            barcodes_2['no_match'] = len(barcodes_2)
            count_array = np.zeros((len(index_tags_1), len(index_tags_2), len(barcodes_2)), dtype = np.int)
            return [0, 1, 1], ['index_tag_1', 'index_tag_2', 'barcode_2'], [{}, {}, {}], [index_tags_1, index_tags_2, barcodes_2], count_array
        elif read_type_dict['barcode'] == ['read_1', 'read_2']:
            index_tags_1 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_1])}
            index_tags_1['multi_match'] = len(index_tags_1)
            index_tags_1['no_match'] = len(index_tags_1)
            barcode_col_1 = amplicon_struct_params['read_1']['genetic_barcode']['barcode_file_column']
            assert barcode_col_1 in barcode_tab.columns, 'read_1 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene_barcode table.'.format(barcode_col_1)
            barcodes_1 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_1])}
            barcodes_1['multi_match'] = len(barcodes_1)
            barcodes_1['no_match'] = len(barcodes_1)
            index_tags_2 = {x:i, for i, x in enumerate(sample_tab.loc[:, index_tag_col_2])}
            index_tags_2['multi_match'] = len(index_tags_2)
            index_tags_2['no_match'] = len(index_tags_2)
            barcode_col_2 = amplicon_struct_params['read_2']['genetic_barcode']['barcode_file_column']
            assert barcode_col_2 in barcode_tab.columns, 'read_2 barcode_file_column "{}" specified in amplicon_struct_file is not present in the gene barcode table.'.format(barcode_col_2)
            barcodes_2 = {x:i for i, x in enumerate(barcode_tab.loc[:, barcode_col_2])}
            barcodes_2['multi_match'] = len(barcodes_2)
            barcodes_2['no_match'] = len(barcodes_2)
            count_array = np.zeros((len(index_tags_1), len(barcodes_1), len(index_tags_2), len(barcodes_2)), dtype = np.int)
            return [0, 0, 1, 1], ['index_tag_1', 'barcode_1', 'index_tag_2', 'barcode_2'], [{}, {}, {}, {}], [index_tags_1, barcodes_1, index_tags_2, barcodes_2], count_array
    

def parse_seqs():


    amplicon_struct_params = get_amplicon_struct_params(config_params)

    # Get all possible common primer/index tag/barcode parameters, then determine
    # how to proceed.
    read_1_params = get_seq_params(amplicon_struct_params, 'read_1')
    read_2_params = get_seq_params(amplicon_struct_params, 'read_2')

    read_type_dict = determine_read_type(read_1_params, read_2_params)

    read_inds, titles, match_dicts, array_ind_dicts, array = initialize_dicts_arrays(read_type_dict, amplicon_struct_params, sample_tab, barcode_tab)
