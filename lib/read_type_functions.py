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
    except (KeyError, TypeError):
        index_tag_start = -1
    try:
        barcode_start = int(amplicon_struct_params[read]['genetic_barcode']['start'])
    except (KeyError, TypeError):
        barcode_start = -1

    return {'common_primer' : { 'start' : common_primer_start , 'end' : common_primer_end, 'seq' : common_primer_seq },
            'index_tag' : { 'start' : index_tag_start },
            'barcode' : {'start' : barcode_start }
            }

def determine_read_type(read_1_params, read_2_params):
    base_case =  {'common_primer' : { 'start' : -1 , 'end' : -1, 'seq' : '' },
                  'index_tag' : { 'start' : -1 },
                  'barcode' : {'start' : -1 }
                  }
    
    required_read_1_list = []
    required_read_2_list = []
    for x in ['common_primer', 'index_tag']:
        for y in base_case[x].keys():
            required_read_1_list.append(read_1_params[x][y] is not base_case[x][y])
            required_read_2_list.append(read_2_params[x][y] is not base_case[x][y])
    
    barcode_read_1_list = []
    barcode_read_2_list = []
    for y in base_case['barcode'].keys():
        barcode_read_1_list.append(read_1_params['barcode'][y] is not base_case['barcode'][y])
        barcode_read_2_list.append(read_2_params['barcode'][y] is not base_case['barcode'][y])
    
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
