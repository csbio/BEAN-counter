#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# Function to parse the config file (in the format "paramter: value")
# Lines beginning with '#' are ignored, as are blank lines

def get_key_val(line):
    
    split_line = line.lstrip().split(':')
    
    if len(split_line) > 2:
        raise Exception('One ":" allowed per line in configuration files')

    key = split_line[0].rstrip()
    val = split_line[1].lstrip().rstrip()

    return key, val



def parse(filename):

    f = open(filename)

    params = {}

    for line in f:
        # print line
        if line.lstrip().startswith('#'):
            continue
        if line.lstrip() == '':
            continue

        param, value = get_key_val(line)
        
        if params.has_key(param):
            raise Exception('Duplicate values for one parameter!')

        params[param] = value

    return params


def parse_species_config(filename):

    f = open(filename)

    spec_params = {}

    while True:
        try:
            line = f.next()
        except:
            break

        # print line
        if line.lstrip().startswith('#'):
            continue
        if line.lstrip() == '':
            continue
        
        if line.startswith('species_ID'):
            # print 'Species_start'
            spec_param, spec_id = get_key_val(line)
            # print spec_param, spec_id

            spec_params[spec_id] = {}
           
            while True:
                line = f.next()
                if line.lstrip().startswith('---'):
                    break

                param, value = get_key_val(line)
                if spec_params[spec_id].has_key(param):
                    raise Exception('Duplicate values for one parameter!')
                spec_params[spec_id][param] = value

    return spec_params
