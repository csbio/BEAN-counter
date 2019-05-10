#!/usr/bin/env python

VERSION='2.6.1'

import yaml
import argparse
import os

# This script must be run from the root of the BEAN-counter source code
# directory. It provides options to bump the major, minor, and patch versions
# of the software within all relevant scripts as well as roll back to the
# latest version in the history.

def parse_version_file(fname):
    with open(fname, 'rt') as f:
        version_dict = yaml.load(f)

    assert 'current' in version_dict, '\nVERSION.yaml must have a key named "current".'
    assert is_valid_version(version_dict['current']), '\nThe version string "{}" is not a valid version (must be "X.Y.Z").'.format(version_dict['current'])

    history = version_dict.get('history')
    if history is None or history == '':
        #print "detected None history"
        history = []
        #print history
    history.sort()
    history = history[::-1]
    invalid_history_strings = [x for x in history if not is_valid_version(x)]
    assert not any(invalid_history_strings), '\nThe current versions in the histroy of VERSION.yaml are invalid:\n' \
            '{}\nPlease correct before moving on.'.format('\n'.join(invalid_history_strings))
    return version_dict['current'], history

def write_version_file(fname, version, history):
    with open(fname, 'wt') as f:
        f.write(yaml.dump({'current': version}, default_flow_style = False))
        f.write(yaml.dump({'history': history}, default_flow_style = False))

def parse_version(x):
    '''
    Assumes x has been stripped of whitespace, if read in from a file for
    example.
    '''
    assert is_valid_version(x), '\nVersion string "{}" is not a valid version string.\n' \
            'Please fix before proceeding'.format(x)
    split_x = [int(xx) for xx in x.split('.')]
    return split_x

def increment_version(vers, Type):
    major, minor, patch = parse_version(vers)
    if Type == 'patch':
        patch += 1
    elif Type == 'minor':
        patch = 0
        minor += 1
    else:
        patch = 0
        minor = 0
        major += 1

    return '.'.join([str(major), str(minor), str(patch)])

def is_valid_version(x):
    split_x = x.split('.')
    if len(split_x) != 3:
        return False
    try:
        for xx in split_x:
            int(xx)
    except ValueError as e:
        return False

    return True


def rollback_version(vers, history):

    if len(history) > 0:
        return history[0], history[1:]
    else:
        return vers, history

def get_py_files(folder):
    return [os.path.join(folder, x) for x in os.listdir(folder) if x.endswith('.py')]

def check_version(fname):
    version = None
    with open(fname, 'rt') as f:
        for line in f:
            if line.startswith('VERSION'):
                if version is not None:
                    print line
                    new_version = line.rstrip().split('=')[1].strip().replace("'", "").replace('"', '')
                    print new_version
                    if version != new_version:
                        version = 'multiple'
                else:
                    version = line.rstrip().split('=')[1].strip().replace("'", "").replace('"', '')
    return version

def replace_version(fname, old, new):
    f = open(fname, 'rt')
    lines = f.readlines()
    for i in range(len(lines)):
        if lines[i].startswith('VERSION'):
            #print lines[i]
            #print old, new
            lines[i] = lines[i].replace(old, new)
            #print lines[i]
    f.close()
    f = open(fname, 'wt')
    f.writelines(lines)
    f.close()

def main(args):
    
    # Get a list of all filenames that will contain the VERSION variable
    fnames = get_py_files('master_scripts') + get_py_files('scripts') + get_py_files('lib')

    # First, scan all files to make sure they have the VERSION script and that
    # all versions are aligned.
    version_check = [check_version(x) for x in fnames]
    assert all(x is not None for x in version_check), '\nThe following files do not contain a version number:\n' \
            '{}'.format('\n'.join(fnames[i] for i,x in enumerate(version_check) if x is None)) + '\n'
    versions = set(version_check)
    assert len(versions) == 1, '\nMultiple versions detected. Please fix before attempting to modify versions\n' \
            'Versions detected:\n{}'.format('\n'.join(versions)) +'\n'
    curr_version = version_check[0]

    # Parse version file and ensure that the versions match
    version, history = parse_version_file('VERSION.yaml')
    assert curr_version == version, '\nVersion given in all *.py files ({}) does not match VERSION.yaml ({}). ' \
            'Please fix before proceeding.'.format(curr_version, version)

    # Determine new version
    if args.major:
        new_version = increment_version(curr_version, 'major')
        history.insert(0, curr_version)
    elif args.minor:
        new_version = increment_version(curr_version, 'minor')
        history.insert(0, curr_version)
    elif args.patch:
        new_version = increment_version(curr_version, 'patch')
        history.insert(0, curr_version)
    elif args.rollback:
        new_version, history = rollback_version(curr_version, history)
  
    while True:
        proceed = raw_input('\nVersion will be modified from {} to {}.\nDo you wish to proceed? (y/n)' \
                ''.format(curr_version, new_version))
        if proceed in ['y', 'n']:
            break

    if proceed == 'y':
        # Gotta modify version file too!
        write_version_file('VERSION.yaml', new_version, history)

        # Once I know what the new version string is, replace it in all of the files!
        for fname in fnames:
            #print fname
            replace_version(fname, curr_version, new_version)

    # Done


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required = True)
    group.add_argument('-1', '--major', action = 'store_true', help = 'Bump major version number')
    group.add_argument('-2', '--minor', action = 'store_true', help = 'Bump minor version number')
    group.add_argument('-3', '--patch', action = 'store_true', help = 'Bump patch version number')
    group.add_argument('-d', '--rollback', action = 'store_true', help = 'Decrement version to last version in the history')

    args = parser.parse_args()

    main(args)




