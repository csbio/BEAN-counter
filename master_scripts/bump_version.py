#!/usr/bin/env python

import yaml

# This script must be run from the root of the BEAN-counter source code
# directory. It provides options to bump the major, minor, and patch versions
# of the software within all relevant scripts as well as roll back to the
# latest version in the history.

def parse_version_file(fname):
    with open(fname, 'rt') as f:
        version_dict = yaml.load(f)

    assert 'current' in version_dict, 'VERSION.yaml must have a key named "current".'
    assert is_valid_version(version_dict['current']), 'The version string "{}" is not a valid version (must be "X.Y.Z").'.format(version_dict['current'])

    history = version_dict.get('history')
    if history is None or history == '':
        history = []
    history = history.sort()
    history = history[::-1]
    invalid_history_strings = [x for x in history if not is_valid_version(x)]
    assert not any(invalid_history_strings), 'The current versions in the histroy of VERSION.yaml are invalid:\n' \
            '{}\nPlease correct before moving on.'.format('\n'.join(invalid_history_strings))
    return version_dict['current'], history

def parse_version(x):
    '''
    Assumes x has been stripped of whitespace, if read in from a file for
    example.
    '''
    split_x = [int(xx) for xx in x.split('.')]
    assert len(split_x) == 3, 'Version string "{}" is not in the proper X.Y.Z format'.format(x)
    return split_x

def increment_version(vers, Type):
    major, minor, patch = parse_version(vers)
    if Type == 'patch':
        patch += 1
    elif type == 'minor':
        patch = 0
        minor += 1
    else:
        patch = 0
        minor = 0
        major += 1

    return '.'.join([str(major), str(minor), str(patch)])

def rollback_version(vers, history):

    if len(history) > 0:
        return history[0]
    else:
        return vers

def main(args):
    pass

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--major', action = 'store_true', help = 'Bump major version number')
    group.add_argument('--minor', action = 'store_true', help = 'Bump minor version number')
    group.add_argument('--patch', action = 'store_true', help = 'Bump patch version number')
    group.add_argument('-d', '--rollback', action = 'store_true', help = 'Decrement version to last version in the history')

    args = parser.parse_args()

    main(args)




