#!/usr/bin/env python



# This script must be run from the root of the BEAN-counter source code
# directory. It provides options to bump the major, minor, and patch versions
# of the software within all relevant scripts as well as roll back to the
# latest version in the history.


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

