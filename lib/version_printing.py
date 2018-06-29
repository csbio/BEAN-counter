#!/usr/bin/env python

VERSION='2.5.1'

import os

# Module for function that adds/updates the "VERSION" file in a specified output directory.

def update_version_file(outdir, version):
    version_file = os.path.join(outdir, 'VERSION')
    header = 'BEAN-counter versions used to generate results in this folder:'
    header_match = False
    version_match = False
    with open(version_file, 'a+') as f:
        for line in f:
            if header == line.rstrip():
                header_match = True
            if version == line.rstrip():
                version_match = True
        if not header_match:
            f.write(header + '\n')
        if not version_match:
            f.write(version + '\n')
    return None


