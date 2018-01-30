#!/usr/bin/env python

import os

# Module for function that adds/updates the "VERSION" file in a specified output directory.

def update_version_file(outdir, version):
    version_file = os.path.join(outdir, 'VERSION')
    match = False
    with open(version_file, 'a+') as f:
        for line in f:
            if version == line.rstrip():
                match = True
                return None
        if not match:
            f.write(version + '\n')
    return None


