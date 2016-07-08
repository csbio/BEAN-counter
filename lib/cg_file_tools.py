#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

import os, sys

# Random tools for manipulating files in the cg pipeline
def create_output_dir(path):

    if not os.path.isdir(path):
        os.makedirs(path)

    return None
