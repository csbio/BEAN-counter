#################################################################
######  Copyright: Regents of the University of Minnesota  ######
#################################################################

# Function to open file with any type of compression, or none at all.
# Assumes it is a text file if the "magic bits" do not match.
# Thanks to Ber from StackOverflow for the inspiration!
# http://stackoverflow.com/questions/13044562/python-mechanism-to-identify-compressed-file-type-and-uncompress

import zipfile, bz2, gzip

# factory function to create a suitable instance for accessing files
def get_compressed_file_handle(filename):
    
    # First, read in the very beginning of the file so we can check
    # its magic bytes.
    f = file(filename, 'rb')
    file_start = f.read(1024)
    f.close()

    # Read in the file using the appropirate compressed file reader
    # based on the magic bytes!

    # First, gzip file (most common)
    if file_start.startswith('\x1f\x8b\x08'):
	return gzip.GzipFile(filename, 'rb')
    # Now for bz2 files
    elif file_start.startswith('\x42\x5a\x68'):
	return bz2.BZ2File(filename, 'r')
    # And zip files
    elif file_start.startswith('\x50\x4b\x03\x04'):
	return zipfile.Zipfile(filename, 'r')
    # And if nothing matches, hopefully it's a text file!
    else:
	return file(filename, 'rt')


