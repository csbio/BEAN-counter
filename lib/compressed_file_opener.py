# Function to open file with any type of compression, or none at all.
# Assumes it is a text file if the "magic bits" do not match.
# Currently all code courtesy of Ber from StackOverflow

class CompressedFile (object):
    magic = None
    file_type = None
    mime_type = None
    proper_extension = None

    def __init__(self, f):
        # f is an open file or file like object
        self.f = f
        self.accessor = self.open()

    @classmethod
    def is_magic(self, data):
        return data.startswith(self.magic)

    def open(self):
        return None

import zipfile

class ZIPFile (CompressedFile):
    magic = '\x50\x4b\x03\x04'
    file_type = 'zip'
    mime_type = 'compressed/zip'

    def open(self):
        return zipfile.ZipFile(self.f)

import bz2

class BZ2File (CompressedFile):
    magic = '\x42\x5a\x68'
    file_type = 'bz2'
    mime_type = 'compressed/bz2'

    def open(self):
        return bz2.BZ2File(self.f)

import gzip

class GZFile (CompressedFile):
    magic = '\x1f\x8b\x08'
    file_type = 'gz'
    mime_type = 'compressed/gz'

    def open(self):
        return gzip.GzipFile(self.f)

class TxtFile (CompressedFile):

    file_type = 'plain'
    mime_type = 'text/plain'

    def open(self):
        return file(self.f)

# factory function to create a suitable instance for accessing files
def get_compressed_file(filename):
    with file(filename, 'rb') as f:
        start_of_file = f.read(1024)
        # f.seek(0)
        f.close()
        for cls in (ZIPFile, BZ2File, GZFile):
            if cls.is_magic(start_of_file):
                return cls(filename)
            
        return TxtFile(filename)

