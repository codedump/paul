#!/usr/bin/python
#
# Copyright (C) 2010 Florin Boariu
#
# This file is part of Paul.
#
# Hooke is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Hooke is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
# Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Hooke.  If not, see
# <http://www.gnu.org/licenses/>.
#
# Based on the IBW loader loadibw module from the Hooke project.
#
#

"""Provides pure Python interface between Elmitec DAT image files
and Numpy arrays.
"""


from paul.base.struct_helper import *
from paul.base.wave import *
from paul.base.errors import *
import numpy, StringIO, curses.ascii

testfile = "C:\\Documents and Settings\\LEEM\\Desktop\\DATA - Permanent\\Pt(111)\\2011-06-01\\Pt111_CO_experiment2_.dat"

__version__ = '0.1'

#
# Converstion BPP -> data type
#
BppTable = {
     8: numpy.int8,
    16: numpy.int16,
    32: numpy.int32,
}

#
# "Data Source" structures. These are basically image meta-data
# structures, called "Data Sources" in the LEEM2k documentation.
#

#
# Unit codes for DS table
#
UnitTable = {
    0: '',
    1:  'V',
    2: 'mA',
    3:  'A',
    4: '^C',
    5: ' K',
    6: 'mV',
    7: 'pA',
    8: 'nA',
    9: 'uA',
}

#
# DS structures
#
DSFloat    = Structure ( name='DSFloat',    fields=[Field('d', 'val')] )
DSFov      = Structure ( name='DSFov',      fields=[Field('c', 'fov', count=17)] )
DSCamExp   = Structure ( name='DSCamExp',   fields=[Field('d', 'ms')] )
DSTitle    = Structure ( name='DSTitle',    fields=[Field('c', 'title', count=17)] )
DSMitutoyo = Structure ( name='DSMitutoyo', fields=[Field('d', 'val1'),
                                                    Field('d', 'val2')] )
DSFovCal   = Structure ( name='DSFovCal',   fields=[Field('c', 'text', count=17),
                                                    Field('d', 'val')] )
DSVarian   = Structure ( name='DSVarian',   fields=[Field('c', 'label', count=17),
                                                    Field('d', 'val'),
						    Field('c', 'unit', count=5)] )
#DSGeneric  = Structure ( name='DSGeneric',  fields=[]] )

#
# MI codes. The metadata buffer is made up of
#    - 1 MI code,
#    - followed by the corresponding MI structure,
#
DataSources = {
 # 0..99: Generic information, with the following format:
 #         - LEEM Module ID (16 bit int?)
 #         - Followed by name (variable-length string
 #         - Followed by unit code (see UnitTable)
 #         - 0 (string termination)
 #         - 1 float (data)
    100: DSMitutoyo,
    101: DSFov,
    102: DSFloat, # old varian
    103: DSFloat, # old varian
    104: DSCamExp,
    105: DSTitle,
    106: DSVarian, 
    107: DSVarian, 
    108: DSVarian, 
    109: DSVarian,
    110: DSFovCal,
    111: DSFloat, # phi, theta
    112: DSFloat,
   #255: Skip!
}

#
# DAT: File Format
#
#  FileHeader7
#  RecipeBlock (if FileHeader7.attachedRecipeSize > 0)
#  ImageHeader
#  Image Data (width * height, 16 bpp, grayscale).
#  ImageMarkup block
#

#
# This is fiel header version 7, as documented by Elmitec.
# Case of field names has been changed to keep the coding style
# consistent (in the original version, upper and lower case
# starting letters were mixed).
#
FileHeader7 = Structure(
    name='FileHeader7',
    fields=[
        Field('c', 'id', help='Version specification', count=20),
        Field('h', 'size', help=''),
        Field('h', 'version', help=''),
        Field('h', 'bitsPerPixel', help='Storage size of the pixels in bits.'),
        Field('c', 'spare0', help='', count=14),
        Field('h', 'imageWidth', help=''),
        Field('h', 'imageHeight', help=''),
        Field('h', 'nrImages', help=''),
        Field('h', 'attachedRecipeSize', help='Size of attached recipe data (up to 128 bytes).'),
        Field('c', 'spare2', help='', count=56),
        ])

#
# 128 of sequencer recipe block. If FileHeader7.attachedRecipeSize > 0,
# then there is a 128 byte recipe attached, regardless of the specified
# size. If specified size is < 128, then only the specified number of bytes
# actually contains useful data.
#
Recipe = Structure(
    name = 'Recipe',
    fields=[
        Field('c', 'recipe', help='Recipe data, see Elmitec docs for format info.', count=128),
        ])

#
# The image header.
#
ImageHeader = Structure(
    name = 'ImageHeader',
    fields = [
        Field('h', 'size', help=''),
        Field('h', 'version', help=''),
        Field('h', 'colorScaleLow', help=''),
        Field('h', 'colorScaleHigh', help=''),
        Field('Q', 'imageTime', help='Standard Windows FILETIME '
                   '(64 bit value with the number of 100ns units since 1. January 1601).'),
        Field('h', 'maskXShift',  help=''),
        Field('h', 'maskYShift',  help=''),
        Field('c', 'useMask', help=''),
        Field('c', 'spare0', help=''),
        Field('h', 'attachedMarkupSize', help='Video files (DAV) never contain markups.'),
        Field('h', 'spin', help=''),
        Field('h', 'leemDataVersion', help='Currently version 2.'),
        # LEEM Data structure field was moved in a separate structure in the Python implementation
        ])

#
# ImageMarkup block (attached to ImageHeader+LeemData if ImageHeader.attacehdMarkupSize > 0
#
ImageMarkup = Structure(
    name = 'ImageMarkup',
    fields = [
        Field ('h', 'data', help='Infos about lines & markers (128 bytes of data, '
                           'present if ImageHeader.attachedMarkupSize > 0', count=128),
        ])

#
# Loads the next LEEM 2000 data block from the file stream 'file'
# and returns it in a usable form:
#  {'key': [value, 'unit', ...]}
#
def getNextDs (io):
    id = struct.unpack('b', io.read(1))[0]
    if (id == 255):
        # buffer ends here, no more meta info
        return {}
    elif (id >= 0 and id < 99):
        # generic structure, variable length
	ds = getGenericDs(io)
	print ds
	getNextDs(io)
    else:
        # typed structure, one of the DS applies
        print "Typed DS: ", DataSources[id].name
	dat = io.read(DataSources[id].size)
        ds = {}
	getNextDs(io)
    return ds

#
# Reads exactly one generic DS from specified buffer,
# returns the data and the buffer.
# The format of a (variable-length) buffer is the following:
#         - Unit name (variable-length string)
#         - Unit code, as an ASCII digit (see UnitTable)
#         - 0 (byte value string termination)
#         - 1 float (data)
#
def getGenericDs (io):
    mod_name = ''
    next_char = io.read(1)
    next_byte = struct.unpack('b', next_char)[0]
    while ((next_byte > 31 and next_byte < 128) and not curses.ascii.isdigit(next_char)):
        ch = struct.unpack('c', next_char)
	mod_name += ch[0]
	next_char = io.read(1)
	next_byte = struct.unpack('b', next_char)[0]
    mod_unit = UnitTable[int(next_char)]
    io.read(1)  # this is the trailing 0 of the name+unit string
    mod_value = struct.unpack('f', io.read(4))[0]
    return { mod_name: [mod_value, mod_unit]}

#
# Loads a DAT image file, returns the data and header information
#
def load(filename):

    if hasattr(filename, 'read'):
        f = filename  # filename is actually a stream object
    else:
        f = open(filename, 'rb')

    try:
        fileh = FileHeader7.unpack_dict_from (buffer(f.read(FileHeader7.size)))

        magic = reduce (lambda x, y: x+y, fileh['id'])
        if (magic != "UKSOFT2001"):
            raise FormatError ("\"%s\" doesn\'t seem to be an Elmitec DAT file (wrong magic)."
                               % filename)
        if (fileh['version'] < 7):
            raise VersionError ('Wrong file format: version %d found, >7 expected'
                                % fileh['version'])

        if (fileh['attachedRecipeSize'] > 0):
            recipe = Recipe.unpack_dict_from (buffer(f.read(Recipe.size)))                                             
        imgh = ImageHeader.unpack_dict_from (buffer(f.read(ImageHeader.size)))
        leemd = buffer(f.read(256))
        raw = buffer(f.read(fileh['imageWidth']*fileh['imageHeight']*fileh['bitsPerPixel']/8))
        data = numpy.ndarray (buffer=raw, shape=[fileh['imageWidth'], fileh['imageHeight']],
                              dtype=BppTable[fileh['bitsPerPixel']])

	io = StringIO.StringIO(buf)
        
    finally:
        if not hasattr(filename, 'read'):
            f.close()

    return data, fileh, imgh, leemd


def save (filename):
    raise NotImplementedError

if __name__ == '__main__':
    """DAT -> ASCII conversion
    """
    import optparse
    import sys

    load (testfile)
