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
# Returns a DS string (variable-length, C-style string).
# Maximum length is the length of the string without the
# trailing 0.
#
def getDsString(io, max=-1):
    s = ''
    next_byte = struct.unpack('B', io.read(1))[0]
    while (next_byte):
	if (not curses.ascii.iscntrl(next_byte)):
            s += chr(next_byte)
	next_byte = struct.unpack('B', io.read(1))[0]

    if (max > -1 and len(s) > max):
	    raise FormatError("String '%s' too long" % s)
    return s

#
# Reads exactly one generic DS from specified buffer,
# returns the data and the buffer.
# The format of a (variable-length) buffer is the following:
#         - Unit name (variable-length string)
#         - Unit code, as an ASCII digit (see UnitTable)
#         - 0 (byte value string termination)
#         - 1 float (data)
#
def getDsGeneric(io):
    unit = getDsString(io)
    mod_name = unit[:-1]
    mod_unit = UnitTable[int(unit[-1:])]
    mod_value = struct.unpack('f', io.read(4))[0]
    return [ { mod_name: [mod_value, mod_unit] } ]

#
# Reads a FovCal chunk. Format:
#   . string + trailing 0
#   . float
#
def getDsFovCal(io):
    name = getDsString(io, max=16)
    if (name == 'none'):
	# make things prettier: for some unknown reason,
        # FOV name is most of the time 'none'.
        name = "FOV cal."
    return [ { name: [struct.unpack('f', io.read(4))[0]] } ]

# Reads cam exposure info (single float)
def getDsCamExp(io):
    return [ { 'Exposure': [struct.unpack('f', io.read(4))[0], 's'] } ]



# Reads Mitutotyo micrometers (2 floats for X and Y axis).
def getDsMitutoyo(io):
    return [ { 'Pos X': [struct.unpack('f', io.read(4))[0], 'um'] },
	     { 'Pos Y': [struct.unpack('f', io.read(4))[0], 'um'] } ]

# Reads title string (Variable-length C string).
def getDsTitle(io):
    return [ { 'Title': [getDsString(io, max=16)] }]


# Phi, Theta (two floats)
def getDsPhiTheta(io):
    return [ { 'Phi':   [struct.unpack('f', io.read(4))[0]] },
	     { 'Theta': [struct.unpack('f', io.read(4))[0]] } ]

# Spin (float)
def getDsSpin(io):
    return [ { 'Spin': [struct.unpack('f', io.read(4))[0]] } ]


# Varian controller info. Format:
#    . Name (string 17)
#    . Unit (string 5)
#    . Float
def getDsVarian(io):
    name = getDsString(io, max=16)
    unit = getDsString(io, max=16)
    return [ { name: [struct.unpack('f', io.read(4))[0], unit] } ]

# Skip entry, called for every entry above DataSources['max'].
# Returns an empty list.
def getDsSkip(io):
    return []

#
# MI codes. The metadata buffer is made up of
#    - 1 MI code,
#    - followed by the corresponding MI structure,
#
DataSources = {
      0: getDsGeneric,
    100: getDsMitutoyo,
    101: DSFov,   # old FOV
    102: DSFloat, # old varian
    103: DSFloat, # old varian
    104: getDsCamExp,
    105: getDsTitle,
    106: getDsVarian,
    107: getDsVarian, 
    108: getDsVarian, 
    109: getDsVarian,
    110: getDsFovCal,
    111: getDsPhiTheta, # phi, theta
    112: getDsSpin,
   'max': 112
}

# some automatic filling of DataSources stuff...
for i in range(100):
    DataSources[i]=DataSources[0]
for i in range(DataSources['max'], 256):
    DataSources[i]=getDsSkip

#
# Loads the next LEEM 2000 data block from the file stream 'file'
# and returns it in a usable form:
#  {'key': [value, 'unit', ...]}
#
def getDs (io):
    dslist = []
    id_str = io.read(1)
    while (len(id_str) > 0):
        id = struct.unpack('B', id_str)[0]
	dslist += DataSources[id](io)
        id_str = io.read(1)
    return dslist

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
        leemd_buf = buffer(f.read(256))
        raw = buffer(f.read(fileh['imageWidth']*fileh['imageHeight']*fileh['bitsPerPixel']/8))
        data = numpy.ndarray (buffer=raw, shape=[fileh['imageWidth'], fileh['imageHeight']],
                              dtype=BppTable[fileh['bitsPerPixel']])
	leemd = getDs(StringIO.StringIO(leemd_buf))
        
    finally:
        if not hasattr(filename, 'read'):
            f.close()

    return data, leemd, fileh, imgh, leemd_buf


def save (filename):
    raise NotImplementedError

if __name__ == '__main__':
    """DAT -> ASCII conversion
    """
    import optparse
    import sys

    load (testfile)
