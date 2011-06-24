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
import numpy

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
# Unit codes for LEEMData table
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
# LEEM2000 module hashes
#
DataSourceTable = {
    100: 1,  # Mitutoyo,           2 floats
    101: 2,  # FOV,                string 17 bytes
    102: 3,  # varian #1.1,        1 float
    103: 4,  # varian #1.2,        1 float
    104: 5,  # cam exposure in ms, 1 float
    105: 6,  # title,              string 17 bytes
    106: 7,  # varian #1.1         string 17 (label), string 5 (units), 1 float
    107: 8,  # varian #1.2         string 17 (label), string 5 (units), 1 float
    108: 9,  # varian #2.1         string 17 (label), string 5 (units), 1 float
    109: 10, # varian #2.2         string 17 (label), string 5 (units), 1 float
    110: 11, # FOV cal. factor     string 17, 1 float
    111: 12, # phi, theta          1 float, 1 float
    112: 13, # spin
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
        Field('c', 'recipe', help='Recipe data, see Elmitec docs for format specification.', count=128),
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
        Field('Q', 'imageTime', help='Standard Windows FILETIME (64 bit value with the number of 100ns units since 1. Januarz 1601).'),
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
# LEEM data structure, originally part of the ImageHeader, is a 256 bytes long
# unsigned char array, constructed as follows:
#   . Source1, Argument1, Source2, Argument2, ..., SourceN, ArgumentN
#
# 0..99:  LEEM tags, as follows:
#           - leem module id
#           - followed by name
#           - followed by 1 ASCII unit code (0=none, 1=V, 2=mA, 3=A, 4=C, 5=K, 6=mV, 7=pA, 8=nA, 9=uA)
#           - \0 (string termination)
#           - data (1 float = 4 bytes)
#
LeemData = Structure(
    name = 'LeemData',
    fields = [
        Field('c', 'modString', help='String containing LEEM2000 module information.', count=100),
        ])

#        Field('f', 'posX', help='Mitutoyo X position value'),
#        Field('f', 'posY', help='Mitutoyo Y position value'),
#        Field('f', 'exposure', help='Camera exposure time in milliseconds.'),
#        Field('c', 'title', help='Title', count=17),
#
#        Field('c', 'pGauge1Label', help='Pressure gauge #1 label.', count=17),
#        Field('c', 'pGauge1Unit', help='Pressure gauge #1 units.',  count=5),
#        Field('f', 'pGauge1Value', help='Pressure gauge #1 value.'),
#
#        Field('c', 'pGauge2Label', help='Pressure gauge #2 label.', count=17),
#        Field('c', 'pGauge2Unit', help='Pressure gauge #2 units.',  count=5),
#        Field('f', 'pGauge2Value', help='Pressure gauge #2 value.'),
#
#        Field('c', 'pGauge3Label', help='Pressure gauge #3 label.', count=17),
#        Field('c', 'pGauge3Unit', help='Pressure gauge #3 units.',  count=5),
#        Field('f', 'pGauge3Value', help='Pressure gauge #3 value.'),
#
#        Field('c', 'pGauge4Label', help='Pressure gauge #4 label.', count=17),
#        Field('c', 'pGauge4Unit', help='Pressure gauge #4 units.',  count=5),
#        Field('f', 'pGauge4Value', help='Pressure gauge #4 value.'),
#
#        Field('c', 'fov', help='Camera FOV.', count=17),
#        Field('f', 'fovFactor', help='FOV calibration factor.', count=17),
#        Field('f', 'phi', help=''),
#        Field('f', 'theta', help=''),
#        Field('f', 'spin', help=''),
#        ])

#
# ImageMarkup block (attached to ImageHeader+LeemData if ImageHeader.attacehdMarkupSize > 0
#
ImageMarkup = Structure(
    name = 'ImageMarkup',
    fields = [
        Field ('h', 'data', help='Infos about lines & markers (128 bytes of data, present if ImageHeader.attachedMarkupSize > 0', count=128),
        ])
    

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
        leemd = LeemData.unpack_dict_from (buffer(f.read(LeemData.size)))

        raw = buffer(f.read(fileh['imageWidth']*fileh['imageHeight']*fileh['bitsPerPixel']/8))
        data = numpy.ndarray (buffer=raw, shape=[fileh['imageWidth'], fileh['imageHeight']],
                              dtype=BppTable[fileh['bitsPerPixel']])
        
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
