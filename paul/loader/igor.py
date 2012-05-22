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
#
# Based on code from the Hooke project, by Trevor King, extended
# to support IBW writing aswell.
#

"""Provides pure Python interface between IGOR Binary
Wave files and Numpy arrays.
"""

# Based on WaveMetric's Technical Note 003, "Igor Binary Format"
#   ftp://ftp.wavemetrics.net/IgorPro/Technical_Notes/TN003.zip
# From ftp://ftp.wavemetrics.net/IgorPro/Technical_Notes/TN000.txt
#   We place no restrictions on copying Technical Notes, with the
#   exception that you cannot resell them. So read, enjoy, and
#   share. We hope IGOR Technical Notes will provide you with lots of
#   valuable information while you are developing IGOR applications.

import logging
log = logging.getLogger (__name__)

from paul.base.struct_helper import *
from paul.base.wave import Wave
from paul.base.errors import *
import os

__version__ = '0.1'


# Numpy doesn't support complex integers by default, see
#   http://mail.python.org/pipermail/python-dev/2002-April/022408.html
#   http://mail.scipy.org/pipermail/numpy-discussion/2007-October/029447.html
# So we roll our own types.  See
#   http://docs.scipy.org/doc/numpy/user/basics.rec.html
#   http://docs.scipy.org/doc/numpy/reference/generated/numpy.dtype.html
complexInt8 = numpy.dtype([('real', numpy.int8), ('imag', numpy.int8)])
complexInt16 = numpy.dtype([('real', numpy.int16), ('imag', numpy.int16)])
complexInt32 = numpy.dtype([('real', numpy.int32), ('imag', numpy.int32)])
complexUInt8 = numpy.dtype([('real', numpy.uint8), ('imag', numpy.uint8)])
complexUInt16 = numpy.dtype([('real', numpy.uint16), ('imag', numpy.uint16)])
complexUInt32 = numpy.dtype([('real', numpy.uint32), ('imag', numpy.uint32)])


# Begin IGOR constants and typedefs from IgorBin.h

# From IgorMath.h
TYPE_TABLE = {       # (key: integer flag, value: numpy dtype)
    0:None,          # Text wave, not handled in ReadWave.c
    1:numpy.complex, # NT_CMPLX, makes number complex.
    2:numpy.float32, # NT_FP32, 32 bit fp numbers.
    3:numpy.complex64,
    4:numpy.float64, # NT_FP64, 64 bit fp numbers.
    5:numpy.complex128,
    8:numpy.int8,    # NT_I8, 8 bit signed integer. Requires Igor Pro
                     # 2.0 or later.
    9:complexInt8,
    0x10:numpy.int16,# NT_I16, 16 bit integer numbers. Requires Igor
                     # Pro 2.0 or later.
    0x11:complexInt16,
    0x20:numpy.int32,# NT_I32, 32 bit integer numbers. Requires Igor
                     # Pro 2.0 or later.
    0x21:complexInt32,
#   0x40:None,       # NT_UNSIGNED, Makes above signed integers
#                    # unsigned. Requires Igor Pro 3.0 or later.
    0x48:numpy.uint8,
    0x49:complexUInt8,
    0x50:numpy.uint16,
    0x51:complexUInt16,
    0x60:numpy.uint32,
    0x61:complexUInt32,
}

# From wave.h
MAXDIMS = 4

# From binary.h
BinHeaderCommon = Structure(  # WTK: this one is mine.
    name='BinHeaderCommon',
    fields=[
        Field('h', 'version', help='Version number for backwards compatibility.'),
        ])

BinHeader1 = Structure(
    name='BinHeader1',
    fields=[
        Field('h', 'version', help='Version number for backwards compatibility.'),
        Field('l', 'wfmSize', help='The size of the WaveHeader2 data structure plus the wave data plus 16 bytes of padding.'),
        Field('h', 'checksum', help='Checksum over this header and the wave header.'),
        ])

BinHeader2 = Structure(
    name='BinHeader2',
    fields=[
        Field('h', 'version', help='Version number for backwards compatibility.'),
        Field('l', 'wfmSize', help='The size of the WaveHeader2 data structure plus the wave data plus 16 bytes of padding.'),
        Field('l', 'noteSize', help='The size of the note text.'),
        Field('l', 'pictSize', default=0, help='Reserved. Write zero. Ignore on read.'),
        Field('h', 'checksum', help='Checksum over this header and the wave header.'),
        ])

BinHeader3 = Structure(
    name='BinHeader3',
    fields=[
        Field('h', 'version', help='Version number for backwards compatibility.'),
        Field('h', 'wfmSize', help='The size of the WaveHeader2 data structure plus the wave data plus 16 bytes of padding.'),
        Field('l', 'noteSize', help='The size of the note text.'),
        Field('l', 'formulaSize', help='The size of the dependency formula, if any.'),
        Field('l', 'pictSize', default=0, help='Reserved. Write zero. Ignore on read.'),
        Field('h', 'checksum', help='Checksum over this header and the wave header.'),
        ])

BinHeader5 = Structure(
    name='BinHeader5',
    fields=[
        Field('h', 'version', help='Version number for backwards compatibility.'),
        Field('h', 'checksum', help='Checksum over this header and the wave header.'),
        Field('l', 'wfmSize', help='The size of the WaveHeader5 data structure plus the wave data.'),
        Field('l', 'formulaSize', help='The size of the dependency formula, if any.'),
        Field('l', 'noteSize', help='The size of the note text.'),
        Field('l', 'dataEUnitsSize', help='The size of optional extended data units.'),
        Field('l', 'dimEUnitsSize', help='The size of optional extended dimension units.', count=MAXDIMS),
        Field('l', 'dimLabelsSize', help='The size of optional dimension labels.', count=MAXDIMS),
        Field('l', 'sIndicesSize', help='The size of string indicies if this is a text wave.'),
        Field('l', 'optionsSize1', default=0, help='Reserved. Write zero. Ignore on read.'),
        Field('l', 'optionsSize2', default=0, help='Reserved. Write zero. Ignore on read.'),
        ])


# From wave.h
MAX_WAVE_NAME2 = 18 # Maximum length of wave name in version 1 and 2
                    # files. Does not include the trailing null.
MAX_WAVE_NAME5 = 31 # Maximum length of wave name in version 5
                    # files. Does not include the trailing null.
MAX_UNIT_CHARS = 3

# Header to an array of waveform data.

WaveHeader2 = Structure(
    name='WaveHeader2',
    fields=[
        Field('h', 'type', help='See types (e.g. NT_FP64) above. Zero for text waves.'),
        Field('P', 'next', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('c', 'bname', help='Name of wave plus trailing null.', count=MAX_WAVE_NAME2+2),
        Field('h', 'whVersion', default=0, help='Write 0. Ignore on read.'),
        Field('h', 'srcFldr', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('P', 'fileName', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('c', 'dataUnits', default=0, help='Natural data units go here - null if none.', count=MAX_UNIT_CHARS+1),
        Field('c', 'xUnits', default=0, help='Natural x-axis units go here - null if none.', count=MAX_UNIT_CHARS+1),
        Field('l', 'npnts', help='Number of data points in wave.'),
        Field('h', 'aModified', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('d', 'hsA', help='X value for point p = hsA*p + hsB'),
        Field('d', 'hsB', help='X value for point p = hsA*p + hsB'),
        Field('h', 'wModified', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('h', 'swModified', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('h', 'fsValid', help='True if full scale values have meaning.'),
        Field('d', 'topFullScale', help='The min full scale value for wave.'), # sic, 'min' should probably be 'max'
        Field('d', 'botFullScale', help='The min full scale value for wave.'),
        Field('c', 'useBits', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('c', 'kindBits', default=0, help='Reserved. Write zero. Ignore on read.'),
        Field('P', 'formula', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('l', 'depID', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('L', 'creationDate', help='DateTime of creation.  Not used in version 1 files.'),
        Field('c', 'wUnused', default=0, help='Reserved. Write zero. Ignore on read.', count=2),
        Field('L', 'modDate', help='DateTime of last modification.'),
        Field('P', 'waveNoteH', help='Used in memory only. Write zero. Ignore on read.'),
        Field('f', 'wData', help='The start of the array of waveform data.', count=4),
        ])

WaveHeader5 = Structure(
    name='WaveHeader5',
    fields=[
        Field('P', 'next', help='link to next wave in linked list.'),
        Field('L', 'creationDate', help='DateTime of creation.'),
        Field('L', 'modDate', help='DateTime of last modification.'),
        Field('l', 'npnts', help='Total number of points (multiply dimensions up to first zero).'),
        Field('h', 'type', help='See types (e.g. NT_FP64) above. Zero for text waves.'),
        Field('h', 'dLock', default=0, help='Reserved. Write zero. Ignore on read.'),
        Field('c', 'whpad1', default=0, help='Reserved. Write zero. Ignore on read.', count=6),
        Field('h', 'whVersion', default=1, help='Write 1. Ignore on read.'),
        Field('c', 'bname', help='Name of wave plus trailing null.', count=MAX_WAVE_NAME5+1),
        Field('l', 'whpad2', default=0, help='Reserved. Write zero. Ignore on read.'),
        Field('P', 'dFolder', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        # Dimensioning info. [0] == rows, [1] == cols etc
        Field('l', 'nDim', help='Number of of items in a dimension -- 0 means no data.', count=MAXDIMS),
        Field('d', 'sfA', help='Index value for element e of dimension d = sfA[d]*e + sfB[d].', count=MAXDIMS),
        Field('d', 'sfB', help='Index value for element e of dimension d = sfA[d]*e + sfB[d].', count=MAXDIMS),
        # SI units
        Field('c', 'dataUnits', default=0, help='Natural data units go here - null if none.', count=MAX_UNIT_CHARS+1),
        Field('c', 'dimUnits', default=0, help='Natural dimension units go here - null if none.', count=(MAXDIMS, MAX_UNIT_CHARS+1)),
        Field('h', 'fsValid', help='TRUE if full scale values have meaning.'),
        Field('h', 'whpad3', default=0, help='Reserved. Write zero. Ignore on read.'),
        Field('d', 'topFullScale', help='The max and max full scale value for wave'), # sic, probably "max and min"
        Field('d', 'botFullScale', help='The max and max full scale value for wave.'), # sic, probably "max and min"
        Field('P', 'dataEUnits', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('P', 'dimEUnits', default=0, help='Used in memory only. Write zero.  Ignore on read.', count=MAXDIMS),
        Field('P', 'dimLabels', default=0, help='Used in memory only. Write zero.  Ignore on read.', count=MAXDIMS),
        Field('P', 'waveNoteH', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('l', 'whUnused', default=0, help='Reserved. Write zero. Ignore on read.', count=16),
        # The following stuff is considered private to Igor.
        Field('h', 'aModified', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('h', 'wModified', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('h', 'swModified', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('c', 'useBits', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('c', 'kindBits', default=0, help='Reserved. Write zero. Ignore on read.'),
        Field('P', 'formula', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('l', 'depID', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('h', 'whpad4', default=0, help='Reserved. Write zero. Ignore on read.'),
        Field('h', 'srcFldr', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('P', 'fileName', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('P', 'sIndices', default=0, help='Used in memory only. Write zero. Ignore on read.'),
        Field('f', 'wData', help='The start of the array of data.  Must be 64 bit aligned.', count=1),
        ])

# End IGOR constants and typedefs from IgorBin.h

# headers for packd files
PackedFileRecordHeader = Structure (
    name = "PackedFileRecord",
    fields = [
        Field ('H','recordType', default=0, help='Record type plus superceeded flag.'),
        Field ('h','version', default=0, help='Version information, record type dependent.'),
        Field ('l','numDataBytes', default=0, help='Size of the record.'),
        ])

PackedFileRecordType = {
    0: 'Unused',            #
    1: 'Variables',         # Contains system numeric variables
    2: 'History',           # Contains the experiment's history as plain text
    3: 'Wave',              # Contains the data for a wave
    4: 'Recreation',        # Contains the experiment's recreation procedures as plain text.
    5: 'Procedure',         # Contains the experiment's main procedure window text as plain text.
    6: 'Unused2',           # 
    7: 'GetHistory',        # Not a real record but rather, a message to go back and read the history text.
    8: 'PackedFile',        # Contains the data for a procedure file or notebook in a packed form.
    9: 'DataFolderStart',   # Marks the start of a new data folder.
    10: 'DataFolderEnd',    # Marks the end of a data folder.
    11: 'Reserved'          # From Igor docs:
                            #   "Igor writes other kinds of records in a packed experiment file, for storing
                            #    things like pictures, page setup records, and miscellaneous settings. The
                            #    format for these records is quite complex and is not described in PTN003.
                            #    If you are writing a program to read packed files, you must skip any record
                            #    with a record type that is not listed above."
}

# invert
PackedFileRecordId = dict((v,k) for k, v in PackedFileRecordType.iteritems())

def need_to_reorder_bytes(version):
    # If the low order byte of the version field of the BinHeader
    # structure is zero then the file is from a platform that uses
    # different byte-ordering and therefore all data will need to be
    # reordered.
    return version & 0xFF == 0

def byte_order(needToReorderBytes):
    little_endian = sys.byteorder == 'little'
    if needToReorderBytes:
        little_endian = not little_endian
    if little_endian:
        return '<'  # little-endian
    return '>'  # big-endian    

def version_structs(version, byte_order):
    if version == 1:
        bin = BinHeader1
        wave = WaveHeader2
    elif version == 2:
        bin = BinHeader2
        wave = WaveHeader2
    elif version == 3:
        bin = BinHeader3
        wave = WaveHeader2
    elif version == 5:
        bin = BinHeader5
        wave = WaveHeader5
    else:
        raise FormatError ('This does not appear to be a valid Igor binary wave file.'
                           ' The version field = %d.' % version)
    checkSumSize = bin.size + wave.size
    if version == 5:
        checkSumSize -= 4  # Version 5 checksum does not include the wData field.
    bin.set_byte_order(byte_order)
    wave.set_byte_order(byte_order)
    return (bin, wave, checkSumSize)

def checksum(buffer, byte_order, oldcksum, numbytes):
    x = numpy.ndarray(
        (numbytes/2,), # 2 bytes to a short -- ignore trailing odd byte
        dtype=numpy.dtype(byte_order+'h'),
        buffer=buffer)
    oldcksum += x.sum()
    if oldcksum > 2**31:  # fake the C implementation's int rollover
        oldcksum %= 2**32
        if oldcksum > 2**31:
            oldcksum -= 2**31
    return oldcksum & 0xffff


def load(filename):
    '''
    Loads an IBW from specified file, returns a Wave() object
    (which is an overloaded ndarray) containing the Igor wave data.

    This function is a wrapper for either wave_read(), which will transparently
    handle loading waves from either proper .IBW files, on disk, or from Igor
    packed files (.PXP or .PXT files), following the Igor path
    format inside the packed file (i.e.: ":this:is:an:igor:path" representing
    a wave called "path", at the path component /this/is/an/igor/)
    '''
    return wave_read(filename)


def wave_read_header(filename):
    '''
    Reads the wave header information. This is useful for quick access to a wave's
    data (i.e. the name?) without having to load the complete wave. Nice for 
    big waves and browsing PXP files.
    '''

    if hasattr(filename, 'read'):
        f = filename  # filename is actually a stream object
    else:
        f = open(filename, 'rb')
    try:

        b = buffer(f.read(BinHeaderCommon.size))
        version = BinHeaderCommon.unpack_dict_from(b)['version']
        needToReorderBytes = need_to_reorder_bytes(version)
        byteOrder = byte_order(needToReorderBytes)
    
        if needToReorderBytes:
            BinHeaderCommon.set_byte_order(byteOrder)
            version = BinHeaderCommon.unpack_dict_from(b)['version']
        bin_struct, wave_struct, checkSumSize = version_structs(version, byteOrder)

        b = buffer(b + f.read(bin_struct.size + wave_struct.size - BinHeaderCommon.size))
        c = checksum(b, byteOrder, 0, checkSumSize)
        if c != 0:
            raise VersionError('load: %s: error in checksum - should be 0, is %d.  '
                               'This does not appear to be a valid Igor binary wave file.'
                               % (filename, c))
        bin_info = bin_struct.unpack_dict_from(b)
        wave_info = wave_struct.unpack_dict_from(b, offset=bin_struct.size)

        #
        # Need to calculate this because we need access to the Structure class
        # that unpacked the header.
        #
        # We are creating some new dictionary items:
        #   'wave_data_size': represents the size of the wave data
        #                     (without tail or header info)
        #   'tail_size': the size of the wdata field in the correct WaveHeader
        #                structure
        #   'tail_data': the corresponding tail data
        #   'byte_order': byte order parameter as returned by byteorder() above
        #
        # and pass those that along with the bin_info dictionary.
        #
        if version in [1,2,3]:
            bin_info['tail_size'] = 16
            bin_info['wave_data_size'] = bin_info['wfmSize'] - wave_struct.size
        else:
            assert version == 5, version
            bin_info['tail_size'] = 4
            bin_info['wave_data_size'] = bin_info['wfmSize'] - (wave_struct.size - bin_info['tail_size'])

        bin_info['tail_data'] = numpy.array(b[-bin_info['tail_size']:])
        bin_info['byte_order'] = byteOrder
        
        return wave_info, bin_info

    finally:
        if not hasattr(filename, 'read'):
            log.debug ("wave_read_header: Closing %s" % filename)
            f.close()


def wave_read_data (f, wave_info, bin_info):
    '''
    Reads main wave data
    '''
        
    # dtype() wrapping to avoid numpy.generic and
    # getset_descriptor issues with the builtin Numpy types
    # (e.g. int32).  It has no effect on our local complex
    # integers.
    t = numpy.dtype(TYPE_TABLE[wave_info['type']])
    assert bin_info['wave_data_size'] == wave_info['npnts'] * t.itemsize, \
        ('%d, %d, %d, %s' % (bin_info['wave_data_size'], wave_info['npnts'], t.itemsize, t))

    if bin_info['version'] == 5:
        shape = [n for n in wave_info['nDim'] if n > 0]
    else:
        shape = (wave_info['npnts'],)
        
    data_b = buffer(buffer(bin_info['tail_data']) + 
                    f.read(bin_info['wave_data_size']-bin_info['tail_size']))
    data = Wave (shape=shape,
                 dtype=t.newbyteorder(bin_info['byte_order']),
                 buffer=data_b,
                 order='F')
    return data



def wave_read_info(f, wave_info, bin_info):
    '''
    Reads the post-data info of a wave and returns an updated bin_info dictionary.
    '''

    version = bin_info['version']

    if version == 1:
        pass  # No post-data information
    elif version == 2:
        # Post-data info:
        #   * 16 bytes of padding
        #   * Optional wave note data
        pad_b = buffer(f.read(16))  # skip the padding
        log.debug("pad: %s, size: %d" % (pad_b, len(pad_b)))
        assert max(pad_b) == 0, pad_b
        bin_info['note'] = str(f.read(bin_info['noteSize'])).strip()
    elif version == 3:
        # Post-data info:
        #   * 16 bytes of padding
        #   * Optional wave note data
        #   * Optional wave dependency formula
        """Excerpted from TN003:
        
        A wave has a dependency formula if it has been bound by a
        statement such as "wave0 := sin(x)". In this example, the
        dependency formula is "sin(x)". The formula is stored with
        no trailing null byte.
        """
        pad_b = buffer(f.read(16))  # skip the padding
        assert max(pad_b) == 0, pad_b
        bin_info['note'] = str(f.read(bin_info['noteSize'])).strip()
        bin_info['formula'] = str(f.read(bin_info['formulaSize'])).strip()
    elif version == 5:
        # Post-data info:
        #   * Optional wave dependency formula
        #   * Optional wave note data
        #   * Optional extended data units data
        #   * Optional extended dimension units data
        #   * Optional dimension label data
        #   * String indices used for text waves only
        """Excerpted from TN003:
        
        dataUnits - Present in versions 1, 2, 3, 5. The dataUnits
              field stores the units for the data represented by the
              wave. It is a C string terminated with a null
              character. This field supports units of 0 to 3 bytes. In
              version 1, 2 and 3 files, longer units can not be
              represented. In version 5 files, longer units can be
              stored using the optional extended data units section of
              the file.

        xUnits - Present in versions 1, 2, 3. The xUnits field
              stores the X units for a wave. It is a C string
              terminated with a null character.  This field supports
              units of 0 to 3 bytes. In version 1, 2 and 3 files,
              longer units can not be represented.

        dimUnits - Present in version 5 only. This field is an
              array of 4 strings, one for each possible wave
              dimension. Each string supports units of 0 to 3
              bytes. Longer units can be stored using the optional
              extended dimension units section of the file.
              """
        bin_info['formula'] = str(f.read(bin_info['formulaSize'])).strip()
        bin_info['note'] = str(f.read(bin_info['noteSize'])).strip()
        bin_info['dataEUnits'] = str(f.read(bin_info['dataEUnitsSize'])).strip()
        bin_info['dimEUnits'] = [
            str(f.read(size)).strip() for size in bin_info['dimEUnitsSize']]
        bin_info['dimLabels'] = []
        for size in bin_info['dimLabelsSize']:
            labels = str(f.read(size)).split(chr(0)) # split null-delimited strings
            bin_info['dimLabels'].append([L for L in labels if len(L) > 0])
        if wave_info['type'] == 0:  # text wave
            bin_info['sIndices'] = f.read(bin_info['sIndicesSize'])    

    return bin_info


def wave_read(filename):
    '''
    Loads an Igor binary wave from the specified file (either a file stream
    or a file name).
    '''
    if hasattr(filename, 'read'):
        f = filename  # filename is actually a stream object
    else:
        f = open(filename, 'rb')
    try:
        
        wave_info, bin_info = wave_read_header (f)
        version = bin_info['version']

        if wave_info['type'] == 0:
            raise NotImplementedError('Text wave')

        data = wave_read_data (f, wave_info, bin_info)
        bin_info_new = wave_read_info (f, wave_info, bin_info)
        bin_info = bin_info_new

        data.info.update (bin_info)
        data.info.update (wave_info)

        # have all the data, now set explicit scaling information
        log.debug ("load: setting scale on %d dimensions" % len(data.shape))
    
        # set the intrinsic scaling information of the wave.
        for mydim in range(0,len(data.shape)):
            igordim = mydim
            log.debug ("load: Setting scaling info for dimension %d here (%d with Igor): %f/%f"
                       % (mydim, igordim, wave_info["sfA"][igordim], wave_info["sfB"][igordim]))
            data.setScale (mydim, wave_info["sfA"][igordim], wave_info["sfB"][igordim])

        return data

    finally:
        if not hasattr(filename, 'read'):
            f.close()


def wave_find (filename='', pack_tree={}):
    '''
    Tries to locate the specified wave. The file name will usually be a mix-up
    of a file system path and an Igor path (inside a PXP), e.g.:
    
        "/home/user/experiment.pxp:path:to:wave"
    
    In that case, we want this function to:
      . look for a file called "experiment.pxp:path:to:wave"
        (just in case it's a legit file)
      . if it doesn't exist, look for a file called experiment.pxp,
        and a path component "path:to:wave" representing a wave,
        return a file stream object to "/home/user/experiment.pxp" and
        seek the stream to the correct position for :path:to:wave.

    The function will return the 'tree' information of the corresponding
    wave (i.e. the wave name, path components inside the packed file,
    offset and size of the wave data, and real path of the packed file).
    '''

    # if a file name is specified, guess the filesystem path / igor path components
    # if no file name specified, just take the first
    if len(filename) == 0:
        pass

    if not os.path.exists(filename):
        real_filename = os.path.join (os.path.dirname(filename), os.path.basename(filename).split(':')[0])
        igor_path = os.path.basename(filename).split(':')[1:]
        print log.debug ("wave_open: Path: filesystem=%s, igor=%s" % (real_filename, igor_path))
    else:
        real_filename = filename
        igor_path = ''

    if len(pack_tree) == 0:
        pack_tree = pack_scan_tree (real_filename)

    tree = pack_tree

    if (len(igor_path) > 0):
        # extract information from the location specified by the igor path
        for idir in igor_path:
            tree = tree[idir]

        if not tree.has_key('offset'):
            log.error ("wave_find: %s (in %s) is not a wave." % (igor_path, real_filename))
            raise IOError ("wave_find: Component '%s' (inside packed file '%s') is not a wave." % (igor_path, real_filename))

        log.debug ("wave_find: Have branch: %s" % tree)
    
        tree['wname'] = igor_path[-1]
        tree['wpath'] = igor_path[:-1]
        tree['file'] = real_filename
        log.debug ("Location components wave '%s': %s" % (igor_path, str(tree)))

    else:
        # if no igor path is specified, then check if the file itself is a wave
        winfo, wbin = wave_read_header (real_filename)
        tree['wname'] = ''.join(winfo.setdefault('bname', os.path.basename(real_filename)))
        tree['wpath'] = real_filename
        tree['offset'] = 0

    return tree


def wave_save(filename):
    raise NotImplementedError


def pack_scan_tree (filename, dbg_folder_prefix=''):
    '''
    Returns the tree structure of an Igor packed file.
    '''

    pack_tree = {}

    if hasattr(filename, 'read'):
        f = filename  # filename is actually a stream object
    else:
        log.debug ("get_pack_tree: Reading %s" % filename)
        f = open(filename, 'rb')
    try:

        num_waves = 0
        
        while True:
            b = buffer(f.read(PackedFileRecordHeader.size))
            if (len(b) < PackedFileRecordHeader.size):
                break

            rec = PackedFileRecordHeader.unpack_dict_from(b)
            rec_type_id = rec['recordType'] & 0x7fff

            if rec_type_id >= 11:
                #log.debug ('skipping block %d (%d bytes)' % 
                #           (rec_type_id, rec['numDataBytes']))
                f.seek (rec['numDataBytes'], 1)

            elif PackedFileRecordType[rec_type_id] == 'DataFolderStart':
                folder_name = str(f.read(rec["numDataBytes"])).split("\0")[0]
                log.debug ("pack_scan_tree: Folder %s/%s (%d bytes)" %
                           (dbg_folder_prefix, folder_name, rec['numDataBytes']))
                # recursively scan the folder
                pack_tree[folder_name] = pack_scan_tree (f, "%s/%s" % (dbg_folder_prefix, folder_name))

            elif PackedFileRecordType[rec_type_id] == 'DataFolderEnd':
                f.seek (rec['numDataBytes'], 1)
                break # exit recursion step

            elif PackedFileRecordType[rec_type_id] == 'Wave':
                wpos1 = f.tell()
                winfo, wbin = wave_read_header (f)
                wpos2 = f.tell()
                wname = ''.join(winfo.setdefault('bname', 'wave%d' % num_waves))
                pack_tree[wname] = {
                    'offset': wpos1, 
                    'size': rec['numDataBytes'],
                    }
                num_waves = num_waves + 1
                log.debug ("pack_scan_tree: Wave %s/%s (%d bytes, offset %d)" %
                           (dbg_folder_prefix, wname, rec['numDataBytes'], wpos1))
                f.seek ((rec['numDataBytes'] - (wpos2-wpos1)), 1)  # skip the res of the wave data

            else:
                # default behavior is to ignore all other fields
                f.seek (rec['numDataBytes'], 1)
        
    finally:
        if not hasattr(filename, 'read'):
            f.close()

    return pack_tree


def pack_unpack (filename, basedir=".", igordir="", packtree=None):
    '''
    Extracts an Igor packed file, specified by 'filename', to the location
    specified by 'basedir'.
    '''
    
    if packtree == None:
        packtree = pack_scan_tree (filename)

    for key in packtree:
        val = packtree[key]
        basedir_new = os.path.join(basedir, key)
        wpath_igor = igordir+":"+key   # complete igor path, inside the PXP
        wname = key                    # igor wave name
        if val.has_key('offset'):
            print "unpacking IBW file %s \t(from %s)" % (basedir_new, wpath_igor)
            if not os.path.isdir (basedir):
                os.makedirs (basedir)
            dst = open(basedir_new+".ibw", "wb")
            src = open (filename, "rb")
            src.seek (val['offset'])
            dst.write (src.read(val['size']))
            src.close()
            dst.close()
        else:
            pack_unpack (filename, basedir_new, wpath_igor, val)



def main_read_ibw(options):
    """
    IBW -> ASCII conversion
    """

    wav = wave_read (options.infile)
    numpy.savetxt(options.outfile, data, fmt='%g', delimiter='\t')

    if options.verbose > 0:
        pprint.pprint(wav)


def main_read_pack(options):
    log.debug ("reading...")
    #pack_tree = pack_scan_tree (options.infile)
    #pprint.pprint (pack_tree)
    #print "found: %s" % wave_find (options.infile+":HO2:gr2")
    #print "found: %s" % wave_find (options.infile, pack_tree)

    pack_unpack (options.infile)





if __name__ == '__main__':
    import optparse
    import sys
    import pprint

    # preparing the logging system
    ch = logging.StreamHandler()
    ch.setLevel (logging.DEBUG)
    log.setLevel (logging.DEBUG)
    log.addHandler (ch)

    # parsing options
    p = optparse.OptionParser(version=__version__)
    p.add_option('-f', '--infile', dest='infile', metavar='FILE',
                 default='-', help='Input IGOR (.ibw, .pxp, .pxt) file.')
    p.add_option('-o', '--outfile', dest='outfile', metavar='FILE',
                 default='-', help='File for ASCII output.')
    options,args = p.parse_args()

    # some defaults...
    if len(args) > 0 and options.infile == None:
        options.infile = args[0]
    if options.infile == '-':
        options.infile = sys.stdin
    if options.outfile == '-':
        options.outfile = sys.stdout

    # testing tasks
    main_read_pack(options)
