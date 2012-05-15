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

# Begin functions from ReadWave.c

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
        raise ValueError('This does not appear to be a valid Igor binary wave file.'
                         ' The version field = %d.\n', version);
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

#
# loads an IBW from specified file, returns a Wave() object
# (which is an overloaded ndarray) containing the Igor wave data.
#
def load(filename):
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
        bin_struct,wave_struct,checkSumSize = version_structs(version, byteOrder)

        b = buffer(b + f.read(bin_struct.size + wave_struct.size - BinHeaderCommon.size))
        c = checksum(b, byteOrder, 0, checkSumSize)
        if c != 0:
            raise ValueError('load: %s: error in checksum - should be 0, is %d.  '
                             'This does not appear to be a valid Igor binary wave file.'
                             % (filename, c))
        bin_info = bin_struct.unpack_dict_from(b)
        wave_info = wave_struct.unpack_dict_from(b, offset=bin_struct.size)
        if wave_info['type'] == 0:
            raise NotImplementedError('Text wave')
        if version in [1,2,3]:
            tail = 16  # 16 = size of wData field in WaveHeader2 structure
            waveDataSize = bin_info['wfmSize'] - wave_struct.size
            # =  bin_info['wfmSize']-16 - (wave_struct.size - tail)
        else:
            assert version == 5, version
            tail = 4  # 4 = size of wData field in WaveHeader5 structure
            waveDataSize = bin_info['wfmSize'] - (wave_struct.size - tail)
        # dtype() wrapping to avoid numpy.generic and
        # getset_descriptor issues with the builtin Numpy types
        # (e.g. int32).  It has no effect on our local complex
        # integers.
        t = numpy.dtype(TYPE_TABLE[wave_info['type']])
        assert waveDataSize == wave_info['npnts'] * t.itemsize, \
            ('%d, %d, %d, %s' % (waveDataSize, wave_info['npnts'], t.itemsize, t))
        tail_data = numpy.array(b[-tail:])

        if version == 5:
            shape = [n for n in wave_info['nDim'] if n > 0]
        else:
            shape = (wave_info['npnts'],)

        data_b = buffer(buffer(tail_data) + f.read(waveDataSize-tail))
        data = Wave (shape=shape,
                     dtype=t.newbyteorder(byteOrder),
                     buffer=data_b,
                     order='F')

        if version == 1:
            pass  # No post-data information
        elif version == 2:
            # Post-data info:
            #   * 16 bytes of padding
            #   * Optional wave note data
            pad_b = buffer(f.read(16))  # skip the padding
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

    finally:
        if not hasattr(filename, 'read'):
            f.close()

    data.info.update (bin_info)
    data.info.update (wave_info)

    log.debug ("load: setting scale on %d dimensions" % len(shape))
    
    # set the intrinsic scaling information of the wave.
    for mydim in range(0,len(shape)):
        igordim = mydim
        log.debug ("load: Setting scaling info for dimension %d here (%d with Igor): %f/%f"
                   % (mydim, igordim, wave_info["sfA"][igordim], wave_info["sfB"][igordim]))
        data.setScale (mydim, wave_info["sfA"][igordim], wave_info["sfB"][igordim])

    return data


def save(filename):
    raise NotImplementedError


if __name__ == '__main__':
    """IBW -> ASCII conversion
    """
    import optparse
    import sys

    p = optparse.OptionParser(version=__version__)

    p.add_option('-f', '--infile', dest='infile', metavar='FILE',
                 default='-', help='Input IGOR Binary Wave (.ibw) file.')
    p.add_option('-o', '--outfile', dest='outfile', metavar='FILE',
                 default='-', help='File for ASCII output.')
    p.add_option('-v', '--verbose', dest='verbose', default=0,
                 action='count', help='Increment verbosity')
    p.add_option('-t', '--test', dest='test', default=False,
                 action='store_true', help='Run internal tests and exit.')

    options,args = p.parse_args()

    if options.test == True:
        import doctest
        num_failures,num_tests = doctest.testmod(verbose=options.verbose)
        sys.exit(min(num_failures, 127))

    if len(args) > 0 and options.infile == None:
        options.infile = args[0]
    if options.infile == '-':
        options.infile = sys.stdin
    if options.outfile == '-':
        options.outfile = sys.stdout

    data,bin_info,wave_info = loadibw(options.infile)
    numpy.savetxt(options.outfile, data, fmt='%g', delimiter='\t')
    if options.verbose > 0:
        import pprint
        pprint.pprint(bin_info)
        pprint.pprint(wave_info)
