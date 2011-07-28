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
# Based on code from the Hooke project, by Trevor King.
#

import array
import sys
import types
import struct
import numpy


class Field (object):
    """Represent a Structure field.

    See Also
    --------
    Structure

    We are also introducing a new field type character, 'S', which
    is not included in the original Python specs. 'S' is a variable
    length null-terminated string. The 'count' argument for this
    character is the maximum size, including the trailing 0.
    """
    def __init__(self, format, name, default=None, help=None, count=1):
        self.format = format # See the struct documentation
        self.name = name
        self.default = None
        self.help = help
        self.count = count
        self.total_count = numpy.prod(count)

class Structure (struct.Struct):
    """Represent a C structure.

    A convenient wrapper around struct.Struct that uses Fields and
    adds dict-handling methods for transparent name assignment.

    See Also
    --------
    Field

    Examples
    --------

    Represent the C structure::

        struct thing {
          short version;
          long size[3];
        }

    As

    >>> from pprint import pprint
    >>> thing = Structure(name='thing',
    ...     fields=[Field('h', 'version'),
                    Field('l', 'size', count=3)])
    >>> thing.set_byte_order('>')
    >>> b = array.array('b', range(2+4*3))
    >>> d = thing.unpack_dict_from(buffer=b)
    >>> pprint(d)
        {'size': array([ 33752069, 101124105, 168496141]), 'version': 1}
    >>> [hex(x) for x in d['size']]
        ['0x2030405L', '0x6070809L', '0xa0b0c0dL']

    You can even get fancy with multi-dimensional arrays.

    >>> thing = Structure(name='thing',
    ...     fields=[Field('h', 'version'),
                    Field('l', 'size', count=(3,2))])
    >>> thing.set_byte_order('>')
    >>> b = array.array('b', range(2+4*3*2))
    >>> d = thing.unpack_dict_from (buffer=b)
    >>> d['size'].shape (3, 2)
    >>> pprint(d)
    {'size': array([[ 33752069, 101124105],
           [168496141, 235868177],
           [303240213, 370612249]]),
     'version': 1}
    """
    def __init__(self, name, fields, byte_order='='):
        # '=' for native byte order, standard size and alignment
        # See http://docs.python.org/library/struct for details
        self.name = name
        self.fields = fields
        self.set_byte_order(byte_order)

    def __str__(self):
        return self.name

    def set_byte_order(self, byte_order):
        """Allow changing the format byte_order on the fly.
        """
        if (hasattr(self, 'format') and self.format != None
            and self.format.startswith(byte_order)):
            return  # no need to change anything
        format = []
        for field in self.fields:
            format.extend([field.format]*field.total_count)
        struct.Struct.__init__(self, format=byte_order+''.join(format).replace('P', 'L'))

    def _flatten_args(self, args):
        # handle Field.count > 0
        flat_args = []
        for a,f in zip(args, self.fields):
            if f.total_count > 1:
                flat_args.extend(a)
            else:
                flat_args.append(a)
        return flat_args

    def _unflatten_args(self, args):
        # handle Field.count > 0
        unflat_args = []
        i = 0
        for f in self.fields:
            if f.total_count > 1:
                data = numpy.array(args[i:i+f.total_count])
                data = data.reshape(f.count)
                unflat_args.append(data)
            else:
                unflat_args.append(args[i])
            i += f.total_count
        return unflat_args
        
    def pack(self, *args):
        return struct.Struct.pack(self, *self._flatten_args(args))

    def pack_into(self, buffer, offset, *args):
        return struct.Struct.pack_into(self, buffer, offset,
                                       *self._flatten_args(args))

    def _clean_dict(self, dict):
        for f in self.fields:
            if f.name not in dict:
                if f.default != None:
                    dict[f.name] = f.default
                else:
                    raise ValueError('%s field not set for %s'
                                     % f.name, self.__class__.__name__)
        return dict

    def pack_dict(self, dict):
        dict = self._clean_dict(dict)
        return self.pack(*[dict[f.name] for f in self.fields])

    def pack_dict_into(self, buffer, offset, dict={}):
        dict = self._clean_dict(dict)
        return self.pack_into(buffer, offset,
                              *[dict[f.name] for f in self.fields])

    def unpack(self, string):
        return self._unflatten_args(struct.Struct.unpack(self, string))

    def unpack_from(self, buffer, offset=0):
        return self._unflatten_args(
            struct.Struct.unpack_from(self, buffer, offset))

    def unpack_dict(self, string):
        return dict(zip([f.name for f in self.fields],
                        self.unpack(string)))

    def unpack_dict_from(self, buffer, offset=0):
        return dict(zip([f.name for f in self.fields],
                        self.unpack_from(buffer, offset)))
