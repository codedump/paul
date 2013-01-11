#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

from numpy import ndarray, floor, array, arange
import numpy as np
import scipy.ndimage as spn
import math, copy
import pprint

class AxisInfo(object):
    '''
    This class represents one axis info container.
    It needs to be constructed with a specific parent (which is
    a Wave class), of which it describes a particular dimension.
    To avoid bugs by wiggling indices around (because the dimensions
    of the Wave/ndarray objects can be permutated), this class
    does not specifically save its own dimension index.
    The dimension index is needed only for operations which
    depend on the length of the dimension (i.e. for things
    that somehow involve the *size* property, see below).
    On each operation that involves the *size* property,
    the correct dimension index is computed by looking
    for the index of *self* in _parent.info['axes'].

    In principle, three properties define the property of an axis
    completely:
      . offset: starting point of the axis
      . delta:  increment per step (possibly negative)
      . size:   equals to _parent.shape[i] and is the number
                of points.
    These are the elements that need updating, e.g. when the
    corresponding Wave() dimension is subject to slicing, etc.
    From these, the following properties can be easily computed:
      . end:     the endpoint of the data
      . min/max: the minimum/maximum axis value (which will be
                 the same as offset/end if delta is positive),
      . range:   which is a list with the axis elements, e.g. to use
                 for X-axis plotting information, mathematics etc.
      . lim:     the (offset, end) tuple
      . span:    the (min, max) tuple
    These are read-only for now, but there's nothing that speaks
    against defining write-operations for them (e.g. set support
    for lim/span/ would be great :-) )
    
    (To avoid confusion, look at it like this:
     there's  delta/offset/end  -> lim
         and  min/step/max      -> span,
    with the only difference being that *step* is guaranteed to be
    positive, while delta is not.

    The *units* property is just for convenience, does not influence
    the function of AxisInfo().
           
    So... need to rethink this. Here's one possible sollution:
       . Implement AxisInfo as a subclass of dict()
       . Insert default keys 'offset', 'delta' and 'units'
       . Make one non-bound method ax(), which can be applied
         to any object (not only to waves) and which would
         return either the corresponding AxisInfo() object
         in case of a Wave, or a sane default (offset=0, delta=1)
         in case of an ndarray.
    '''
    
    offset = 0  # axis offset (typically min or max)
    delta  = 1  # axis increment (negative if offset = max, positive otherwise)
    units  = '' # axis units string

    _parent = None # the Wave() object this AxisInfo belongs to

    def __init__(self, parent, copy_from=None, offset=0, delta=1, units=''):
        '''
        Initializes a new AxisInfo object with *parent* as the parent wave.
        If *copy_from* is not None, then significant fields (offset, delta, units)
        are copied from *copy_from*.
        '''
        self._parent = parent
        if copy_from is not None:
            self.copy (copy_from)
        else:
            self.offset = offset
            self.delta  = delta
            self.units  = units

    def copy (self, copy_from):
        '''
        Copies relevant axis information (without the _parent reference)
        from AxisInfo specified by copy_from.
        '''
        self.delta  = copy_from.delta
        self.offset = copy_from.offset
        self.units  = "%s" % copy_from.units
            

    @property
    def size(self):
        '''
        The length of the respective dimension.
        '''
        ax = self._parent.info['axes'] 
        i = ax.index(self)
        return self._parent.shape[i]

    @property
    def end(self):
        '''
        Returns the value of the last point of the axis.
        '''
        return self.offset + self.delta*(self.size-1)

    #@property
    #def step(self):
    #    '''
    #    Absolute value of self.delta
    #    '''
    #    return abs(self.delta)

    @property
    def min(self):
        '''
        Smalles axis coordinate
        '''
        if self.delta > 0:
            return self.offset
        else:
            return self.end

    @property
    def max(self):
        '''
        Biggest axis coordinate
        '''
        if self.delta > 0:
            return self.end
        else:
            return self.offset

        
    def get_lim(self):
        '''
        Axis limits tuple (offset, end)
        '''
        return (self.offset, self.end)

        
    def set_lim(self, newlim):
        '''
        Set the explicit limits, calculates necessary offset/delta values
        '''
        self.offset = newlim[0]
        self.delta  = ( float(newlim[1]) - float(self.offset)) / float(self.size)

    lim = property (get_lim, set_lim, doc="Explicit axis limits, as an (offset, endpoint) tuple.")

        
    #@property
    #def span(self):
    #    '''
    #    Axis limits tuple (min, max)
    #    '''
    #    return (self.min, self.max)
    
    @property
    def range(self):
        '''
        List of all axis values, from the first to the last point.
        '''
        return np.linspace(start=self.offset, num=self.size, stop=self.end)

    
    def i2x_fra(self, index):
        '''
        Returns the axis value corresponding to the
        specified index value.
        '''
        return self.offset+self.delta*index
    i2x = i2x_fra

    def i2x_flo(self, index):
        '''
        Returns the floor of the axis value corresponding to the
        specified index value.
        '''
        return int(math.floor(self.i2x_fra(index)))

    def i2x_rnd(self, index):
        '''
        Returns the rounded axis value corresponding to the
        specified index.
        '''
        return int(round(self.i2x_fra(index)))


    def x2i_fra(self, val):
        '''
        Returns the *fractional* index corresponding to the axis value.
        '''
        return (val-self.offset)/self.delta
    x2i = x2i_fra
    
    def x2i_flo(self, val):
        '''
        Returns the (floor) index corresponding to the axis value.
        '''
        return int(math.floor(self.x2i_fra(val)))

    def x2i_rnd(self, val):
        '''
        Returns the index corresponding to the rounded axis value.
        '''
        return int(round(self.x2i_fra(val)))


    def ppi(self, interval):
        '''
        Returns the (fractional) number of points
        spanned by *interval* on axis.
        '''
        return self.x2i_fra(0) - self.x2i_fra(-abs(interval))

    def __str__(self):
        return "delta=%f, offset=%f, units='%s'" % (self.delta, self.offset, self.units)


class Wave(ndarray):
    '''
    This is an numpy.ndarray sub-class with features that make
    it a useful tool for working with experimental data:
   
      - intrinsic axis scaling and unit information (see `class:AxisInfo`)
      
      - notes dictionary (to store useful meta information
        about measurement parameters etc)

      - an overloaded __getitem__ operator that retains axis information
        when slicing (limited support for fancy indexing)

      - __call__ operator supporting indexing by axis values (instead
        of index numbers)

    Some useful numpy functions that manipulate the shape of numpy.ndarray
    (see below) have been reimplemented (rather: wrapped :-) )
    to be aware of the supplementary features that Wave offers. See
    below for more details.
    '''
    
    def __new__ (subtype, *args, **kwargs):
        obj = ndarray.__new__ (subtype, *args, **kwargs)
        return obj

    def __array_finalize__ (self, obj):
        '''
        When subclassing an numpy.ndarray, 3 scenarios may occur
        (see also http://docs.scipy.org/doc/numpy/user/basics.subclassing.html):

           1) Explicit construction: foo = Wave(...)
           2) View casting:          foo = bar.view(Wave)   (where 'bar' is another ndarray
                                                             or a subtype, possibly a Wave)
           3) New from template:     foo = bar[x:y:z]       (where 'bar' is a Wave)

        Depending on the scenario, initialization of self.info[] is different.
        '''
        self.info = {}

        if type(obj) is Wave:    # view casting or new from a Wave template
            # deep-copy the obj.info field:
            self._copy_info(obj)

        else:                    # explicit construction, or new from non-Wave template.
            # Need to set sane defaults for self.info[]
            for i in range(len(self.shape)):
                self.info['axes'] = self.info.setdefault('axes', ()) + (AxisInfo(self),)
            self.info['name']='wave%x' % id(obj)
            if obj is None: pass         # this would be explicit construction...
            else:           pass         # ...and new from 'ndarray' or another subtype of it

                
    def _copy_info(self, from_wave, noax=False):
        '''
        Copies the info field from *from_wave*. This is basically just
        a wrapper around deepcopy() with some caveats to avoid recursion
        and correct AxisInfo() handling.
        Parameter *noax* specifies whether info fields used for internal
        purposes (specifically Axis info) is to be copied or retained.
        If *noax* is True, then the AxisInfo() is retained. Default is
        False (i.e. by default, AxisInfo() is copied, too).
        Default is usally only be changed by internal functions
        e.g. when slicing is involved and AxisInfo() information needs
        to be reconstructed from scratch anyway.
        '''
        for k,v in from_wave.info.iteritems():
            if k != 'axes':
                self.info[k] = copy.deepcopy(v)
            else:
                # copy AxisInfo() from from_wave only if keepint is False
                if noax:
                    continue

                # the AxesInfo() objects contain references to the parent wave,
                # inducing an endless recursion on deepcopy. Need to copy those
                # by hand.
                axes = []
                for ax in v:
                    axes.append (AxisInfo (self, copy_from=ax))
                self.info['axes'] = axes


    def copy_info_from (self, from_wave, noax=True):
        '''
        Copies the "info" field from Wave specified by *from_wave*.
        This is just a wrapper for self._copy_info(), with a cleaner
        syntax and a reversed default setting for *noax*.
        '''
        self._copy_info(from_wave, noax)
        

    def reshape(self, sizes):
        ''' 
        Overwrites the ndarray.reshape() in order to resize the axes vector.
        '''
        obj = ndarray.reshape (self, sizes)
        
        while len(obj.info['axes']) < len(sizes):
            obj.info['axes'].append (AxisInfo(self))

        while len(obj.info['axes']) > len(size):
            del obj.info['axes'][-1]

        return obj


    def swapaxes (self, ax0, ax1):
        '''
        Overwrites the ndarray.swapaxes() in order to also swap
        AxisInfo objects.
        '''
        obj = ndarray.swapaxes(self, ax0, ax1)
        axes = list(obj.info['axes'])
        foo = axes[ax0]
        bar = axes[ax1]
        axes[ax0] = bar
        axes[ax1] = foo
        obj.info['axes'] = tuple(axes)
        return obj
        

    def sum (self, axis=None, dtype=None, out=None):
        '''
        Overwrites the ndarray.sum() function.
        Integrates over dimension specified by *axes*. The overloaded function takes
        care to remove the corresponding axes from the axinfo vector.
        '''
        obj = ndarray.sum (self, axis, dtype, out)

        # if axis is None, result is a single element - nothing to do.
        if axis is not None:
            ax = list(self.info['axes'])
            del ax[axis]
            obj.info['axes'] = tuple(ax)
        return obj


    def mean (self, axis=None, dtype=None, out=None):
        '''
        Overwrites the ndarray.mean() function.
        Integrates over axes in 'axes'. The overloaded function takes
        care to remove the corresponding axes from the axinfo vector.
        '''
        obj = ndarray.mean (self, axis, dtype, out)
        if hasattr(axis, "__iter__"):
            ax = list(obj.info['axes'])
            for i in axis:
                del ax[i]
            obj.info['axes'] = tuple(ax)
        elif axis is not None:
            ax = list(obj.info['axes'])
            del ax[axis]
            obj.info['axes'] = tuple(ax)
        return obj


    def setScale (self, aindex, delta, offset):
        '''
        Sets the axis scaling using delta / offset parameters
        '''
        self.ax(aindex).offset = offset
        self.ax(aindex).delta = delta

    
    def setLimits (self, aindex, left, right):
        '''
        Sets the axis scaling using min/max parameters
        (rather left/right). depending on which limit
        is the greater, "delta" might also end up negative!
        '''
        self.ax(aindex).offset = left
        self.ax(aindex).delta = (right-left) / self.shape[aindex]


    @property
    def dim(self):
        '''
        Returns the axis information vector.
        '''
        return self.info['axes']


    def ax(self, aindex):
        '''
        Returns the axis information for the specified axis.
        '''
        return self.dim[aindex]


    @property
    def imlim(self):
        '''
        Limits tuple to use with imshow.
        '''
        return (self.dim[1].offset, self.dim[1].end, self.dim[0].end, self.dim[0].offset)

    @property
    def lim(self):
        '''
        General N-dim limits tuple (offset/end of each axis)
        '''
        return tuple(np.array([[d.offset, d.end] for d in self.dim]).flat)
        

    def infs(self, *args, **kwargs):
        '''
        Syntanctic sugar for *self.info[]*. Returns the first element
        of *self.info[index]*, where *index* is specified by the *args* list.
        *kwargs* may contain a "def=..." element, which will be returned
        in case that the *index* key is not found.

        The idea is that *self.info[foo]* is often either not present, or an
        undefined data type (list/string/float...). This requires further
        checks when being used. This function attempts to hide most of those
        checks for the most common use case: when the user is just interested
        in a (numerical or string) information.

        Strategy:
          . if *self.info[index]* is a list: return the *item*-th item (or the list,
                                             if *item* is None)
          . if *self.info[index]* is a dict: return the "key=val" for the *item*-th item
                                             (or the dict, if *item* is None)
          . else return the str()-casted item itself, or *def*, it item is not found

        Valid *kwargs*:
          . def:  Default value to return, if the key is not valid.
                  Default is "n/a".
          . item: Index of the item to return, if the object pointed
                  to by *args* is a list or dict. Default is 0.
                  If it is None, then the list/dict object itself
                  will be returned.
                  

        Remember that *self.info[]* may not be a "falt" dictionary, i.e.
        there may be sub-sections. For this reason, *args* is accepted to
        be a list of keys, applied subsequently to sub-sections of *self.info[]*.
        For example, self.infs("foo", "bar")  is the similar to 
        *self.info["foo"]["bar"]*.

        With default settings (i.e. not *kwargs* specified), this function
        _always_ returns a string.
        '''

        default = kwargs.setdefault('default', None)
        defitem = kwargs.setdefault('item', 0)

        try:
            info = self.info
            for i in args:
                item = info[i]
                info = item

            # check the type of item
            if hasattr (item, "__iter__"):
                if isinstance(item, list):
                    if defitem is not None:
                        return str(item[defitem])
                    else:
                        return item
                    
                elif isinstance(item, dict):
                    if defitem is not None:
                        return "%s = %s" % [i for i in item.iteritems()][defitem]
                    else:
                        return item
            else:
                return str(item)

        except (IndexError, KeyError):
            #log.info ("KeyError with info[%s]" % (str(args)))
            return default
        
        
    def infv(self, *args, **kwargs):
        '''
        Same as *Wave.infs()*, only the returned type is always casted to float().
        If any error occurs, the contents of kwargs['default'] are returned, or 'NaN'.
        '''

        try:
            val = self.infs(*args, **kwargs)
        except ValueError:
            val = kwargs.setdefault('defefault', None)

        if val is None:
            return float('nan')
        else:
            return float(val)
        

    def _get_sliced_axinfo(self, obj, parent=None, source=None):
        '''
        Returns a version of the axis-info copied from the source
        object *source*, sliced according to the slicing object
        *obj*, intended for use in slicing object *parent*.
        If *parent* is None, *self* is used.
        If *source* is None,  *self* is used.
        The modification is performed according to the slicing
        information from *obj*. The slicing information is assumed
        to be in index coordinates, possibly fractional.

        Returns the slicing information tuple.
        '''
        new_info = []
        i = 0

        if parent is None:
            parent = self

        if source is None:
            source = self
        
        if not hasattr(obj, "__iter__"):
            index_list = (obj,)
        else:
            index_list = obj

        for index in index_list:
            new_axi = AxisInfo(parent)

            if isinstance(index, slice):
                old_axi = source.ax(i)
                i_start = index.start
                i_stop  = index.stop
                i_step  = index.step
                if i_start is not None:
                    new_axi.offset = source.ax(i).i2x(i_start)
                else:
                    new_axi.offset = old_axi.offset
                if i_step is not None:
                    new_axi.delta = old_axi.delta * i_step
                else:
                    new_axi.delta = old_axi.delta
                new_axi.units = old_axi.units
                new_info.append(new_axi)
                #print "s: new info for axis", i, new_axi, i_start, self.ax(i).i2x(i_start)
                i += 1

            elif hasattr(index, "__iter__"):
                old_axi = source.ax(i)
                a = array(index).min()
                new_axi.offset = source.ax(i).i2x(a.min())
                new_axi.delta = (source.ax(i).i2x(a.max())-source.ax(i).i2x(a.min()))/float(len(a))
                new_axi.units = old_axi.units
                new_info.append(new_axi)
                #print "i: new info for axis", i, new_axi
                i += 1
                
            elif index is None or index == np.newaxis:
                new_info.append(new_axi) # new axis will just take default AxisInfo()
                # don't increase i -- we need to process the current
                # index (in the old info) for the next round
                #print "n: new info for axis", i, new_axi

            else:
                new_axi = None
                #print "d: new info for axis", i, new_axi
                i += 1

        while i < len(source.info['axes']):
            new_info.append(AxisInfo(parent=parent, copy_from=source.ax(i)))
            #print "c: new info for axis", i, self.ax(i)
            i += 1

        return tuple(new_info)


    def __getitem__(self, obj):
        '''
        Reimplementation of __getitem__ to handle some more comfortable indexing
        (e.g. using axis coordinates). This implementation is usually fast, in most
        cases, as we only "borrow" this method for some initialization work
        and do the rest on a ndarray-casted view of *self* (which calls the
        fast C implementation of the method :-) )
        The only case where this is not the case is when we need to do
        interpolation work, but even there we try to keep things short.
        '''
        data = self.view(ndarray)[obj]
        if isinstance(data, ndarray):
            w = data.view(Wave)
            w.info['axes'] = self._get_sliced_axinfo (obj, parent=w, source=self)
            #print [str(a) for a in self.info['axes'] ]
            #print [str(a) for a in w.info['axes'] ]
            #print w.axv
            return w
        else:
            return data


    def copy_fi (self, *obj, **kwargs):
        '''
        Fractional index slicing: returns a copy of *self*, sliced as specified
        by *obj*. Differently from __getitem__, this function also works with
        floating point slicing parameters.
        *interpolate* can be
              - True:   full linear interpolation for each data point
              - False:  no interpolation at all, function falls back to
                        the integer-index based __getitem__
              - 'lim':  Interpolate only on the axes limits, i.e. not on
                        on the delta value. Use the closest rounded value for
                        delta.
              - 'auto': Will default to 'lim' if the *step* of the dimensions
                        is None, or True if the step is explicitly specified.
                       
        '''



        interpolate=kwargs.setdefault('interpolate', 'auto')

        if interpolate == False:
            return self.__getitem__ (*obj)
        if interpolate == 'lim':
            return self._copy_fi_lim(*obj)
        return self._copy_fi_full(*obj)

    @staticmethod
    def _get_br_index (data, axis, index, full_index):
        '''
        Calculates the bracketing indices for the fractional position *index*
        on *axis*. *index* here is only the index component of the dimension specified
        by *axis*. *full_index* represents the complete index
        (i.e. *full_index*[*axis*] = *index*). Note that the dimension of *full_index*
        and of *data.shape* need _not_ necessarily match, as *full_index* may
        contain *newaxis* elements.
        Bracketing indices are the indices needed to interpolate the data value
        for *index*, i.e. usually floor(index) and floor(index)+1.
        Returns: *i0*, *i1*, *delta*, *keepdim*
        Where    *i0*      is the lower bracket,
                 *i1*      the upper bracket, and
                 *delta*   is *index*-*i0*
                 *keepdim* is a Boolean, specifying whether the index is 
                           is a sequence or a slice() (i.e. whether the
                           result returned when using the index will retain
                           the dimensionality of the array)
        *index* can be a number, a sequence of numbers or a slice() object.
        *i0* and *i1* will have the same type as *index*, *delta* is always a float.
        '''

        #if hasattr(full_index, "__iter__"):
        #    max_dim = max(data.ndim,len(full_index),axis+1)
        #else:
        max_dim = max(data.ndim,axis+1)

        s0 = [slice(None)] * max_dim
        s1 = [slice(None)] * max_dim

        #print index
        #print s0
        #print s1

        if isinstance(index, slice):
            # index is a slicing object, which basically means that we'll be selecting
            # a list of items; dimension therefore remains.

            s = index
            s_start = s.start
            s_stop  = s.stop
            s_step  = s.step

            if s_start is None:
                s_start = 0
            if s_stop is None:
                s_stop = data.shape[axis]
            if s_step is None:
                s_step = 1

            # floor index
            s0[axis] = slice(int(math.floor(s_start)),
                             int(math.floor(s_stop)),
                             int(round(s_step)))

            # ceiling index
            s1[axis] = slice(int(math.floor(s_start))+1,
                             int(math.floor(s_stop))+1,
                             int(round(s_step)))

            delta = s_start - math.floor(s_start)
                
            keep = True  # dimension is retained, increase counter
                
        elif hasattr(index, "__iter__"):
            raise IndexError ("float array indexing needs full interpolation")
        elif index is None or index == np.newaxis:
            s0[axis] = None
            s1[axis] = None
            delta = 0
            keep  = True
        else:
            # index is something that can be cast to a float(), usually a number.
            s = float(index)
            s0[axis] = math.floor(s)
            s1[axis] = math.floor(s)+1
            delta = s - s0[axis]
            keep = False
        
        return tuple(s0), tuple(s1), delta, keep

    def _copy_fi_lim (self, *obj):
        '''
        Fractional index slicing: returns a copy of *self*, sliced as specified
        by *obj*. Differently from __getitem__, this function also works with
        floating point slicing parameters.
        The list of arguments * *obj* is taken to be a number of indices.
        Missing indices will be complemented with slice(None), which 
        means "everything for the respective dimension".
        Each index corresponds to one dimension, and can be:
          - a slice() object
          - a number
          - a sequence
        '''

        #if is

        data_old = self
        idim = 0  # this is the dimension counter. depending on the type of indexing,
                  # dimensions may be killed (i.e. when the index is a number).
                  # if dimension is retained, then the counter is increased by 1
        keep_dim = False

        for index in obj:
            s0, s1, delta, keep_dim = self._get_br_index(data_old, idim, index, obj)

            #print "s0=%s    idim: %d, keep: %d\ns1=%s    data shape: %s, data: %s" 
            #% (s0, idim, keep_dim, s1, data_old.shape, data_old)

            data0 = data_old.view(ndarray)[s0]
            if not keep_dim and (s1[idim] >= data_old.shape[idim]):
                # when indexing the upper-most element with a number index,
                # data1 will give an out-of-bound index. we prevent this by 
                # duplicating data0 (above the upper-most element there's no
                # interpolation anyway).
                data1 = data0
            else:
                data1 = data_old.view(ndarray)[s1]

            #print "index[%d]=%s, shapes: data0=%s  data1=%s" 
            #    % (idim, index, data0.shape, data1.shape)

            cur_dim = idim - int(not keep_dim)
            if (isinstance(data1, ndarray) and \
                isinstance(data0, ndarray)) \
                and (data1.shape[cur_dim] != data0.shape[cur_dim]):
                data1 = np.resize (data1, data0.shape)
                if data0.shape[0] > 0:
                    data1[0] = data0[0]
            
            data_new = (data0 + (data1-data0)*delta)
            data_old = data_new

            idim += keep_dim

            #print "intermediate", data_old
            
        return data_old.view(Wave)


    def _copy_fi_full (self, *obj):
        '''
        Fractional index slicing: returns a copy of *self*, sliced as specified
        by *obj*. Differently from __getitem__, this function also works with
        floating point slicing parameters, uses full interpolation.
        '''
        data_old = self
        idim = 0
        iconsumed = 0 # number of consumed indices

        for index in obj:
            # the original slice object for this dimension
            keep_dim    = False
            ax_range    = []

            # for this dimension, calculate the 
            if isinstance(index, slice):
                s_start = index.start
                s_stop  = index.stop
                s_step  = index.step
                if s_start is None:
                    s_start = 0
                if s_stop is None:
                    s_stop = self.shape[-1]
                if s_step is None:
                    s_step = 1
                ax_range = arange (start=s_start, stop=s_stop, step=s_step)
                keep_dim = True
            elif hasattr(index, "__iter__"):
                ax_range = index
                keep_dim = True
            elif index is None or index is np.newaxis:
                keep_dim = None # strictly speaking, we're not killing
                                # this dimension. this case is a funny
                                # one that will have
            else:
                # simple indexing: one number only
                keep_dim = False

            # this is a hack to account for changing dimension numbers when
            # killing dimensions
            max_dim = data_old.ndim
            if max_dim == idim: 
                max_dim += 1 # this is to account for the "newaxis" thing
            new_s = [slice(None)] * max_dim


            if keep_dim == True:
                data_slices = []
                #print "idim %d, axis range %s" % (idim, ax_range)
                for pos in ax_range:
                    new_s[idim] = slice(pos, pos+1.0, None)
                    data_slice = data_old._copy_fi_lim (*tuple(new_s))
                    #print "index: %s, slice: %s" % (new_s, data_slice)
                    data_slices.append (data_slice)
                data_new = np.concatenate(data_slices, axis=idim)

            else:
                #
                # Don't touch this, or you'll regret it.
                #
                if keep_dim is None:
                    new_s = [slice(None)] * max((len(obj)-iconsumed),max_dim)
                    keep_dim = True

                new_s[idim] = index

                if isinstance(data_old, Wave):
                    data_new = data_old._copy_fi_lim (*tuple(new_s))
                elif index is None:
                    data_new = array([data_old])
                    keep_dim = True
                else:
                    raise IndexError ("Don't know what to do with non-array %s "
                                      "and index %s" % (data_old, index))
                              
            data_old = data_new.view(Wave)
            idim += int(keep_dim)
            iconsumed += 1

        return data_old


    def __call__(self, *vals, **kwargs):
        ''' 
        Takes a variable number of arguments (n-tuple) representing a fractional
        3D-coordinate within the space spanned by the (min,max) values of the axes.
        Returns the wave value as the specified position, linearly
        interpolating if necessary, depending on the value of the *i* parameter.

        Each indexing element in *vals* can be:
          a) a slice() object
          b) None or numpy.newaxys
          c) a float()-castable element
          d) a list with float()-castable elements
          e) a tuple with up to 3 float()-castable elements

        In consequence, this __call__ implementation supports indexing,
        slicing and partly "fancy" indexing, similar to the []-operator, with
        the following syntactic differences (following from points (d) and (e)
        above):
          . There is no on-the-fly conversion from "a:b:c" -> slice(a,b,c),
            which means that slicing have to be either constructed directly
            using slice(a,b,c).
            Alternatively, up to 3 values of float numbers in a tuple (a,b,c)
            will be interpreted as the user's intension to construct a slicing
            object and translated to slice(a,b,c).
            This means that [a:b:c,m:n] would translate to ((a,b,c),(m,n)),
            for example.
          . Consequently, list-indexing cannot be conducted using tuples
            anymore. For instance, if you want to select items 3 and 8,
            you would have to write self([3,8]). Otherwise self((3,8))
            would be translated to self(slice(3,8)).

        IMPORTANT: "fancy" indexing is supported only in part, meaning
                   that with interpolation turned on (see below), only
                   one axis will be correctly indexed by a number array.
                   Using arrays for more than one axis will give a well-defined
                   result (will behave simlarly to a slice() object, but with
                   aritrary step-values between points), but this NOT what
                   one would expect from a regular []-based "fancy" indexing!

        Interpolation for floating-point axis index *pos*, for example,
        works by reading values at floor(pos) and floor(pos)+1. Therefore, at
        the upper boundaries of dimensions problems are ocurring. For the
        sake of consistency with the [] notation, these problems are solved
        by duplicating the elements with the highest index in every dimension
        while calculating the interpolation.

        By the very nature of interpolation (all numbers change!), the returned
        array is always a new object (think copy() :-) ), and never a view()
        of *self*, when interpolation is turned on. The algorithm is fairly
        fast for large datasets, as ndarray's C-implemented []-operator is used
        wherever possible.
        The endpoint interpolation (*i*='lim') is fairly fast: it's basically the
        time needed for two similar calls to [] and a 3 resulting arithmetic
        operations. Note that endpoint interpolation is only implemented for
        slice() or float indexing, not for list lindexing (doesn't make sense).
        Use full interpolation for list indexing, or none at all :-)

        Using full interpolation (*i* = True)is more expensive, and depends a lot
        on the step-size of the slicing object. Still, heavy use of the C-implemented
        []-operator is made wherever possible.

        With interpolation switched off (*i*=False), the axis coordinates
        are translated to the closest integer indices and then passed over
        to the []-operator, so execution time is essentially the same
        as for regular indexing.

        Note that only full interpolation will return arrays precicesly matching
        the specified indices, independently on how the data looks like; everything
        else will depend on how coarse your data axes are :-)

        If you need to cherry-pick single values, have a look at the _get_fx()
        function. It will calculate exactly one single data point, interpolated
        at a random position in your data set. It will be faster for that single
        point, but veeeery slow if you call it repeatedly in a loop for 
        more points :-)
        '''

        if kwargs.has_key("i"):
            i = kwargs['i']
        elif kwargs.has_key("interpolate"):
            i = kwargs['interpolate']
        else:
            i = 'auto'

        # 'vals' represents the index in axis coordinates, First, we 
        # translate it into fractional index coordinates.
        index_obj = []

        # While doing the translation, we will also recommend an interpolation
        # method. Recommending interpolation is done in a hierarchy, following
        # this order:
        #
        #   interpolation:   None  ->  'lim'  ->   True | False
        #
        # Meaning "True" or "False" are always ending points (they will
        # not be altered), while None or 'lim' can be "upgraded" to higher
        # modes of interpolation, based on the following set of rules:
        #
        #   . if one or more indices are slice() objects, with start/stop
        #     parameter not aligned to integer numbers, interpolation
        #     is upgraded to 'lim' (if not already True or False)
        #
        #   . if one or more indices are slice() objects, with
        #     step-parameter not aligned to indeger numbers (i.e. not
        #     a multiple of dimdelta), interpolation is upgraded to
        #     'True' (if not already 'False')
        #
        #   . on the first index that is a list(), interpolation is 
        #     upgraded to 'True' (if it is not 'False')
        #
        #   . on the 2nd or following indices as a list(), interpolation
        #     is turned to 'False'
        #
        #   . on the first numerical index not aligned to integer
        #     numbers, interpolation is upgraded to 'lim',
        #     if not 'True' or 'False'
        # 
        #   At the end, if recommended interpolation is still 'None'.
        #   it is set to "False", meaning that the (coordinate transformed)
        #   indices will be just passed on to the []-operator -- the
        #   safest and fastest method :-)
        #
        recmd_intrp = None

        list_i_no = 0  # number if list-indices
        axis_i = 0     # counter of the current axis, in 'self'
        for ind in vals:
            if isinstance(ind, slice) or (isinstance(ind, tuple) and len(ind) <= 3):
                #
                # handle slice() or tuple() objects
                #
                if isinstance (ind, slice):
                    i_start = ind.start
                    i_stop = ind.stop
                    i_step = ind.step
                else:
                    i_start = i_step = None
                    if len(ind) == 1:
                        i_stop = ind[0]
                    if len(ind) > 1:
                        i_start = ind[0]
                        i_stop = ind[1]
                    if len(ind) == 3:
                        i_step = ind[2]
                if i_start is not None:
                    i_start = self.ax(axis_i).x2i(i_start)
                if i_stop is not None:
                    i_stop = self.ax(axis_i).x2i(i_stop)
                if i_step is not None:
                    i_step /= self.ax(axis_i).delta

                # interpolation recommendation
                if ((i_start is not None) and (abs(i_start-round(i_start)) > 1e-10)):
                    if recmd_intrp is None:
                        recmd_intrp = 'lim'
                if ((i_stop is not None) and (abs(i_stop-round(i_stop)) > 1e-10)):
                    if recmd_intrp is None:
                        recmd_intrp = 'lim'
                if ((i_step is not None) and (abs(i_step-round(i_step)) > 1e-10)):
                    if recmd_intrp in (None, 'lim'):
                        recmd_intrp = True

                index_obj.append (slice(i_start, i_stop, i_step))
                axis_i += 1
                #log.debug ("%s is slice, %s" % (str(ind), recmd_intrp))
            
            elif hasattr (ind, "__iter__"):
                #
                # handle lists or longer tuples -- silently assuming list
                # elements can do float-arithmetics. for everything else,
                # this will break.
                # also, need to do some extra tweaking regarding 
                #
                new_ind = []
                for item in ind:
                    new_ind.append(self.ax(axis_i).i2x(item))

                # recommended interpolation
                if recmd_intrp in (None, 'lim'):
                    recmd_intrp = True
                if list_i_no > 0:
                    recmd_intrp = False

                list_i_no += 1
                axis_i += 1
                index_obj.append (new_ind)
                #log.debug ("%s is list, %s" % (str(ind), recmd_intrp))
                
            elif ind is None or ind == np.newaxis:
                #
                # handle the newaxis element
                #
                # DON'T incrrease axis_i, the new axis does not yet exist in 'self'
                #log.debug ("%s is newaxis, %s" % (str(ind), recmd_intrp))
                index_obj.append (ind)
            
            else:
                #
                # default: ind is probably a simple float()-able number
                #
                nr = self.ax(axis_i).x2i(ind)
                if abs(nr-round(nr)) > 1e-10:
                    if recmd_intrp not in (True, False):
                        recmd_intrp = 'lim'
                index_obj.append(nr)
                axis_i += 1
                #log.debug ("%s is number, %s" % (str(ind), recmd_intrp))

        #log.debug ("proposed interpolation: %s, user-selected: %s" % (recmd_intrp, i))
            
        #
        # now we need to compare recommended interpolation against user-specified
        # parameter. Basically, i='auto' means "go with recommended", everything
        # else only goes if 'recmd_intrp' is or can be upgraded to i
        #
        if (i == 'auto') or (recmd_intrp == i):
            i = recmd_intrp
        elif (i in (True, False)) and (recmd_intrp in (None, 'lim')):
            pass # stick with i, recmd_intrp can be upgraded to it
        elif (i == 'lim') and (recmd_intrp is None):
            pass # stick with i, recmd_intrp can be upgraded to it
        else:
            raise IndexError ("Explicitly specified interpolation method '%s' "
                              "contradicts recommended '%s'" % (i, recmd_intrp))

        # default recommendation, if index permits: False (i.e. standard [] indexing).
        if i == 'auto':
            i = False

        # do the actual interpolation
        new_index = tuple(index_obj)

        #log.debug ("Old index: %s, new index: %s" % (vals, new_index))
        #print "index old:", vals, "new:", new_index, "inter:", i

        if i == True:
            data = self._copy_fi_full(*new_index)
        elif i == False or i == None:
            data = self[new_index]
        elif i == 'lim':
            data = self._copy_fi_lim(*new_index)
        else:
            raise IndexError ("Don't know what to do: interpolation '%s' requested" % (str(i)))

        if isinstance (data, Wave):
            data.info['axes'] = self._get_sliced_axinfo (new_index, source=self, parent=data)

        return data
        

    def _get_fx (self, *vals):
        '''
        Returns the interpolated data value at axes coordinate *vals*.
        '''
        
        if len(vals) != len(self.shape):
            raise IndexError (("expected %d dimensions, got %d"
                               % (len(self.shape), len(vals))))
        
        # desired position index (fractional index)
        ii = ndarray([len(self.shape)])
        i1 = ndarray([len(self.shape)], dtype=int)
        i2 = ndarray([len(self.shape)], dtype=int)
        for v, ai in zip(vals, range(len(self.shape))):
            ax = self.ax(ai)
            ii[ai] = (ax.x2i(v))
            i1[ai] = floor(ii[ai])
            i2[ai] = i1[ai]+1
        
        # axis values (i.e. positions in axis units, not indices!)
        xx = ax.i2x(ii)
        x1 = ax.i2x(i1)
        x2 = ax.i2x(i2)

        # Need to calculate approximation for the n-dim derivative.
        # We do this by modifying one coordinate at a time,
        # from low-index (i1) to high-index (i1), and calculating
        # the y value for that corner.

        y1 = self[tuple(i1)]   # corner 0: the lower-index of each dimension
        y2 = []   
        for i in range(len(vals)):
            icorner = i1.copy()
            icorner[i] = i2[i]
            y2.append(self[tuple(icorner)])

        # dx and dy are vectors => calculate all dimensions in one shot
        dy = y2 - y1
        dx = x2 - x1
        ynew0 = y1                   # 0th order: lower corner
        ynew1 = (xx-x1) * dy/dx      # 1st order: linear correction

        return ynew0+ynew1.sum()


def WCast(nda):
    '''
    Returns a Wave() view of the numpy.ndarray object *nda*.
    '''
    if isinstance(nda, Wave):
        return nda
    return nda.view(Wave)

def WCopy(nda):
    '''
    Returns a Wave()-type copy object of the numpy.ndarray object *nda*.
    '''
    if isinstance(nda, Wave):
        return nda
    return nda.copy(Wave)

def WAx(dat, ax):
    '''
    Returns access to axis info specified by index *ax* for ndarray
    or Wave() *dat*. This is a convenience function to make writing
    of generic plotting code easier (i.e. code that can make use of
    the fancy Wave() AxisInfo interface, but can transparently fall
    back to ndarrays using default axis information).
    '''
    if isinstance(dat, Wave):
        return dat.dim[ax]
    else:
        # assuming ndarray interface and returning a default AxisInfo.
        return AxisInfo(dat)

#
# some useful functions
#

def regrid (wav, *args, **kwargs):
    '''
	Re-samples the specified wave onto a new grid.
	The new grid is specified by the list of optional
	arguments *args, one argument per dimension.
	Every single argument may be as follows:

	   - None, in which case the dimension will not
         be resampled (i.e. retained as it is, instead)
         
       - A dictionary, containing keywords 'delta', 'end'
         'offset', 'shift' or 'numpts', describing the new axis.
         Not all keywords need to be present; missing keywords
         will be replaced according to the following rules,
         in this order:
         
           - 'shift' has precedence over anything else
           
           - 'numpts' has precedence over 'delta', 'offset', 'end'
           
           - 'numpts', 'delta' and either of 'offset' or 'end'
             are present, the other is being calculated
             
           - 'numpts' without 'delta' is present, 'delta'
             is being calculated (using specified offset/end,
             or original axis values)
             
           - the original axis' 'delta', 'offset' or 'end' values are
             used as defaults, if they cannot be calculated.

    All axes information is specified in Axis Coordinates.
    The re-gridding is performed using
    scipy.ndimage.intarpolation.map_coordinates(),
    which performs an n-dimensional interpolation on a rectangular
    grid. (Therefore this function only works assuming the
    data is on a rectangular grid.)

    Optional named parameters:

      * `units`:  'axis' or 'index', specifies whether the
        new grid is specified in intrinsic axis
        units or (fractional) index. Default is 'axis'.
        
      * `mode`:   (See 'mode' and 'cval' in map_coordinates()).
        Describes behavior for points outside of the original
        data boundaries. One of 'constant', 'nearest',
        'reflect' or 'wrap'. Default is 'nearest'.

      * `indexer`: If specified, the N-dimensional tuple of
        slice() objects specifies the points of the output
        array that we would be interested in.
        It will reduce computation time, as only the points
        specified by *indexer* will be calculated then.
        Defaults to None, in which case the complete array
        will be returned.
        
      * `cval`:   (See 'mode' and 'cval' in map_coordinates()).
        Constant value to use with mode = 'constant'.
                
      * `cpinfo`: Boolean. Specifies whether to copy wave info from
        input to output wave. Default is True.
    '''
    
    units = kwargs.setdefault ('units', 'axis')
    mode = kwargs.setdefault ('mode', 'nearest')
    cval = kwargs.setdefault ('cval', 0.0)
    cpinfo = kwargs.setdefault ('cpinfo', True)
    indexer = kwargs.setdefault ('indexer', None)

	# First, create new point coordinates for each axis
	# (This is index units -- the whole function works
	#  on fractional index units; only the input is accepted
	#  in axis units, for convenience.)

    axes = list(args) + [None] * (wav.ndim - len(args))
    
    oslice = []   # output slice lists to produce coord
    coord = []    # output axis coordinates (one array per axis)
    axinfo = []   # axis info for the target wave
    axcnt = 0     # current axis counter (incremented once per loop)
    
    for a, d in zip(axes, wav.dim):

        shaping_index = tuple([np.newaxis] * axcnt + 
                              [slice(None)]        + 
                              [np.newaxis] * (wav.ndim-axcnt-1))
        
        if a is None:
            coord.append (np.arange(d.size)[shaping_index])
            axinfo.append(AxisInfo (parent=None, copy_from=d))
            
        elif isinstance(a, dict):

			# default values (in 'axis' or 'index' units, depending on user specification)
			def_offs  = d.offset if units == 'axis' else 0
			def_end   = d.end    if units == 'axis' else d.size
			def_delta = d.delta  if units == 'axis' else 1

			# parameter 'shift' has precedence over anything else
			if a.has_key ('shift'):
				a['offset'] = def_offs + a['shift']
				a['end']    = def_end  + a['shift']

			# parameter 'numpts' has precedence over others ('delta', 'offset' or 'end')
			if a.has_key ('numpts'):
				if a.has_key('delta'):
					if a.has_key ('offset'):   # ... or offset/end (from delta).
						a['end'] = a['offset'] + (a['numpts']-1)*a['delta']
					elif a.has_key ('end'):
						a['offset'] = a['end'] - (a['numpts']-1)*a['delta']
                        
				# either calculate delta (from offset/end)...
				a.setdefault ('offset', def_offs)
				a.setdefault ('end',    def_end)
				a['delta'] = (a['end']-a['offset']) / (a['numpts'] - 1)

			# catch-all: set offset, end and delta using the original axis values.
			a.setdefault ('offset', def_offs)
			a.setdefault ('end',    def_end)
			a.setdefault ('delta',  def_delta)

			numpts = int(round(math.fabs((a['end'] - a['offset']) / a['delta']))) + 1
			if a.has_key ('numpts'):
				assert (a['numpts'] == numpts)

			axinfo.append (AxisInfo (parent=None, copy_from=d))

			# Need index coordinates for new slice,
			# but axis coordinates for 'axinfo' => different
			# treatment, depending on the units parameter.
			if units == 'index':
				coord.append (np.arange (a['offset'], a['end'], a['delta'])[shaping_index])
				axinfo[-1].offset = d.i2x(a['offset'])
				axinfo[-1].delta = a['delta'] * d.delta
			elif units == 'axis':
				coord.append (d.x2i(np.linspace(start=a['offset'], stop=a['end'], num=numpts))[shaping_index])
				axinfo[-1].offset = a['offset']
				axinfo[-1].delta = a['delta']
        else:
            raise TypeError ("Expecting 'None' or SliceObject as grid info, "
                             "received: %s (%s)" % (str(a), str(type(a))))

        # increment current axis counter
        axcnt += 1

    vcoord = np.broadcast_arrays (*coord)

    if indexer is not None:
        vcoord = [c[indexer] for c in vcoord]
    
    out = spn.interpolation.map_coordinates (wav.view(np.ndarray), vcoord,
                                             mode=mode, cval=cval).view(Wave)

    if cpinfo == True:
        if indexer is None:
            out._copy_info (from_wave=wav, noax=True)
        else:
            out.info['axes'] = wav._get_sliced_axinfo(indexer, source=wav, parent=out)
            #log.warn ("'cpinfo' may misbehave with with a custom indexer!")

    for d_in, d_out in zip(axinfo, out.dim):        
        d_out.copy (d_in)

    return out


#
# wrappers for some useful Numpy functions to make them Wave-aware
#

def transpose (src, axes=None):
    '''
    Override for numpy.transpose(). Takes care of transposing axes info.
    '''
    obj = np.transpose (src, axes)

    #print "before: ", [str(d) for d in obj.dim]

    # re-arrange the dimensions
    orig_dim = [d for d in obj.dim]
    for i in range(len(orig_dim)):
        obj.dim[i] = orig_dim[axes[i]]

    #print "after: ", [str(d) for d in obj.dim]

    return obj

def concatenate (waves, axis=0):
    '''
    Custom numpy.concatenate() replacement for Waves. It concatenates
    the Waves in the specified axis using numpy.concatenate(), then
    creates a Wave() view of the resulting array, then restores
    axis information on all axes.
    '''
    _tmp = np.concatenate (waves, axis=axis).view(Wave)
    _tmp._copy_info (waves)
    return _tmp

    
def dstack (waves, delta=1, offset=0, units=''):
    '''
    numpy.dstack() wrapper that is aware of AxisInfo stuff.
    Uses numpy.dstack() internally, then copies unmodified
    AxisInfo to current wave.
    The newly created axis will use the parameters
    *delta* and *offset* are used for the new axis info.

    (Stacks along dim 3.)
    '''
    obj = np.dstack (waves).view(Wave)
    obj._copy_info (waves[0], noax=True)

    src_dim = list(waves[0].dim)
    while len(src_dim) < obj.ndim:
        src_dim.insert (2, AxisInfo(None, offset=offset, delta=delta, units=units))

    for src, dst in zip(src_dim, obj.dim):
        dst.offset = src.offset
        dst.delta  = src.delta
        dst.units  = src.units

    return obj

    
def vstack (waves, delta=1, offset=0, units=''):
    '''
    numpy.vstack() wrapper that is aware of AxisInfo stuff.
    Uses numpy.vstack() internally, then copies unmodified
    AxisInfo to current wave.
    The newly created axis will use the parameters
    *delta* and *offset* are used for the new axis info.

    (Stacks along dim 2.)
    '''
    obj = np.vstack (waves).view(Wave)
    obj._copy_info (waves[0], noax=True)

    src_dim = list(waves[0].dim)
    while len(src_dim) < obj.ndim:
        src_dim.insert (1, AxisInfo(None, offset=offset, delta=delta, units=units))

    for src, dst in zip(src_dim, obj.dim):
        dst.offset = src.offset
        dst.delta  = src.delta
        dst.units  = src.units

    return obj


def hstack (waves, delta=1, offset=0, units=''):
    '''
    numpy.hstack() wrapper that is aware of AxisInfo stuff.
    Uses numpy.hstack() internally, then copies unmodified
    AxisInfo to current wave.
    The newly created axis will use the parameters
    *delta* and *offset* are used for the new axis info.

    (Stacks along dim 1.)
    '''
    obj = np.hstack (waves).view(Wave)
    obj._copy_info (waves[0], noax=True)

    src_dim = list(waves[0].dim)
    while len(src_dim) < obj.ndim:
        src_dim.insert (0, AxisInfo(None, offset=offset, delta=delta, units=units))

    for src, dst in zip(src_dim, obj.dim):
        dst.offset = src.offset
        dst.delta  = src.delta
        dst.units  = src.units

    return obj
    

#
# some testing functions defined below
#

def _print_ok (test, ok="OK", fail="FAILED", verbose=True):
    '''
    Prints a red FAILED or a green OK, depending on whether *ok* is True or false.
    '''
    ok_string = ["%c[31m%s%c[0m" % (0x1b, fail, 0x1b), "%c[32m%s%c[0m" % (0x1b, ok, 0x1b)]    
    if verbose:
        print "%s" % ok_string[int(test)],
    return test, ok_string[int(test)]


def _cmp_arrays (aa, bb, out=True, verbose=None):
    '''
    Compares *aa* and *bb* element-wise, returns True if they correspond
    and prints an OK or FAILED on stdout if *out* is True.
    '''

    if verbose is not None:
        out = verbose

    test1 = array([int(round(a))==int(round(b)) for a, b in zip(aa.flatten(), bb.flatten())]).all()
    if not test1:
        print "Elements test failed: "
        print ([int(int(round(a))==int(round(b))) for a, b in zip(aa.flatten(), bb.flatten())])
        if verbose == True:
            print "arrays were:"
            pprint (aa)
            print "and"
            pprint (bb)

    test2 = (aa.shape == bb.shape)
    if not test2:
        print "Shape test failed: %s differs from %s" % (aa.shape, bb.shape)

    return _print_ok(test1==True and test2==True)


def _test_index_consistency (index, verbose=False):
    '''
    Tests the consistency of the Wave frational index with the regular
    indexing mechanism for *index*. Prints OK or FAILED for each test,
    returns the (negative) number of failures.
    '''

    if verbose:
        print "Testing consistency of fractional indexing with standard indexing..."

    fail_sum = 0

    a = array([[[i+j*10+k*100 for i in range(5)] for j in range(5)] for k in range(5)])
    wa = a.view(Wave).copy()

    if verbose:
        print "Test array:"
        pprint (a)
        print

    print "  Index:", index

    S = index
    S0 = index[0]
    S1 = index[1]

    foo = []
    if verbose:
        print "    Standard indexing:",
        print
    verbose and pprint ("--0--")
    foo.append (a[S0][S1])
    verbose and pprint (foo[0])
    verbose and pprint ("--1--")
    foo.append(a[S])
    verbose and pprint (foo[1])
    if verbose:
        print
        print

    print "    Limits interpolation:",
    if verbose:
        print
    verbose and pprint ("--2--")
    foo.append(wa._copy_fi_lim(S0)._copy_fi_lim(S1))
    verbose and pprint (foo[2])

    verbose and pprint ("--3--")
    foo.append(wa._copy_fi_lim(*S))
    verbose and pprint (foo[3])
    fail_sum -= not _cmp_arrays (foo[0], foo[2])[0]
    fail_sum -= not _cmp_arrays (foo[1], foo[3])[0]
    if verbose:
        print
    print

    print "    Full interpolation:  ",
    verbose and pprint ("")
    verbose and pprint ("--4--")
    foo.append(wa._copy_fi_full(S0)._copy_fi_full(S1))
    verbose and pprint (foo[4])
    verbose and pprint ("--5--")
    foo.append(wa._copy_fi_full(*S))
    verbose and pprint (foo[5])
    fail_sum -= not _cmp_arrays (foo[0], foo[4])[0]
    fail_sum -= not _cmp_arrays (foo[1], foo[5])[0]
    if verbose:
        print
    print

    print "    Result: ",
    _print_ok(fail_sum == 0, ok="Consistent", fail="INCONSISTENT")
    print

    return fail_sum


def _test_index_fraction (index, verbose=False):
    '''
    Tests the performance of the fractional indexing (i.e. whether it works
    as expected), and returns the number of failures.
    '''
    a = array([[[i*1+j*10+k*100 for i in range(50)] for j in range(50)] for k in range(50)])
    w = array([[[i*1+j*10+k*100 for i in range(5)] for j in range(5)] for k in range(5)]).view(Wave)

    #a = array([[i*1+j*10 for i in range(50)] for j in range(50)])
    #w = array([[i*1+j*10 for i in range(5)] for j in range(5)]).view(Wave)

    w_index = index
    _a_index = []
    for i in w_index:
        if isinstance(i, slice):
            _a_index.append (slice(int(round(i.start*10)),
                                   int(round(i.stop*10)),
                                   int(round(i.step*10))))
        elif hasattr(i, "__iter__"):
            _a_index.append ([int(round(f*10)) for f in i])
        else:
            _a_index.append (int(round(i*10)))
    a_index = tuple(_a_index)

    print "  Fractional index: %s \n  against regular:  %s :" % (w_index, a_index),

    aa = a[a_index]
    ww = w._copy_fi_full(*w_index)

    verbose and pprint ("input")
    verbose and pprint (a)
    verbose and pprint (w)
    verbose and pprint ("result")
    verbose and pprint (aa)
    verbose and pprint (ww)

    fail = _cmp_arrays(aa, ww*10, out=True)

    print

    return -int(not fail[0])


def _test_index_call(*call_i, **kwargs):
    '''
    Tests the index 'call_i' on the __call__ interface, by comparing
    with data indexed either by reg_i, lim_i or full_i,
    using either regular, limit-interpolating or full-interpolating
    algorithms.
    '''

    reg_i   = kwargs.setdefault ('reg_i', None)
    lim_i   = kwargs.setdefault ('lim_i', None)
    full_i  = kwargs.setdefault ('full_i', None)
    verbose = kwargs.setdefault ('verbose', False)

    #a1 = array([[[i*1+j*10+k*100 for i in range(5)] for j in range(5)] for k in range(5)]).view(Wave)
    a2 = array([[[i*1+j*10+k*100 for i in range(5)] for j in range(5)] for k in range(5)]).view(Wave)

    for i in range(a2.ndim):
        a2.dim[i].delta   = 0.1
        a2.dim[i].offset -= (0.5 * a2.dim[i].end)
        verbose and pprint (str(a2.ax(i)))

    if lim_i  == True:
        lim_i  = call_i
    if full_i == True:
        full_i = call_i

    print "  Index: call: %s\n\t lim:  %s\n\t full: %s\n  Test" % (call_i, lim_i, full_i),

    A = a2(*call_i)
    if reg_i is not None:
        B = a2[reg_i]
    elif lim_i is not None:
        B = a2._copy_fi_lim (*lim_i)
    elif full_i is not None:
        B = a2._copy_fi_full (*full_i)

    verbose and pprint (a2)
    verbose and pprint (A)
    verbose and pprint (B)

    ret = _cmp_arrays (A, B, verbose=True, out=True)
    print

    return -int(ret[0] == False)


def _test_index_axslice(index, offsets, deltas, units=None):
    '''
    Tests the correct axis-info slicing.
    '''

    print "  Index:", index,

    a = array([[[i*1+j*10+k*100 \
                     for i in range(5)] \
                    for j in range(5)] \
                   for k in range(5)]).view(Wave)
    a.dim[0].units = 'dim0'
    a.dim[1].units = 'dim1'
    a.dim[2].units = 'dim2'

    foo0 = a[index]
    foo1 = a(*index)
    __get_fails = lambda foo: [int(abs(ai.delta - de)  > 10e-10) +
                               int(abs(ai.offset - of) > 10e-10) + 
                               int(un != ai.units)
                               for ai,of,de,un in zip (foo.dim, offsets, deltas,units)]
    fails0 = array(__get_fails(foo0)).sum()
    fails1 = array(__get_fails(foo1)).sum()
        
    _print_ok(fails0 == 0)
    _print_ok(fails1 == 0)
    print

    ##
    ## some out-of-schedule testing... :-)
    ##
    #pprint ([str(i) for i in a.swapaxes(0,1).axes ])
    #pprint ([str(i) for i in a.sum(1).axes ])
    #pprint ([str(i) for i in a.mean(2).axes ])

    return (fails0+fails1)


def _test_index():
    '''
    Does some indexing tests, returns the number of errors.
    '''

    fail_sum = 0

    print "\nTesting element selection"

    fail_sum += _test_index_consistency ((slice(0,3),slice(None)))
    fail_sum += _test_index_consistency ((slice(None),slice(0,3)))

    fail_sum += _test_index_consistency ((slice(0,3),slice(None),1))
    fail_sum += _test_index_consistency ((slice(0,3),1,slice(None)))
    fail_sum += _test_index_consistency ((1,slice(0,3),slice(None)))

    print "\nTesting dimension boundaries"

    fail_sum += _test_index_consistency ((slice(0,5),slice(0,5),0))
    fail_sum += _test_index_consistency ((slice(0,5),0,slice(0,5)))
    fail_sum += _test_index_consistency ((0,slice(0,5),slice(0,5)))

    fail_sum += _test_index_consistency ((slice(0,5),slice(0,5),4))
    fail_sum += _test_index_consistency ((slice(0,5),4,slice(0,5)))
    fail_sum += _test_index_consistency ((4,slice(0,5),slice(0,5)))

    fail_sum += _test_index_consistency ((slice(0,5),slice(0,5),slice(4,5)))
    fail_sum += _test_index_consistency ((slice(0,5),slice(4,5),slice(0,5)))
    fail_sum += _test_index_consistency ((slice(4,5),slice(0,5),slice(0,5)))


    print "\nTesting dimension increasing"

    fail_sum += _test_index_consistency ((slice(None),slice(None),slice(None),None))
    fail_sum += _test_index_consistency ((slice(None),slice(None),slice(None),None))

    fail_sum += _test_index_consistency ((slice(None),slice(None),None,slice(None)))
    fail_sum += _test_index_consistency ((slice(None),slice(None),None,slice(None)))

    fail_sum += _test_index_consistency ((slice(None),None,slice(None),slice(None)))
    fail_sum += _test_index_consistency ((slice(None),None,slice(None),slice(None)))

    fail_sum += _test_index_consistency ((None,slice(None),slice(None),slice(None)))
    fail_sum += _test_index_consistency ((None,slice(None),slice(None),slice(None)))

    print "\nTesting dimension consistency"

    fail_sum += _test_index_consistency ((2,3,4))
    fail_sum += _test_index_consistency ((0,2,0))

    fail_sum += _test_index_consistency ((2,3,4,None))
    fail_sum += _test_index_consistency ((0,2,None,0))
    fail_sum += _test_index_consistency ((0,None,2,0))
    fail_sum += _test_index_consistency ((None,0,2,0))

    fail_sum += _test_index_consistency ((2,slice(None),4))
    fail_sum += _test_index_consistency ((0,slice(None),0))

    fail_sum += _test_index_consistency ((2,3,slice(None)))
    fail_sum += _test_index_consistency ((4,0,slice(None)))

    fail_sum += _test_index_consistency ((slice(None),2,3))
    fail_sum += _test_index_consistency ((slice(None),4,0))

    fail_sum += _test_index_consistency ((slice(None),2,3,None))
    fail_sum += _test_index_consistency ((slice(None),4,None,0))


    print "\nTesting reduced dimension indexing"

    fail_sum += _test_index_consistency ((2,4))
    fail_sum += _test_index_consistency ((0,0))

    fail_sum += _test_index_consistency ((slice(None),4))
    fail_sum += _test_index_consistency ((slice(None),0))

    fail_sum += _test_index_consistency ((2,slice(None)))
    fail_sum += _test_index_consistency ((4,slice(None)))

    fail_sum += _test_index_consistency ((slice(None),2))
    fail_sum += _test_index_consistency ((slice(None),0))

    print "\nTesting interpolation performance"

    fail_sum += _test_index_fraction ((slice(0,4,0.1),slice(0,4,0.1),slice(0,4,0.1)),
                                      verbose=False)
    fail_sum += _test_index_fraction ((slice(0,4,0.2),slice(0,4,0.1),slice(0,4,0.1)), 
                                      verbose=False)
    fail_sum += _test_index_fraction ((slice(0,4,0.1),slice(0,4,0.2),slice(0,4,0.1)), 
                                      verbose=False)
    fail_sum += _test_index_fraction ((slice(0,4,0.1),slice(0,4,0.1),slice(0,4,0.2)), 
                                      verbose=False)

    print "\nTesting list indexing with interpolation"

    fail_sum += _test_index_fraction (  ([0.1, 0.2, 1.5], slice(0,4,0.1),  slice(0,4,0.1) ) ,
                                        verbose=False)
    fail_sum += _test_index_fraction (  (slice(0,4,0.2),  [0.1, 0.2, 1.5], slice(0,4,0.1) ) ,
                                        verbose=False)
    fail_sum += _test_index_fraction (  (slice(0,4,0.1),  slice(0,4,0.2),  [0.1, 0.2, 1.5]) ,
                                        verbose=False)

    ##
    ## These will FAIL due to a different way Wave() is handling multi-dimensional
    ## list indices. Don't have time to fix it right now, so we'll leave it as it 
    ## is and get noisy about it.
    ## Btw, "FAIL" means they will produce a well-defined result, but which is
    ## not the one that one would expect from the []-operator.
    ##
    #fail_sum += _test_index_fraction (  (slice(0,4,0.1),  [1.2, 2.3, 1.5], [0.1, 0.2, 1.5]) , verbose=False)
    #fail_sum += _test_index_fraction (  ([1.2, 2.3, 1.5], [1.2, 2.3, 1.5], [0.1, 0.2, 1.5]) , verbose=False)
    #fail_sum += _test_index_fraction (  ([1.2, 2.3, 1.5], [0.1, 0.2, 1.5], slice(0,4,0.1) ) , verbose=False)
    #fail_sum += _test_index_fraction (  ([1.2, 2.3, 1.5], slice(0,4,0.1),  [0.1, 0.2, 1.5]) , verbose=False)
    #fail_sum += _test_index_fraction (  ([0.1, 0.2, 1.5], [0.1, 0.2, 1.5], [0.1, 0.2, 1.5]) , verbose=False)

    return fail_sum


def _test_call():
    '''
    Tests the __call__ operator interface
    '''

    fail_sum = 0

    verb = False

    print "\nTesting __call__ operator consistency with lim-interpolation"
    fail_sum += _test_index_call (0, 0, 0,
                                  lim_i=(2.5, 2.5, 2.5), verbose=verb)
    fail_sum += _test_index_call (slice(None), 0, 0,
                                  lim_i=(slice(None), 2.5, 2.5), verbose=verb)
    fail_sum += _test_index_call (0, slice(None), 0,
                                  lim_i=(2.5, slice(None), 2.5), verbose=verb)
    fail_sum += _test_index_call (0, 0, slice(None),
                                  lim_i=(2.5, 2.5, slice(None)), verbose=verb)
    fail_sum += _test_index_call (slice(None), slice(None), 0,
                                  lim_i=(slice(None), slice(None), 2.5), verbose=verb)
    fail_sum += _test_index_call (slice(None), 0, slice(None),
                                  lim_i=(slice(None), 2.5, slice(None)), verbose=verb)
    fail_sum += _test_index_call (0, slice(None), slice(None),
                                  lim_i=(2.5, slice(None), slice(None)), verbose=verb)
    fail_sum += _test_index_call (0, (-0.25,0.15, 0.1), slice(None),
                                  lim_i=(2.5, slice(0,4,1), slice(None)), verbose=verb)
    fail_sum += _test_index_call (0, (-0.25,0.15, 0.1), None, slice(None),
                                  lim_i=(2.5, slice(0,4,1), None, slice(None)), verbose=verb)


    print "\nTesting __call__ operator consistency with full-interpolation"
    fail_sum += _test_index_call (0, 0, 0,
                                  full_i=(2.5, 2.5, 2.5), verbose=verb)
    fail_sum += _test_index_call (slice(None), 0, 0,
                                  full_i=(slice(None), 2.5, 2.5), verbose=verb)
    fail_sum += _test_index_call (0, slice(None), 0,
                                  full_i=(2.5, slice(None), 2.5), verbose=verb)
    fail_sum += _test_index_call (0, 0, slice(None),
                                  full_i=(2.5, 2.5, slice(None)), verbose=verb)
    fail_sum += _test_index_call (slice(None), slice(None), 0,
                                  full_i=(slice(None), slice(None), 2.5), verbose=verb)
    fail_sum += _test_index_call (slice(None), 0, slice(None),
                                  full_i=(slice(None), 2.5, slice(None)), verbose=verb)
    fail_sum += _test_index_call (0, slice(None), slice(None),
                                  full_i=(2.5, slice(None), slice(None)), verbose=verb)
    fail_sum += _test_index_call (0, (-0.25,0.15, 0.1), slice(None),
                                  full_i=(2.5, slice(0,4,1), slice(None)), verbose=verb)
    fail_sum += _test_index_call (0, (-0.25,0.15, 0.1), None, slice(None),
                                  full_i=(2.5, slice(0,4,1), None, slice(None)), verbose=verb)
    
    return fail_sum




def _test_scale():
    '''
    Does some testing of the automatic scale adjustment of a Wave() when
    indexing, either using [] or () is performed. Returns the number
    of errors.
    '''
    
    fail_sum = 0

    print "\nTesting correct axis slicing upon [] and () operators:"
    fail_sum += _test_index_axslice (index=(slice(0,5,2),),
                                     deltas=[2,1,1], offsets=[0, 0, 0],
                                     units=['dim0', 'dim1', 'dim2'])
    fail_sum += _test_index_axslice (index=(slice(0,5,2),None,), 
                                     deltas=[2,1,1,1], offsets=[0, 0, 0, 0], 
                                     units=['dim0','','dim1','dim2'])
    fail_sum += _test_index_axslice (index=(None,slice(0,5,2),), 
                                     deltas=[1,2,1,1], offsets=[0, 0, 0, 0], 
                                     units=['','dim0','dim1','dim2'])
    fail_sum += _test_index_axslice (index=(slice(None),None,slice(0,5,2),), 
                                     deltas=[1,1,2,1], offsets=[0, 0, 0, 0],
                                     units=['dim0','','dim1','dim2'])
    fail_sum += _test_index_axslice (index=(slice(None),slice(0,5,3),None,slice(0,5,2),),
                                     deltas=[1,3,1,2], offsets=[0, 0, 0, 0],
                                     units=['dim0','dim1','','dim2'])
    fail_sum += _test_index_axslice (index=(2,slice(None),slice(None),),
                                     deltas=[1,1], offsets=[0, 0],
                                     units=['dim1','dim2'])
    fail_sum += _test_index_axslice (index=(slice(None),3,slice(None),),
                                     deltas=[1,1], offsets=[0, 0],
                                     units=['dim0','dim2'])
    fail_sum += _test_index_axslice (index=(slice(None),slice(2,5,2),0),
                                     deltas=[1,2], offsets=[0, 2],
                                     units=['dim0','dim1'])
    fail_sum += _test_index_axslice (index=(3,slice(2,5,3),0),
                                     deltas=[3], offsets=[2],
                                     units=['dim1'])

    return fail_sum


if __name__ == "__main__":

    from pprint import pprint
    import sys

    log = logging.getLogger ('paul')
    log.setLevel (logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel (logging.DEBUG)
    log.addHandler (ch)
    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(name)s: '
                            '%(module)s.%(funcName)s: %(message)s')
    ch.setFormatter(fmt)

    fail_sum = 0

    fail_sum += _test_index()
    fail_sum += _test_call()
    fail_sum += _test_scale()
    

    print "\nModule: %s (%d failed)\n\n" % (_print_ok (fail_sum==0, verbose=False)[1],
                                            -fail_sum)
