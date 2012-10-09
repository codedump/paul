#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

from numpy import ndarray, floor, array, arange
import numpy as np
import math, copy

class AxisInfo:
    '''
    This class represents one axis info container.
    The concept probably needs to be reconsidered for a number of reasons:
       (1) I'd like the axis info to be editable together with the
           wave.info[] block seamlessly with other wave.info[] sections.
           For this, the offset/delta/units parameters from AxisInfo
           should be accessible in a dictionarly-like manner
           (and iterable like a dictionary -- maybe derrive AxisInfo
           from dict() could be an option?)
       (2) Currently, AxisInfo is used fully by the Wave object,
           which itself is ndarray based. Working with Waves
           in APIs designed for ndarray, while retaining ndarray
           compatibility, involves asking the object whether it's
           an instance of Wave or not; if yes, use AxisInfo stuff,
           otherwise fall back to ndarray menchanisms (basically
           point indices). This makes programming cumbersome and
           mostly defeats the purpose of AxisInfo (which is making
           life with intrinsically scaled axes easier).
           Therefore, most of the Wave.ax*() functions need
           to be unbound methods (they're currently bound),
           and deliver best restuls for Wave
           instances, and usable results for ndarrays.
       (3) Beyond the basic information (offset, delta, units)
           there is more information often used, that can be
           created from (offset, delta, units), but for which
           there should exist functions to conveniently create it.
           These functions are now spread across two classes:
           One is the AxisInfo.*() stuff (like x2i(), i2x(), ppi()...)
           The other is Wave.ax*() stuff (axMin(), axMax(), axEnd()...).
           On one hand, this is ugly. It would be nice to have
           a consistent interface: everything axes-related either
           inside, or outside AxisInfo.
           On the other hand, some things (like Wave.axMax()
           and axEnd()) need information from the Wave instance
           which it's applied to, which is not avaliable locally
           to AxisInfo (like the dimension size).
           
    So... need to rethink this. Here's one possible sollution:
       . Implement AxisInfo as a subclass of dict()
       . Insert default keys 'offset', 'delta' and 'units'
       . Add a reference to the parent Wave/ndarray to AxisInfo
       . Move all axis related stuff inside AxisInfo
       . Make one non-bound method ax(), which can be applied
         to any object (not only to waves) and which would
         return either the corresponding AxisInfo() object
         in case of a Wave, or a sane default (offset=0, delta=1)
         in case of an ndarray.
       . (Optionally, keep the Wave.ax*() members for backwards
          compatibility, as wrappers to the corresponding AxisInfo
          members.)
       . On slicing operations, i.e. for ndarray views, we need to
         (could?) reset the .delta parameter. This could work
         on simple slicing (i.e. when selecting every N-th slice),
         but will break on more advanced slicing (i.e. on boolean
         indexing).
           
    '''
    
    offset = 0  # axis offset (typically min or max)
    delta  = 1  # axis increment (negative if offset = max, positive otherwise)
    units  = '' # axis units string
    
    def i2x(self, index):
        '''
        Returns the axis value corresponding to the
        specified index value.
        '''
        return self.offset+self.delta*index


    def _x2i(self, val):
        '''
        Returns the *fractional* index corresponding to the axis value.
        '''
        return (val-self.offset)/self.delta

    
    def x2i(self, val):
        '''
        Returns the (floor) index corresponding to the axis value.
        '''
        return int(math.floor(self._x2i(val)))


    def fx2i(self, val):
        '''
        Returns the floor index corresponding to the axis value.
        '''
        return int(math.floor(self._x2i(val)))


    def rx2i(self, val):
        '''
        Returns the index corresponding to the rounded axis value.
        '''
        return int(round(self._x2i(val)))


    def ppi(self, interval):
        '''
        Returns the (fractional) number of points
        spanned by *interval* on axis.
        '''
        return self._x2i(0) - self._x2i(-abs(interval))

    def __str__(self):
        return "delta=%f, offset=%f, units='%s'" % (self.delta, self.offset, self.units)


#
# ndarray with some supplementary information:
#   . per-axis scaling and unit information
#   . notes dictionary
#
class Wave(ndarray):
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
            self.info = copy.deepcopy(obj.info)
        
            ## insert AxisInfo altering code here (e.g. to adjust
            ## the delfa/offset to the new view parameters)
            # ...

        else:                    # explicit construction, or new from non-Wave template.
            # Need to set sane defaults for self.info[]
            for i in range(len(self.shape)):
                self.info['axes'] = self.info.setdefault('axes', ()) + (AxisInfo(),)
            self.info['name']='wave%x' % id(obj)
            if obj is None: pass         # this would be explicit construction...
            else:           pass         # ...and new from 'ndarray' or another subtype of it



    def reshape(self, sizes):
        ''' 
        Overwrites the ndarray.reshape() in order to resize the axes vector.
        '''
        obj = ndarray.reshape (self, sizes)
        obj.info['axes'].resize (len(sizes))
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
        Integrates over axes in 'axes'. The overloaded function takes
        care to remove the corresponding axes from the axinfo vector.
        '''
        obj = ndarray.sum (self, axis, dtype, out)
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


    def ax(self, aindex=-1):
        '''
        Returns the axis information for the specified axis.
        This is a kind of a special function, as it's supposed
        to work not only for Waves, but for any object that's
        derived from numpy.ndarray. For generic ndarray objects
        which don't have intrinsic axis information, the function
        generates harmless defaults (offset=0, delta=1). This
        way, one can use Wave.ax* calls on any kind of ndarray
        without breaking stuff. The only difference to Wave() is
        that with Wave(), it will have a slightly more meaningful
        behavior :-)
        NOTE: This won't work as a member function, need to take this
              out of the object name space.
        '''
        if isinstance(self, Wave):
            axi = self.info['axes']
        elif isinstance(self, ndarray):
            axi = tuple([ AxisInfo() for i in range(self.ndim) ])
            log.debug ("%s is not a Wave, using AxisInfo defaults." % self)
        else:
            raise NotImplementedError ("Intrinsic axis information only available for Wave or ndarrays")

        if (aindex < 0):
            return axi
        else:
            return axi[aindex] 

    def axInfo(self, aindex=-1): # legacy
        log.debug ("Deprecated.")
        return self.ax(aindex)


    def axMin(self, aindex=0):
        '''
        Returns the minimum axis value.
        '''
        return min(self.axOff(aindex), self.axEnd(aindex))


    def axMax(self, aindex=0):
        '''
        Returns the minimum axis value.
        '''
        return max(self.axOff(aindex), self.axEnd(aindex))


    def axOff (self, aindex=0):
        '''
        Returns the axis offset.
        '''
        return self.ax(aindex).offset

    def axOffset (self, ai=0):  # legacy
        log.debug ("Deprecated.")
        return self.axOff(ai)


    def axDelta (self, aindex=0):
        '''
        Returns the axis increment.
        '''
        return self.ax(aindex).delta


    def axEnd (self, aindex=0):
        '''
        Returns the axis endpoint (offset + dim*delta), opposite of offset.
        '''
        return self.ax(aindex).offset + self.ax(aindex).delta*self.shape[aindex]

    def axEndpoint (self, ai):   # legacy
        log.debug ("Deprecated.")
        return self.axEnd(ai)


    def axLim(self):
        '''
        Returns a tuple with the axis limits,
        in the format (ax[0].min, ax[0].max, ax[1].min, ax[1].max...)
        '''
        #l = ()
        #for a in range(len(self.ax())):
        #    l += ((self.axMin(a), self.axMax(a)))
        #return l
        return tuple([ (self.axMin(i), self.axMax(i)) for i in range(len(self.ax())) ])


    def imLim (self):
        '''
        same as axLim(), only this one is especially for images,
        to use with Matplotlib's imshow(). This means that:
          . it only works with 2 dimensions
          . axes are switched (axis 1 represents left-right boundaries,
                               axis 0 represents top-bottom boundaries)
       
        The imshow() limit tuple format is: (left, right, top, bottom).
        '''
        return (self.axOff(1), self.axEnd(1), self.axEnd(0), self.axOff(0))


    def imgLim (self):  # legacy
        log.debug ("Deprecated.")
        return self.imLim()


    def i2x(self, aindex, pindex):
        '''
        For the specified axis 'aindex', returns the axis value corresponding
        to the specified point index 'pindex'.
        '''
        return self.ax(aindex).i2x(pindex)
    fx = i2x

    def x2i(self, aindex, pt):
        '''
        For the specified axis 'aindex', returns the index corresponding
        to the specified point 'pt' on the axis.
        '''
        return self.ax(aindex).x2i(pt)
    fi = x2i


    def ri2x(self, aindex, pindex):
        '''
        For the specified axis 'aindex', returns the *rounded axis value*, i.e. the
        axis value corresponding
        to the specified point index 'pindex'.
        '''
        return self.ax(aindex).ri2x(pindex)
    rx = ri2x


    def rx2i(self, aindex, pt):
        '''
        For the specified axis 'aindex', returns the *rounden index*,
        i.e. the index corresponding
        to the specified point 'pt' on the axis.
        '''
        return self.ax(aindex).rx2i(pt)
    ri = rx2i

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
#        print "item index: ", obj
        return self.view(ndarray)[obj]


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



    def _copy_fi_lim (self, *obj):
        '''
        Fractional index slicing: returns a copy of *self*, sliced as specified
        by *obj*. Differently from __getitem__, this function also works with
        floating point slicing parameters.
        '''

        data_old = self

        for i in range(len(obj))[::-1]:
            # the original slice object for this dimension
            if isinstance(obj[i], slice):
                s = obj[i]
            else:
                s = slice(obj[i])

            s_start = s.start
            s_stop  = s.stop
            s_step  = s.step

            if s_start is None:
                s_start = 0
            if s_stop is None:
                s_stop = data_old.shape[-1]
            if s_step is None:
                s_step = 1

            # floor index
            s0    = [slice(None)] * data_old.ndim
            s0[i] = slice(int(math.floor(s_start)),
                          int(math.floor(s_stop)),
                          int(round(s_step)))

            # ceiling index
            s1    = [slice(None)] * data_old.ndim
            s1[i] = slice(int(math.floor(s_start))+1,
                          int(math.floor(s_stop))+1,
                          int(round(s_step)))

            print
            print s0

            delta = s_start - math.floor(s_start)
            data0 = data_old.view(ndarray)[s0]
            data1 = data_old.view(ndarray)[s1]

            if data1.shape[i] != data0.shape[i]:
                data1 = np.resize (data1, data0.shape)
                data1[0] = data0[0]
            
            data_new = (data0 + (data1-data0)*delta)
            data_old = data_new

            print "intermediate"
            pprint (data_old)
            
        return data_old.view(Wave)


    def _copy_fi_full (self, *obj):
        '''
        Fractional index slicing: returns a copy of *self*, sliced as specified
        by *obj*. Differently from __getitem__, this function also works with
        floating point slicing parameters, uses full interpolation.
        '''
        data_old = self

        for i in range(len(obj)):
            # the original slice object for this dimension
            if isinstance(obj[i], slice):
                s = obj[i]
            else:
                s = slice(obj[i])

            # full interpolation: need to construct a whole new array, slice
            # by slice, at interpolated N*step positions
            data_slices = []
            for pos in arange (start=s.start, stop=s.stop, step=s.step):
                new_s    = [slice(0,dim_max,None) for dim_max in data_old.shape]
                new_s[i] = slice(pos, pos+1.0, None)
                data_slice = data_old._copy_fi_lim (*tuple(new_s))
                data_slices.append (data_slice)
            data_new = np.concatenate(data_slices, axis=i)
            data_old = data_new.view(Wave)

        return data_old


    def test(self,obj):
        print obj


#    def __call__(self, *vals):
#        ''' 
#         () operator for Wave instaces.
#         Takes a variable number of arguments (n-tuple) representing a fractional
#         3D-coordinate within the space spanned by the (min,max) values of the axes.
#         Returns the wave value as the specified position, linearly
#         interpolating if necessary.
#        '''
#        
#        if len(vals) != len(self.shape):
#            raise IndexError (("expected %d dimensions, got %d"
#                               % (len(self.shape), len(vals))))
#        
#        # desired position index (fractional index)
#        ii = ndarray([len(self.shape)])
#        i1 = ndarray([len(self.shape)], dtype=int)
#        i2 = ndarray([len(self.shape)], dtype=int)
#        for v, ai in zip(vals, range(len(self.shape))):
#            ax = self.ax(ai)
#            ii[ai] = (ax._x2i(v))
#            i1[ai] = floor(ii[ai])
#            i2[ai] = i1[ai]+1
#        
#        # axis values (i.e. positions in axis units, not indices!)
#        xx = ax.i2x(ii)
#        x1 = ax.i2x(i1)
#        x2 = ax.i2x(i2)
#
#        # Need to calculate approximation for the n-dim derivative.
#        # We do this by modifying one coordinate at a time,
#        # from low-index (i1) to high-index (i1), and calculating
#        # the y value for that corner.
#
#        y1 = self[tuple(i1)]   # corner 0: the lower-index of each dimension
#        y2 = []   
#        for i in range(len(vals)):
#            icorner = i1.copy()
#            icorner[i] = i2[i]
#            y2.append(self[tuple(icorner)])
#
#        # dx and dy are vectors => calculate all dimensions in one shot
#        dy = y2 - y1
#        dx = x2 - x1
#        ynew0 = y1                   # 0th order: lower corner
#        ynew1 = (xx-x1) * dy/dx      # 1st order: linear correction
#
#        return ynew0+ynew1.sum()

    # 
    # inverse of __call__: sets the values of the wave by evaluating
    # the specified object (function?) in all dimensions at the points
    # specified by self's axis specifications.
    # 
    #def set(obj):
    #    #for 
    #    return self


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



if __name__ == "__main__":

    from pprint import pprint

    log = logging.getLogger ('paul')
    log.setLevel (logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel (logging.DEBUG)
    log.addHandler (ch)
    fmt = logging.Formatter('%(asctime)s %(levelname)s: %(name)s: %(module)s.%(funcName)s: %(message)s')
    ch.setFormatter(fmt)

    a = array([[i+j*10 for i in range(5)] for j in range(5)])
    wa = a.view(Wave).copy()
    

    pprint (a)

#    s = (slice(0,5,0.6),slice(0,5,0.5))
    s = (slice(0,2,None),slice(1,3,None))

    S = (slice(None),slice(2,4))
    pprint ((a[S[0]])[S[1]])
    print "--"
    pprint ((wa._copy_fi_lim(S[0]))._copy_fi_lim(S[1]))
    print "--"
    pprint (wa._copy_fi_lim(*S))
    import sys
    sys.exit(0)

    foo0 = a[s[0]][slice(None),s[1]]
    foo1 = a[s]
    foo2 = wa._copy_fi_lim(s[0])._copy_fi_lim(None,s[1])
    foo3 = wa._copy_fi_lim(*s)
    foo4 = wa._copy_fi_full(s[0])._copy_fi_full(None,s[1])
    #foo5 = wa._copy_fi_full(*s)

    print
    
    print "Regular 0:    "
    pprint (foo0)
    print "Regular 1:    "
    pprint (foo1)

    print

    print "Fractional 2: "
    pprint (foo2)
    print "Fractional 3: "
    pprint (foo3)

    print

    print "Interpolated 4: "
    pprint (foo4)
    print "Interpolated 5: "
    pprint (foo5)
