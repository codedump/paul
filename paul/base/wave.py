import logging
log = logging.getLogger (__name__)

from numpy import ndarray, floor, array
import math

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
        Returns the index corresponding to the axis value
        '''
        return int(math.floor(self._x2i(val)))


    def ppi(self, interval):
        '''
        Returns the (fractional) number of points
        spanned by *interval* on axis.
        '''
        return self._x2i(0) - self._x2i(-abs(interval))


#
# ndarray with some supplementary information:
#   . per-axis scaling and unit information
#   . notes dictionary
#
class Wave(ndarray):
    def __new__ (subtype, *args, **kwargs):
        obj = ndarray.__new__ (subtype, *args, **kwargs)
        log.debug ('new with class %s' % subtype)
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
            for key, val in obj.info.iteritems():
                self.info[key] = val
        
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



    # overwrites the ndarray.reshape() in order to resize the axes vector 
    def reshape(self, sizes):
        self.info['axes'].resize (len(sizes))
        ndarray.reshape (self, sizes)
    

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


    def axMin(self, aindex):
        '''
        Returns the minimum axis value.
        '''
        return min(self.axOff(aindex), self.axEnd(aindex))


    def axMax(self, aindex):
        '''
        Returns the minimum axis value.
        '''
        return max(self.axOff(aindex), self.axEnd(aindex))


    def axOff (self, aindex):
        '''
        Returns the axis offset.
        '''
        return self.ax(aindex).offset

    def axOffset (self, ai):  # legacy
        log.debug ("Deprecated.")
        return self.axOff(ai)


    def axDelta (self, aindex):
        '''
        Returns the axis increment.
        '''
        return self.ax(aindex).delta


    def axEnd (self, aindex):
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


    def x2i(self, aindex, pt):
        '''
        For the specified axis 'aindex', returns the index corresponding
        to the specified point 'pt' on the axis.
        '''
        return self.ax(aindex).x2i(pt)


    def __call__(self, *vals):
        ''' 
         () operator for Wave instaces.
         Takes a variable number of arguments (n-tuple) representing a fractional
         3D-coordinate within the space spanned by the (min,max) values of the axes.
         Returns the wave value as the specified position, linearly
         interpolating if necessary.
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
            ii[ai] = (ax._x2i(v))
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

    # 
    # inverse of __call__: sets the values of the wave by evaluating
    # the specified object (function?) in all dimensions at the points
    # specified by self's axis specifications.
    # 
    #def set(obj):
    #    #for 
    #    return self
