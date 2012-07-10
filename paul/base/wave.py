import logging
log = logging.getLogger (__name__)

from numpy import ndarray, floor, array
import math

class AxisInfo:
    offset = 0 # axis offset (typically min or max)
    delta = 1  # axis increment (negative if offset = max, positive otherwise)
    units = '' # axis units string
    
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
    def __new__ (subtype, shape, dtype=None, buffer=None, offset=0, strides=[], order='C'):
        obj = ndarray.__new__ (subtype, shape, dtype, buffer, offset, strides, order)
        for i in range(0,len(shape)):
            obj.info.setdefault('axes', []).append(AxisInfo())
        obj.info['name']='wave%x' % id(obj)
        return obj

    def __array_finalize__ (self, obj):
        if obj is None:
            # this means that object was created explicitly:  foo=Wave()
            pass

        # otherwise, if this is a view casting (ndarray.view()) or
        # from a template type (ndarray[x:y:z]), we need to set 'info'
        # and 'ax' to sane default values.
        self.info = getattr (obj, 'info', {})

    # overwrites the ndarray.reshape() in order to resize the axes vector 
    def reshape(self, sizes):
        self.info['axes'].resize (len(sizes))
        ndarray.reshape (self, sizes)
    

    def setScale (self, aindex, delta, offset):
        '''
        Sets the axis scaling using delta / offset parameters
        '''
        self.axInfo(aindex).offset = offset
        self.axInfo(aindex).delta = delta

    
    def setLimits (self, aindex, left, right):
        '''
        Sets the axis scaling using min/max parameters
        (rather left/right). depending on which limit
        is the greater, "delta" might also end up negative!
        '''
        self.axInfo(aindex).offset = left
        self.axInfo(aindex).delta = (right-left) / self.shape[aindex]


    def axInfo(self, aindex=-1):
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
        '''
        if not isinstance(self, ndarray):
            raise NotImplementedError ("Intrinsic axis information only available for Wave or ndarrays")

        if not isinstance(self, Wave):
            axi = [ AxisInfo for i in (0, self.ndim) ]
        else:
            axi = self.info['axes']
            
        if (aindex < 0):
            return axi
        else:
            return axi[aindex]

    # returns the minimum axis value
    def axMin(self, aindex):
        return min(self.axInfo(aindex).offset,
                   self.axInfo(aindex).offset+self.axInfo(aindex).delta*self.shape[aindex])

    # returns the minimum axis value
    def axMax(self, aindex):
        return max(self.axInfo(aindex).offset,
                   self.axInfo(aindex).offset+self.axInfo(aindex).delta*self.shape[aindex])

    # returns the axis offset
    def axOffset (self, aindex):
        return self.axInfo(aindex).offset

    # returns the axis increment
    def axDelta (self, aindex):
        return self.axInfo(aindex).delta

    # returns the axis endpoint (offset + dim*delta), opposite of offset.
    def axEndpoint (self, aindex):
        return self.axInfo(aindex).offset + self.axInfo(aindex).delta*self.shape[aindex]

    # Returns a tuple with the axis limits,
    # in the format (ax[0].min, ax[0].max, ax[1].min, ax[1].max...)
    # Useful for imshow()'ing a 2D wave :-)
    def axLim(self):
        l = ()
        for a in range(len(self.ax)):
            l += ((self.axMin(a), self.axMax(a)))
        return l

    # same as axLim(), only this one is especially for images,
    # to use with Matplotlib's imshow(). This means that:
    #   . it only works with 2 dimensions
    #   . axes are switched (axis 1 represents left-right boundaries,
    #                        axis 0 represents top-bottom boundaries)
    #
    # The imshow() limit tuple format is: (left, right, top, bottom).
    def imgLim (self):
        return (self.axOffset(1), self.axEndpoint(1), self.axEndpoint(0), self.axOffset(0))

    def i2x(self, aindex, pindex):
        '''
        For the specified axis 'aindex', returns the axis value corresponding
        to the specified point index 'pindex'.
        '''
        return self.axInfo(aindex).i2x(pindex)

    def x2i(self, aindex, pt):
        '''
        For the specified axis 'aindex', returns the index corresponding
        to the specified point 'pt' on the axis.
        '''
        return self.axInfo(aindex).x2i(pt)

    # 
    # () operator for Wave instaces.
    # Takes a variable number of arguments (n-tuple) representing a fractional
    # 3D-coordinate within the space spanned by the (min,max) values of the axes.
    # Returns the wave value as the specified position, linearly
    # interpolating if necessary.
    #
    def __call__(self, *vals):
        
        if len(vals) != len(self.shape):
            raise IndexError (("expected %d dimensions, got %d"
                               % (len(self.shape), len(vals))))
        
        # desired position index (fractional index)
        ii = ndarray([len(self.shape)])
        i1 = ndarray([len(self.shape)], dtype=int)
        i2 = ndarray([len(self.shape)], dtype=int)
        for v, ai in zip(vals, range(len(self.shape))):
            ax = self.axInfo(ai)
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
    def set(obj):
        #for 
        return self
