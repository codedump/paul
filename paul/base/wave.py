from numpy import ndarray, floor, array
import math

class AxisInfo:
    offset = 0 # axis offset (typically min or max)
    delta = 1  # axis increment (negative if offset = max, positive otherwise)
    units = '' # axis units
    
    # returns the axis value corresponding to the
    # specified index value
    def i2x(self, index):
        return self.offset+self.delta*index

    # Returns the *fractional* index corresponding to the axis value.
    def _x2i(self, val):
        return (val-self.offset)/self.delta

    # returns the index corresponding to the axis value
    def x2i(self, val):
        return int(math.floor(self._x2i(val)))


#
# ndarray with some supplementary information:
#   . per-axis scaling and unit information
#   . notes dictionary
#
class Wave(ndarray):
    ax = []      # axis info vector
    info = {}    # meta information map  'key': [value, 'unit', ...]

    def __init__(self, shape, dtype=None, buffer=None, offset=0, strides=[], order='C'):
        ndarray.__init__(shape, dtype, buffer, offset, strides, order)
        #ndarray.__new__(shape, dtype, buffer, offset, strides, order)
        self.ax = []
        self.info = {}
        for i in range(0,len(shape)):
            self.ax.append(AxisInfo())

    # overwrites the ndarray.reshape() in order to resize the axes vector 
    def reshape(self, sizes):
        self.ax.resize (len(sizes))
        ndarray.reshape (self, sizes)
    
    # sets the axis scaling using delta / offset parameters
    def setScale (self, aindex, delta, offset):
        self.ax[aindex].offset = offset
        self.ax[aindex].delta = delta

    # sets the axis scaling using min/max parameters
    # (rather left/right). depending on which limit
    # is the greater, "delta" might also end up negative!
    def setLimits (self, aindex, left, right):
        self.ax[aindex].offset = left
        self.ax[aindex].delta = (right-left) / self.shape[aindex]

    # returns the minimum axis value
    def axMin(self, aindex):
        return min(self.ax[aindex].offset,
                   self.ax[aindex].offset+self.ax[aindex].delta*self.shape[aindex])

    # returns the minimum axis value
    def axMax(self, aindex):
        return max(self.ax[aindex].offset,
                   self.ax[aindex].offset+self.ax[aindex].delta*self.shape[aindex])

    # Returns a tuple with the axis limits,
    # in the format (ax[0].min, ax[0].max, ax[1].min, ax[1].max...)
    # Useful for imshow()'ing a 2D wave :-)
    def axLim(self):
        l = ()
        #print self.ax
        for a in range(len(self.ax)):
            l += ((self.axMin(a), self.axMax(a)))
        return l

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
            ii[ai] = (self.ax[ai]._x2i(v))
            i1[ai] = floor(ii[ai])
            i2[ai] = i1[ai]+1
        
        # axis values (i.e. positions in axis units, not indices!)
        xx = self.ax[ai].i2x(ii)
        x1 = self.ax[ai].i2x(i1)
        x2 = self.ax[ai].i2x(i2)

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
