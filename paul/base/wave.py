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
    def min(self, aindex):
        return min(self.ax[aindex].offset,
                   self.ax[aindex].offset+self.ax[aindex].delta*self.shape[aindex])

    # returns the minimum axis value
    def max(self, aindex):
        return max(self.ax[aindex].offset,
                   self.ax[aindex].offset+self.ax[aindex].delta*self.shape[aindex])

    # Returns a tuple with the axis limits,
    # in the format (ax[0].min, ax[0].max, ax[1].min, ax[1].max...)
    # Useful for imshow()'ing a 2D wave :-)
    def lim(self):
        l = ()
        #print self.ax
        for a in range(len(self.ax)):
            l += ((self.min(a), self.max(a)))
        return l

    # () operator for Wave instaces.
    # Takes a real n-tuple representing a fractional 3D-coordinate within
    # the space spanned by the (min,max) values of the axes.
    # Returns the wave value as the specified position, linearly
    # interpolating if necessary.
    def __call__(self, vals):


        ## FIXME: check if this also works for negative deltas!
        
        # val is a n-tuple (number of axes)
        # initialize some variables
        # i1 = ()
        # i2 = ()   # integer intex above ii
        # xx = ()   # desired x-value (should actually be exactly the same as 'val' :-)
        # x1 = ()   # x-value belonging to i1
        # x2 = ()   # x-value belonging to i2

        # desired position index (fractional index)
        ii = ndarray([len(self.shape)])
        i1 = ndarray([len(self.shape)], dtype=int)
        i2 = ndarray([len(self.shape)], dtype=int)
        print ii
        for v, ai in zip(vals, range(len(self.shape))):
            ii[ai] = (self.ax[ai]._x2i(v))
            i1[ai] = floor(ii[ai])
            i2[ai] = i1[ai]+1


        #_i1 = floor(ii)    # integer index below ii (data type is float!)
        #_i2 = i1+1         # integer index below ii0

        # same as i1 and i2, only with integer arguments (need these for i2x() calls
        #for i in range(len(self.shape)):
        #    _i1[i] = int(i1[i])
        #    _i2[i] = int(i2[i])
        
        # axis values
        xx = self.ax[ai].i2x(ii)
        x1 = self.ax[ai].i2x(i1)
        x2 = self.ax[ai].i2x(i2)


        print "Indices:    ", ii, i1, i2
        print "Axis values:", xx, x1, x2

        # Need to calculate approximation for the n-dim derivative.
        # We do this by modifying one coordinate at a time,
        # from low-index (i1) to high-index (i1), and calculating
        # the y value for that corner.
        
        # This is the block that needs to be changed:
        _i1 = ()
        _i2 = ()
        for a, b in zip(i1,i2):
            _i1 += (a,)
            _i2 += (b,)
        y1 = self[_i1]
        y2 = self[_i2]
        print "Y values:  ", y1, y2
        dy = y2 - y1  # this is a single value, should become a vector
                      # of y-values, one for each dimesion direction

        dx = x2 - x1  # this is an n-tuple of dx values
        print "Deltas:     ", dx, dy

        # 0th order approximation: value at ii is the same as the value at i1
        ynew0 = y1

        # 1st order approximation: linear interpolation (for each dimension)
        ynew1 = (xx-x1) * dy/dx
        print "Correction: ", ynew1

        return ynew0+ynew1.sum()
