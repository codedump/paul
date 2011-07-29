from numpy import ndarray

class AxisInfo:
    offset = 0 # axis offset (typically min or max)
    delta = 1  # axis increment (negative if offset = max, positive otherwise)
    units = '' # axis units

#
# ndarray with some supplementary information:
#   . per-axis scaling and unit information
#   . notes dictionary
#
class Wave(ndarray):
    axis = []      # axis info vector
    info = {}    # meta information map  'key': [value, 'unit', ...]

    def __init__(self, shape, dtype=None, buffer=None, offset=0, strides=[], order='C'):
        ndarray.__init__(shape, dtype, buffer, offset, strides, order)
        self.axis = []
        self.info = {}
        for i in range(0,len(shape)):
            self.axis.append(AxisInfo())

    # overwrites the ndarray.reshape() in order to resize the axes vector 
    def reshape(self, sizes):
        self.ax.resize (len(sizes))
        ndarray.reshape (self, sizes)
    
    # sets the axis scaling using delta / offset parameters
    def setScale (self, aindex, delta, offset):
        self.axis[aindex].offset = offset
        self.axis[aindex].delta = delta

    # sets the axis scaling using min/max parameters
    # (rather left/right). depending on which limit
    # is the greater, "delta" might also end up negative!
    def setLimits (self, aindex, left, right):
        offset = left
        delta = (right-left) / self.shape[aindex]

    # returns the minimum axis value
    def min(self, aindex):
        return min(self.axis[aindex].offset,
                   self.axis[aindex].offset+self.axis[aindex].delta*self.shape[aindex])

    # returns the minimum axis value
    def max(self, aindex):
        return max(self.axis[aindex].offset,
                   self.axis[aindex].offset+self.axis[aindex].delta*self.shape[aindex])
