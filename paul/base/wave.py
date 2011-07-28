from numpy import ndarray

class AxisInfo:
    offset = 0 # axis offset (typically min or max)
    delta = 0  # axis increment (negative if offset = max, positive otherwise)
    units = '' # axis units

#
# ndarray with some supplementary information:
#   . per-axis scaling and unit information
#   . notes dictionary
#
class Wave(ndarray):
    ax = []                  # axis vector

    def __init__(self, shape):
        ndarray.__init__(shape)
        self.ax = AxisInfo()
        for i in range(0,len(shape)):
            self.ax[i] = AxisInfo()
        info = { 'notes': '' }   # info map. the field 'notes'
                                 # is implicitly created and will initally
                                 # hold IBW notes. any other fields will be
                                 # converted to a 'key=value' string
                                 # representation when saved to IBW.

    #def __init__(self, data):
        #ndarray.__init__(self, ndarray(data).shape)
        #self.copy(data)

    # overwrites the ndarray.reshape() in order to resize the axes vector 
    #def reshape(self, sizes):
    #    self.ax.resize (len(sizes))
    #    ndarray.reshape (self, sizes)
    #
    
    # sets the axis scaling using delta / offset parameters
    def setScale (self, aindex, delta, offset):
        self.ax[aindex].offset = offset
        self.ax[aindex].delta = delta

    # sets the axis scaling using min/max parameters
    # (rather left/right). depending on which limit
    # is the greater, "delta" might also end up negative!
    def setLimits (self, aindex, left, right):
        offset = left
        delta = (right-left) / self.shape[aindex]
