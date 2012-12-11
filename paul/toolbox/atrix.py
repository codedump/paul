#!/usr/bin/python

from paul.base.wave import Wave
import numpy as np
import copy

#
# Array Tricks:
#
# ...is a common namespace for numpy.ndarray
# manimpulations too simple to go anywhere else, but comlpicated
# enough not to want to type them by hand all the time :-)
#

def ncomp (iwave, axis=0, step=1, intg=-1, norm=False):
    '''
    This is a mix of ndarray.compress() and step-wise slicing.
    For an N-dim array, selects every 'step'-th (N-1)-dim sub-array
    and integrates over 'intg' such slices. This function is e.g. useful
    to create waterfall diagrams of 2D arrays, where the waterfall
    graphs (single 1D arrays) are integrated over 'step' slices.
    Returns a new N-dim ndarray (or Wave?), but with dimension
    size in 'axis' direction reduced.
    If 'norm' is True, then each slice will be divided by
    the surface below the graph.
    '''

    # create a new wave, same as input wave, but with a reduced
    # size in dimention 'axis'. retain axis info, alter it
    # to match the new dimension.
    new_shape = list(iwave.shape)
    new_shape[axis] = iwave.shape[axis] / step
    owave = None
        
    # ignore the points that don't align well with 'step'
    max_i = new_shape[axis] * step

    if intg == -1:
        intg = step

    nval = 1.0 / intg
    for s in range(0, intg):
        # marker:  select every N-th point ...but ignore bogus points at the end
        marker = [ ((((step+i-s) % step) == 0) and (i < max_i))*1 for i in range (0, iwave.shape[axis]) ]
        tmp = iwave.compress (marker, axis) * nval
        if owave is None:
            owave = tmp.copy()
        else:
            owave += tmp

    owave.dim[axis].delta *= step
        
    if norm:
        for s in owave.swapaxes(0,axis):    # scaling and normalizing
            s /= s.sum()

    return owave


