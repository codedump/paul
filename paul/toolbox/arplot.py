#!/usr/bin/python

import matplotlib.collections as mpl_col
import matplotlib.transforms as mpl_trans

import numpy as np

import itertools

#
# Some Matplotlib plotting helpers, specifically for ARPES
#

def plot_bz2d (ax, bz, repeat=(1, 1), kvec=None, rotation=0,
               labels=None, color=(.2, .7, .7, .7)):
    '''
    Plots the Brollouin Zone (BZ) projection specified by points
    *bz* as a line onto axis *ax*, rotated by *rotation* degrees.
    If *kvec* is not None, it is expected to be a list of (x,y)
    tuples representing the projections of the *kvec* vectors
    by which to repeat the BZ *repeat* times in each direction.
    '''

    kspace = []

    # if *repeat* is non-enumerable, we assume it's
    # something castable to an int, and use that in both directions.
    if not hasattr (repeat, "__len__"):
        repeat = (int(repeat), int(repeat))
    
    if kvec is None:
        kvec = [(0,0), (0,0)]

    comb = [r for r in itertools.product (range(-repeat[0], repeat[0]+1),
                                          range(-repeat[1], repeat[1]+1))]
    for c in comb:
        kspace.append (np.array(kvec[0])*c[0]+np.array(kvec[1])*c[1]+bz)

    lines = mpl_col.LineCollection (kspace, color=color)

    trans = mpl_trans.Affine2D().rotate_deg (rotation) + ax.transData
    lines.set_transform (trans)

    ax.add_collection (lines)

    return lines

