#!/usr/bin/python

from matplotlib.collections import LineCollection
import itertools
import numpy as np

# BZ projection along 001, as a (X,Y) line; repeat by kx*ky.
uru2si2_001 = [
 ( .9004750999, -.6241998107),
 ( .9004750999,  .6241998107),
 ( .6241998107,  .9004750999),
 (-.6241998107,  .9004750999),
 (-.9004750999,  .6241998107),
 (-.9004750999, -.6241998107),
 (-.6241998107, -.9004750999),
 ( .6241998107, -.9004750999),
 ( .9004750999, -.6241998107),
]

# kx, ky and ky vectors in 3D; when using,
# ignore dimensions as required
uru2si2_kvec = [
  (2*1.5246749107, 0),
  (1.5246749107, 1.5246749107)
]

# distance: 1.5246749107

def plot_bz2d (ax, bz, repeat=(1, 1), kvec=None, rotate=0):
    '''
    Plots the Brollouin Zone (BZ) projection specified by points
    *bz* as a line onto axis *ax*, rotated by *rotate* degrees.
    If *kvec* is not None, it is expected to be a list of (x,y)
    tuples representing the projections of the *kvec* vectors
    by which to repeat the BZ *repeat* times in each direction.
    '''

    kspace = [bz]
    
    if kvec is not None:
        comb = [r for r in itertools.product (range(-repeat[0], repeat[0]+1),
                                              range(-repeat[0], repeat[1]+1))]
        
        for c in comb:
            c = np.array(c)
            rep = (kvec[0]*c + kvec[1]*c)
            new_bz = rep + np.array(bz)
            kspace.append (np.array(kvec[0])+bz)
            kspace.append (np.array(kvec[1])+bz)
            #kspace.append(new_bz)
        

    lines = LineCollection (kspace)
    ax.add_collection (lines)

    return lines

def decorate (*args, **kwargs):
    ax = kwargs['can'].axes
    ax.set_aspect (1)
    ax.axhline (0, ls=':', color='black')
    ax.axvline (0, ls=':', color='black')

    ax.set_xlim (-3, 3)
    ax.set_ylim (-3, 3)

    bz = plot_bz2d (ax, uru2si2_001, repeat=(2, 2),
                  kvec=[(v[0], v[1]) for v in uru2si2_kvec[0:2]])
    
