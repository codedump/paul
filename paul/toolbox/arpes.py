#!/usr/bin/python

import numpy as np
import paul.base.wave as w
import logging

log = logging.getLogger (__name__)

'''
Module with tools for analysis of Angle Resolved Photoelectron
Spectroscopy (ARPES) data.
'''

def e(mrel=1.0, Ebind=0.0, klim=((-1.0), (-1.0)), pts=(100, 100)):
    '''
    Returns the parabolic dispersion of an electron bound at
    energy *Ebind*, with relative mass *mrel* (1.0 being
    the mass of the free electron).
    *pts* contains the number of points in k-direction. If it's
    a single number, a 1D object wull be generated. If it's a
    tuple of at least 2 numbers, then a 2D dispersion is returned
    (i.e. a paraboloid).
    The dispersion is returned for k-values specified by *klim*.
    The semification of *klim* depends on the dimensionality
    expected. For 1D results (i.e. 1D pts value):
      1) *klim* is expected to be an iterable with 2 elements
      2) if *klim* is a single number, then (-*klim*, *klim*) is assumed
    For 2D results:
      4) *klim* is supposed to be a 2x2 iterable.
      5) if *klim* is a number, then ((-*klim*, *klim*), (-*klim*, *klim*)) is assumed
      6) if *klim* is ((nr), (nr)), then ((-nr, nr), (-nr, nr)) is assumed
    A Wave object is returned with proper intrinsic scaling.
    '''
    
    out = np.zeros(pts)
    axx = np.zeros(pts)
    axy = np.zeros(pts)

    # make 'pts' an iterable, now that 

    if len(out.shape) == 1:
        if not hasattr(klim, "__iter__") or len(klim) == 1:
            klim = (-klim, klim)
        if not hasattr(pts, "__iter__") or len(pts) == 1:
            pts = (pts,)
        axx = np.arange(start=klim[0], stop=klim[1], step=(klim[1]-klim[0])/pts[0])
    elif len(out.shape) == 2:
        if not hasattr(klim, "__iter__"):
            klim = ((-klim, klim), (-klim, klim))
        else:
            new_klim = []
            for l in klim:
                if not hasattr(l, "__iter__") or len(l) == 1:
                    new_klim.append (-l, l)
                else:
                    new_klim.append (tuple(l))
                klim = new_klim
        
        if not hasattr(pts, "__iter__") or len(pts) == 1:
            pts = (pts,pts)
        axx = np.arange(start=klim[0][0], stop=klim[0][1], step=(klim[0][1]-klim[0][0])/pts[0])
        axy = np.arange(start=klim[1][0], stop=klim[1][1], step=(klim[1][1]-klim[1][0])/pts[1])

    else:
        log.error ("Expecting 1D or 2D request")

    log.debug ("axx: %s, axy: %s, klim=%s, pts=%s" % (axx, axy, klim, pts))
