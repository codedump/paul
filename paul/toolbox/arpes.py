#!/usr/bin/python

import numpy as np
import scipy as sp
from pprint import pprint
import paul.base.wave as w
import logging

log = logging.getLogger (__name__)

'''
Module with tools for analysis of Angle Resolved Photoelectron
Spectroscopy (ARPES) data.
'''

def e(mrel=1.0, Ebind=0.0, kpos=0.0, klim=1.0, pts=10):
    '''
    Returns the parabolic dispersion of an electron bound at
    energy *Ebind*, with relative mass *mrel* (1.0 being
    the mass of the free electron), centered at *kpos* momentum.
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

    A Wave object is returned with proper intrinsic scaling. The object is 
    always 2D, regardless of the user input. But if the user input requested
    a 1D wave, then the 2nd dimension will have only 1 entry.
    '''

    # Make sure 'pts' and 'klim' have sane values.
    # We will work with full 2D pts / klim layout.
    if not hasattr(pts, "__iter__"):
        pts = (pts, 1)

    if len(pts) == 1:
        pts = (pts[0], 1)

    if not hasattr(klim, "__iter__"):
        klim = ( (min(-klim, klim), max(-klim, klim)),
                 (min(-klim, klim), max(-klim, klim)) )
    else:
        new_klim = []
        for l in klim:
            if not hasattr(l, "__iter__") or len(l) == 1:
                new_klim.append ( (min(-l, l),
                                   max(-l, l)) )
            else:
                new_klim.append (tuple(l))
        klim = tuple(new_klim)

    # these are the X and Y axes arrays (1D)
    axx = sp.linspace(klim[0][0], klim[0][1], pts[0])
    axy = sp.linspace(klim[1][0], klim[1][1], pts[1])
    kx = axx[:,np.newaxis]
    ky = axy[np.newaxis,:]

    #log.debug ("axx: %s, axy: %s, klim=%s, pts=%s" % (axx, axy, klim, pts))

    me   = 9.10938215e-31        # free electron mass [kg]
    hbar = 1.054571628e-34       # in    [kg m^2 / s]
    eV   = 1.60217646e-19        # conversion factor J -> eV [kg m^2 / s^2]

    wav = w.Wave (shape=pts)
    out = np.zeros(pts)

    out = Ebind + 1.0/eV * (( (kx**2+ky**2) - (kpos**2))*1.0e20) * (hbar**2.0) / (2.0*mrel*me)

    #for i in range(len(wav.shape)):
    #    wav.setLimits (i, klim[i][0], klim[i][1])

    return out


if __name__ == "__main__":

    fmt = logging.Formatter('%(levelname)s: %(funcName)s: %(message)s')
    ch  = logging.StreamHandler()
    ch.setFormatter(fmt)
    ch.setLevel (logging.DEBUG)
    log.addHandler (ch)
    log.setLevel (logging.DEBUG)

    pprint (e(pts=(5, 4), mrel=-2))

    #e(pts=5, klim=1.0)
    #e(pts=5, llim=(-1.0, 0.5))
    #e(pts=(4, 5), klim=((1.0), (1,0)))
