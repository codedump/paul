#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

import numpy as np
import scipy as sp
from pprint import pprint
import paul.base.wave as w


'''
Module with tools for analysis of Angle Resolved Photoelectron
Spectroscopy (ARPES) data.
'''

def e(mrel=1.0, ebind=0.0, kpos=0.0, klim=1.0, pts=100, out=None):
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

    # Make sure 'pts' and 'klim' have sane values. We will work with full 2D pts / klim layout.

    dim = 2

    # if input is 1D, add a 2nd dimension
    if not hasattr(pts, "__iter__"):
        pts = (pts, 1)
        dim = 1
    if len(pts) == 1:
        pts = (pts[0], 1)
        dim = 1

    if not hasattr(klim, "__iter__"):
        if dim == 2:
            klim = ( (min(-klim, klim), max(-klim, klim)),
                     (min(-klim, klim), max(-klim, klim)) )
        else:
            klim = ( (min(-klim, klim), max(-klim, klim)),
                     (0.0, 0.0))

    else:
        new_klim = []
        if dim == 2:
            for l in klim:
                if not hasattr(l, "__iter__") or len(l) == 1:
                    new_klim.append ( (min(-l, l),
                                       max(-l, l)) )
                else:
                    new_klim.append (tuple(l))
        elif dim == 1:
            new_klim = [klim, (0.0, 0.0)]
        klim = tuple(new_klim)

    # these are the X and Y axes arrays (1D)
    axx = sp.linspace(klim[0][0], klim[0][1], pts[0])
    axy = sp.linspace(klim[1][0], klim[1][1], pts[1])
    kx = axx[:,np.newaxis]
    ky = axy[np.newaxis,:]

    #print ("axx: %s, axy: %s\nklim=%s, pts=%s" % (axx, axy, klim, pts))

    me   = 9.10938215e-31        # free electron mass [kg]
    hbar = 1.054571628e-34       # in    [kg m^2 / s]
    eV   = 1.60217646e-19        # conversion factor J -> eV [kg m^2 / s^2]

    out = ebind + 1.0/eV * (( (kx**2+ky**2) - (kpos**2))*1.0e20) * (hbar**2.0) / (2.0*mrel*me)

    if type(out) is w.Wave:
        wav = out
    else:
        wav = out.view(w.Wave)

    for i in range(len(wav.shape)):
       wav.setLimits (i, klim[i][0], klim[i][1])

    return wav


def hybridize(wlist, V=0.0, count=1):
    '''
    Hybridizes bands from *wlist* using the coupling matrix *V*.
    *V* is a NxN matrix, where N = len(wlist).
    Returns a list hybridized bands, corresponding to:

       hi/hj = 1/2 * (wlist[i] + wlist[j] +/- sqrt((wlist[i]-wlist[j])**2 + 4*abs(v)))

    If *count* is specified, then the procedure will be repeated *count*
    times (default is 1). This is intended for multi-band hybridization,
    where one pass may not be enough.
    '''
    # if only one value is specified, construct a coupling matrix out of it
    if type(V) == float:
        V = np.matrix([[V for i in wlist] for j in wlist]) / 2

    V -= np.diag(np.diag(V))

    V2  = V + V.T - 2*np.diag(np.diag(V))      # make sure matrix is symmetric
    V2 -= np.diag(np.diag(V))                  # remove diagonal elements

    ''' This is what we'll be operating on'''
    hlist = list(wlist)
    
    '''
    Hybridizing bands more than once brings some normalization problems:
    the interaction potential V will be applied to the complete band.
    The formula has most drastic consequences for the crossover-points,
    but it actually affects the _whole_ band. Now, if we repeatedly
    hybridize two bands, then they will get pushed appart, even if they
    don't cross anymore.
    
    To avoid this, we need to normalize the interaction potential by the
    number of times that we're going to hybridize. (Mind the sqrt() -- this
    is because the potential needs to be halvened under the sqrt().)
    '''
    norm = 1.0/np.sqrt(count)

    '''
    Next problem arises when hybridizing more than 2 bands: multi-band
    hybridization is performet in a each-with-every-other kind of loop.
    So, we need to keep the interaction potential even smaller.
    The (missing) factor 1/2 accounts for the fact that our algorithm,
    for symmetry reasons, would hybridize a i->j and j->i.
    There is no sqrt() here, because this normalization factor is
    '''
    norm *= 1.0/(len(wlist)*(len(wlist)-1)/2.0)

    
    V2 *= norm  # The actual normalization

    '''
    Now, V2 is a traceless, symmetric matrix containing coupling
    factors for the bands in 'hlist', correctly normalized. Let's go!
    '''
    for t in range(count):
        for i in range(len(hlist)):
            for j in range(len(hlist)):
                if i == j:
                    continue
                h1 = 0.5 * (hlist[i] + hlist[j]) + np.sqrt( (0.5*(hlist[j]-hlist[i]))**2 + V2[i,j]**2 )
                h2 = 0.5 * (hlist[i] + hlist[j]) - np.sqrt( (0.5*(hlist[j]-hlist[i]))**2 + V2[i,j]**2 )
                hlist[i] = h1
                hlist[j] = h2

    return hlist



if __name__ == "__main__":

    log = logging.getLogger ("paul")

    fmt = logging.Formatter('%(levelname)s: %(funcName)s: %(message)s')
    ch  = logging.StreamHandler()
    ch.setFormatter(fmt)
    ch.setLevel (logging.DEBUG)
    log.addHandler (ch)
    log.setLevel (logging.DEBUG)

    #foo = e(pts=(5, 4), mrel=-2)
    foo = e()
    pprint (foo)
    print foo.info['axes']

    #e(pts=5, klim=1.0)
    #e(pts=5, llim=(-1.0, 0.5))
    #e(pts=(4, 5), klim=((1.0), (1,0)))
