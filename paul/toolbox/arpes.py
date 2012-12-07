#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

import numpy as np
import scipy as sp
import scipy.interpolate as spi
from pprint import pprint
import paul.base.wave as w

from itertools import product

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
    --------------------------------------------------------------
    ACHTUNG: This implementation is broken, gives too high gaps
             (by a factor of somewhere around ~2).
             Besides, this implementation is not physically
             sound (although the results, apart from the
             wrong factor, seem plausible).
             A physically correct implementation would
             involve constructing a Hamiltonian matrix
             with the bands on the diagonals and the
             interaction potentials off-diagonal, and
             diagonalizing that matrix in order to gain
             the hybridized bands.
    --------------------------------------------------------------
    
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


def norm_by_noise (data, dim=0, xpos=(None, None), ipos=(None, None), copy=True):
    '''
    Normalizes 1D sub-arrays obtained from a 2D ndarray
    along axis *dim* by the values integrated along
    *dim* in range specified by *pos*.
    
    In ARPES, this is a useful feature for normalizing    
    spectra obtained in a synchrotron facility, which usually
    have significant amount of 2nd-order intensity above
    the Fermi level. Usually, *dim* would be set to represent
    the energy axis, and *pos* should be set to a range
    well above the Fermi level.

    Parameters:
      - *data* is the (multi-dimensional) ndarray containing
        the data. Data poins _will be modified_, so be sure to
        operate on a copy if you wish to retain original data.
        If the data has more than 2 dimensions, the array
        is rotated such that *dim* becomes the last axis, and
        norm_by_noise() iterated on all elements of data,
        recursively until dimension is 2.
    
      - *ipos* is expected to be a tuple containing a
        (from, to) value pair in 
        If *ipos* is 'None', a background auto-detection is
        attempted using gnd_autodetect().
        If *ipos* is a tuple, but either of the *ipos* elements
        is None, the element is replaced by the beginning
        or ending index, respectively.

      - *xpos* same as *ipos*, only positions are specified
        in axis coordinates of the axis specified by *dim*.
        If present, *xpos* takes precedent over *ipos*.

      - *dim* is the axis representing the array to be, i.e.
        the one to be normalized.

      - *copy* if True (the default), then the normalization will
        be performed on a copy of the original array, the original
        data remaining unchanged. Otherwise the original array
        will be modified.

    '''

    # rotate axes such that dim is last axis
    dim2 = len(data.shape)-1
    if dim2 != dim:
        _data = data.swapaxes(dim, dim2)
    else:
        _data = data

    # we'll be working on Waves all along -- this is
    # because we want to retain axis scaling information
    if copy == True:
        data2 = _data.copy(w.Wave)
    else:
        #if isinstance(data, w.Wave):
        #data2 = _data
        #else:
        data2 = _data.view(w.Wave)
        data2.setflags(write=True)

    # for more than 2 dimensions: work recursively through all dimensions
    if len(data2.shape) > 2:
        for d in data2:
            norm_by_noise (d, dim=(dim2-1), ipos=ipos, xpos=xpos, copy=False)
            
    elif len(data2.shape) == 2:
        # translate everything to index coordinates,
        # xpos has precedence over ipos
        index = [data2.dim[1].x2i_rnd(x) if x is not None
                 else i if i is not None
                 else f
                 for x, i, f in zip(xpos, ipos, (0, data2.shape[1]-1))]

        #print "index: ", index, "shape: ", data2.shape
        for d in data2:
            d /= (d[index[0]:index[1]].sum() / (index[1]-index[0]))
            
    else:
        raise ValueError ("Wrong dimension for normalizing along %d: %s" % (dim, data.shape))

    return data2.swapaxes(dim2, dim) if dim2 != dim else data2


def deg2k (*args, **kwargs):
    '''
    Converts a Wave from the natural coordinates of an ARPES
    measurement (energy / degrees) into k-space (inverse space)
    coordinates.

    Following parametes:
      - (unnamed):    The 3D data, either as a Wave or as an ndarray,
                      containing the deg_tilt*deg_detector*E dependent data
                      (i.e. the intensity in dependence of energy, the tilt
                      angle, and the detector angle). If data is a Wave,
                      axis information is extracted from the wave scaling.
                      'e', 'd' and 't' parameters below override internal
                      Wave data. If data is an ndarray, then 'e', 'd' and 't'
                      parameters are mandatory.
      - axes:         Combination of the letters 'e', 'd' and 't'
                      describing the current axes layout
                      in terms of (e)nergy, (d)etector or (t)ilt.
                      Default is 'edt'.
      - energy, e:    Values to go along with the energy axis.
      - detector, d:  Values for the detector axis.
      - tilt, t:      Values for the tilt axis.
      - hv            The photon energy at which the data was measured.
      
    '''

    # Strategy:
    #  . create a new data view
    #  . tilt data such that axes configuration
    #    is (energy, detector, tilt)
    #  . apply corrections (offsets, and possibly increments)
    #  . ...


    # parameter helpers
    _param = lambda k0, k1, d: \
      kwargs[k0] if kwargs.has_key(k0) \
      else (kwargs[k1] if kwargs.has_key(k1) else d)

    # rotate data into 'edt' axis configuration
    axes = kwargs.setdefault('axes', 'edt')
    idata = np.transpose(args[0], (axes.find('e'), axes.find('d'), axes.find('t')))
    
    E     = _param ('energy',   'e', idata.dim[0].range)
    deg_d = _param ('detector', 'd', idata.dim[1].range)
    deg_t = _param ('tilt',     't', idata.dim[2].range)

    C     = 0.51232 * np.sqrt(E)
    dsin  = np.sin(deg_d*np.pi/180.0)
    tsin  = np.sin(deg_t*np.pi/180.0)

    # transforamtion helpers
    _d2k = lambda d, c: c*np.sin(d*np.pi/180.0)
    
    #_d2kx = lambda d, t, hv: 0.51232*np.sqrt(hv)*np.sin(d*np.pi/180.0)
    #_d2ky = lambda d, t, hv: 0.51232*np.sqrt(hv)*np.sin(t*np.pi/180.0)

    
    # input data coordinates (degrees)
    #icoord_d = ((e, d, t) for e, d, t in product (E, deg_d, deg_t))

    # input data coordinates (k-space)
    icoord_k = np.array([[c, d*c, t*c] for c, d, t in product (C, dsin, tsin)])

    # output array, with proper scales
    E_lim = (min(E), max(E))
    kx_lim = _d2k (np.array([deg_d[0], deg_d[-1]]), E_lim[1])
    ky_lim = _d2k (np.array([deg_t[0], deg_t[-1]]), E_lim[1])
    
    ocoord_k = np.array([[c, kx, ky] for c, kx, ky in \
                product (C, np.linspace(start=kx_lim[0], stop=kx_lim[1], num=len(deg_d)),
                            np.linspace(start=ky_lim[0], stop=ky_lim[1], num=len(deg_t))
                        )
            ])
               
    
    _out = spi.griddata(icoord_k, idata.flat, ocoord_k, method='linear')

    print "finished"

    if isinstance (idata, w.Wave):
        odata = _out.view(w.Wave)
        odata.dim[0].lim = E_lim
        odata.dim[1].lim = kx_lim
        odata.dim[2].lim = ky_lim
    else:
        odata = _out
    
    return odata
    


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
