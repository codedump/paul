#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

import numpy as np
import scipy as sp
import scipy.interpolate as spi
import scipy.ndimage.interpolation as spni
import scipy.ndimage as spimg
from pprint import pprint
import paul.base.wave as wave
from paul.toolbox.atrix import ncomp

from itertools import product
import pprint
import math

'''
Module with tools for analysis of Angle Resolved Photoelectron
Spectroscopy (ARPES) data.
'''


#
# Some theory and simulation helpers (related to solid state
# physics in general, or ARPES in particular: model
# dispersions, band hybridizations etc.)
#

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

    # Make sure 'pts' and 'klim' have sane values. We will work with
    # full 2D pts / klim layout.

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

    if type(out) is wave.Wave:
        wav = out
    else:
        wav = out.view(wave.Wave)

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
                h1 = 0.5 * (hlist[i] + hlist[j]) + \
                    np.sqrt( (0.5*(hlist[j]-hlist[i]))**2 + V2[i,j]**2 )
                h2 = 0.5 * (hlist[i] + hlist[j]) - \
                    np.sqrt( (0.5*(hlist[j]-hlist[i]))**2 + V2[i,j]**2 )
                hlist[i] = h1
                hlist[j] = h2

    return hlist


#
# Data manipulation helpers (normalization, coordinate transformation etc)
#

def norm_by_noise (data, axis=0, xpos=(None, None), ipos=(None, None),
                   copy=True, smooth=None, stype='gauss', field=False):
    '''
    Normalizes 1D sub-arrays obtained from a N-dimensional ndarray
    along *axis* by the values integrated along
    *axis* in range specified by *pos*.
    
    In ARPES, this is a useful feature for normalizing    
    spectra obtained in a synchrotron facility, which usually
    have significant amount of 2nd-order intensity above
    the Fermi level. Usually, *dim* would be set to represent
    the energy axis, and *pos* should be set to a range
    well above the Fermi level.

    

    Parameters:
      - `data`:  the (multi-dimensional) ndarray containing
        the data. Data poins _will be modified_, so be sure to
        operate on a copy if you wish to retain original data.
        If the data has more than 2 dimensions, the array
        is rotated such that *dim* becomes the last axis, and
        norm_by_noise() iterated on all elements of data,
        recursively until dimension is 2.
    
      - `ipos`: expected to be a tuple containing a
        (from, to) value pair in 
        If *ipos* is 'None', a background auto-detection is
        attempted using gnd_autodetect().
        If *ipos* is a tuple, but either of the *ipos* elements
        is None, the element is replaced by the beginning
        or ending index, respectively.

      - `xpos`: same as *ipos*, only positions are specified
        in axis coordinates of the axis specified by *dim*.
        If present, *xpos* takes precedent over *ipos*.

      - `axis`: is the axis representing the array to be, i.e.
        the one to be normalized.

      - `copy`: if True (the default), then the normalization will
        be performed on a copy of the original array, the original
        data remaining unchanged. Otherwise the original array
        will be modified.

      - `smooth`: Smoothing factor to apply to the data. Can
        be either None (default), 'auto', a number, or an
        (N-1) long tuple of numbers. If it's a tuple, one
        number per dimension is assumed.
        See also: PRO TIP below.

        The meaning of *smooth* depends on the parameter *stype*
        *stype* == 'spline':
          The number is the factor by which the intensity field
          is down-sampled for smoothing (smoothing is done by
          up-sampling the intensity field again to its original value,
          once it has been down-sampled).
          If it is 'auto', for each dimension the proper factor
          is guesst such that after downsampling, approximately
          40 ('ish :-) data points remain.
          None (the default) means no intensity smoothing.
          Down- and up-sampling are performed using 3rd degree
          splines from (scipy.ndimage.map_coordinates).
          Produces polinomial artefacts if intensity distribution
          is very uneven, or data is very noisy.

        *stype* == 'gauss': Smooting is done by a convolution
          of the intensity map with an (N-1)-dimensional Gauss
          profile. In that case, *smooth* contains the Sigma
          parameters of the Gaussian in each dimension.
          Produces artefacts at the border of the image.

      - `stype`: smooth type, either one of 'spline' or 'gaussian'.
        Specifies the type of smoothing.

      - `field`: if True, then the smooth field in original,
        in its down-, and its up-sampled version will also be
        returned (useful for debugging and data quality
        estimates). See also: PRO TIP below.
        
        
        PRO TIP: intensity smoothing can create very strange
        artefacts when dealing with low-intensity / noisy data!
        In that case, either go with the default (i.e. no
        smoothing -- the eye does a better job of "smoothing"
        noisy data anyway), or enable the *field* option and
        _check_ _your_ _normalization_ _fields_ _manually_! ;-)
        In most cases, gaussian smoothing works best with
        experimental data, i.e. produces the most predictible
        amount of artefacts.
    '''

    # rotate axes such that dim is the first
    _data = data.swapaxes(0, axis)

    # we'll be working on Waves all along -- this is
    # because we want to retain axis scaling information
    if copy == True:
        data2 = _data.copy(wave.Wave)
    else:
        data2 = _data.view(wave.Wave)
        data2.setflags(write=True)

    # translate everything to index coordinates,
    # xpos has precedence over ipos
    index = [data2.dim[0].x2i_rnd(x) if x is not None
             else i if i is not None
             else f
             for x, i, f in zip(xpos, ipos, (0, data2.shape[0]-1))]

    # Calculate (N-1)-dim normalization values field. Note that the field
    # itself will be normalized by the number of composing elements along
    # axis. This way, the normalized area will be, by definition, roughly ~1.0
    # Later we can substract 1.0 from the data to have a well defined zero-level :-)
    _norm_field   = data2[index[0]:index[1]].sum(0) / (index[1]-index[0])

    if smooth is not None and stype == 'spline':
        # Smoothing "hack": resample the intensity map twice:
        #  a) first to lower resolution (1/smooth), to get rid of noise
        #  b) then back to original resolution (a * smooth), to get the
        #     proper size again.
        # Use splines for down- and up-sampling, so as little information
        # gets lost as possible

        # a decent estimate: use something like 40(-ish) points per dimension,
        if smooth == 'auto':
            smooth = (np.array(_norm_field.shape) / 40).astype('int') + 3

        # smooth is interpreted as an (N-1)-dim tuple specifying for
        # each dimension the step size (in data points) of smoothing).
        if not hasattr(smooth, "__len__"):
            smooth = np.array([smooth] * _norm_field.ndim)
        else:
            smooth = np.array(smooth)
            if _norm_field.ndim != len(smooth):
                raise ValueError ("Need a smoothing tuple of length %d for an %d-dim array." 
                                  % (_norm_field.ndim))

                
        
        # down-sampling
        _cmp_coord = np.indices((np.array(_norm_field.shape)/smooth) + np.ones(smooth.shape))
        for i, s in zip (_cmp_coord, smooth):
            i *= s
        _cmp_field = spni.map_coordinates (_norm_field, _cmp_coord,
                                           order=3, mode='nearest')[np.newaxis]
            
        # up-sampling (expand again to original size)
        _smooth_coord = np.indices(_norm_field.shape).astype('float32')
        for i, s in zip (_smooth_coord, smooth):
            i /= s
        _smooth_field = spni.map_coordinates (_cmp_field[0], _smooth_coord,
                                              order=3, mode='nearest')

        ## Apply correct scaling to _smooth_field and _tmp_field
        ## (this is just for debugging purposes).
        #b = _cmp_field.view(wave.Wave)
        #for d1, d0, s in zip(b.dim[1:], _norm_field.dim, smooth):
        #    d1.offset = d0.offset
        #    d1.delta = d0.delta * float(s)
        #_cmp_field = b

    elif smooth is not None and stype == 'gauss':

        # auto-select smooth parameter:
        # sigma somewhere close to 1%
        # (i.e. FWHM somewhere at ~2.4%)
        if smooth == 'auto':
            smooth = [] 
            for d in _norm_field.shape:
                sigma = math.floor(float(d)/100.0)
                smooth.append (sigma)
        _smooth_field = spimg.filters.gaussian_filter (_norm_field, smooth)


    else:
        _smooth_field = _norm_field
        
    
    data2 /= np.abs(_smooth_field[np.newaxis])
    data2 -= 1.0

    if field:         # for debugging of code and data... ;-)
        _smooth_wave = _smooth_field.view(wave.Wave)
        for d1,d0 in zip(_smooth_wave.dim[1:], _norm_field.dim):
            d1.lim = d0.lim
        return data2.swapaxes(axis, 0), _norm_field, _smooth_wave

    else:
        return data2.swapaxes(axis, 0)


def get_ref2d_profile (refdata, axis=0, steps=None, width=None, ipos=None, xpos=None):
    '''
    Returns the intensity profile of the spectrum and returns a smoothened
    version of the intensity profile.

    Parameters:
      - refdata: The 2D reference data (ndarray or Wave)
      - axis:    Axis along which to build the profile
                 (the resulting profile will have length=shape[axis].
      - ipos:    Indices along (not axis) between which to integrate.
      - xpos:    Same as *ipos*, only positions are specified in axis
                 coordinates. Works only if input is a Wave(), takes
                 precedence over *ipos* if both are specified.
      - steps:   Number of steps to sustain the profile while spline-
                 smoothing. Must be smaller than shape[axis].
                 The smaller the number, the more drastic the smoothing.
                 None is aequivalent to *step* of shape[axis] means no
                 smoothing at all.
      - width:   Alternative way of specifying the smoothing. Reduces
                 data resolution in (not *axis*) direction by the factor
                 1/*width* by integrating over *width* values, then
                 interpolates the missing values using splines.
    
    Retuns: an 1D 
    '''

    if refdata.ndim != 2:
        raise ValueError ("Wrong dimension %d, expecting 2D data." % refdata.ndim)

    ref = refdata.swapaxes(0,axis)
    
    if xpos is not None:
        if not isinstance (refdata, wave.Wave):
            raise ValueError ("Expecting 'Wave' container with parameter 'xpos'.")
        ipos = (ref.dim[0].x2i(xpos[0]), ref.dim[0].x2i(xpos[1]))

    if ipos is None:
        ipos = (0, ref.shape[1])

    #
    # smoothing will be done by reshaping the array
    #
    if steps is None:
        steps = ref.shape[1]
    width = math.floor(ref.shape[1] / steps)

    _tmp = ref[ipos[0]:ipos[1]].sum(0)[0:steps*width].view(np.ndarray)
    _xin  = np.arange(len(_tmp))
    _xout = np.arange(ref.shape[1])


    # output data -- copying to retain Wave() information
    out = ref[0].copy()

    #print width, steps, width*steps, _xout.shape
    #print _xin.reshape((width,steps)).sum(0).shape
    #print _xin.shape, _xout.shape
    #print _xin[::width]
    #return  _xin[::width], _tmp.reshape((steps,width)).sum(1)/width

    out[:] = spi.UnivariateSpline (_xin[::width],
                                   _tmp.reshape((steps,width)).sum(1)/width)(_xout)[:]

    return out
    


def deg2ky (*args, **kwargs):
    '''
    Converts a 3D wave from the natural coordinates of an ARPES
    measurement (energy / degrees / tilt) into k-space (inverse space)
    coordinates. The 3D wave is assumed to be put together by stacking
    2D waves (ARPES scans in energy x detector coordinates) measured
    at the same photon energy but different tilt angles along
    the 3rd dimension (tilt coordinate).

    It is useful for transforming a kx*ky Fermi Surface maps (single-slices or
    multi-slices, extending into the energy direction).

    For a similar function to be applied to a single 2D arpes slice,
    see *deg2ky_single()*.

    For a similar function to be applied to 2D scans measured
    at different photon energies, resulting in ky*kz scans, see *deg2kz()*.

    Parametes to go with *kwargs*:
    
    - `(unnamed)`:    The 3D data, either as a Wave or as an ndarray,
      containing the deg_tilt*deg_detector*E dependent data
      (i.e. the intensity in dependence of energy, the tilt
      angle, and the detector angle). If data is a Wave,
      axis information is extracted from the wave scaling.
      'e', 'd' and 't' parameters below override internal
      Wave data. If data is an ndarray, then 'e', 'd' and 't'
      parameters are mandatory.
      
    - `eoffs`: Offset to be added to the specified energy axis.
      Typically, if the energy axis is specified with respect to
      EF=0 (Fermi level), then *eoffs* is the value of the Fermi 
      level.
      
    - `axes`:         Combination of the letters 'e', 'd' and 't'
      describing the current axes layout
      in terms of (e)nergy, (d)etector or (t)ilt.
      Default is 'edt'.
      
    - `energy`, `e`: Values to go along with the energy axis.
      If *eoffs* is also specified, it will be added to the
      values of this parameter.
    
    - `detector`, `d`: Values for the detector axis.
    
    - `tilt`, `t`: Values for the tilt axis.
    
    - `degree`:       The spline degree to use for interpolation.
      Default is 3.
      
    - `fill`:         Value to fill invalid data points. Defaults to min(data).


    See notes for **deg2kz()** for supplementary information
    and implementation details.      
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

    if args[0].ndim != 3:
        raise ValueError ("Input has to be a 3D array of values. "
                          "Use numpy.newaxis for proper casting of 2D arrays!")

    # rotate data into 'edt' axis configuration
    axes = kwargs.setdefault('axes', 'edt')
    idata = wave.transpose(args[0], (axes.find('e'), axes.find('d'), axes.find('t')))
    
    E      = _param ('energy',   'e', idata.dim[0].range)
    ideg_d = _param ('detector', 'd', idata.dim[1].range)
    ideg_t = _param ('tilt',     't', idata.dim[2].range)
    eoffs  = _param ('eoffs',    'Ef', 0.0)
    fill   = _param ('fill',     'fill', idata.min())
    degree = _param ('degree',   'deg', 3)
    Phi    = _param ('Phi',      'phi', 4.3523)

    print "Preparing grid... ",
    
    odata  = idata.copy()

    # Energy offset -- usually the value of the Fermi level,
    # if E is specified with respect to the Fermi level Ef.
    # From here, the E axis will show use absolute energies
    # (without the work function Phi)
    E += eoffs
    
    # This is the maximum kinetic energy available *in* *the* *data*
    # (which is usually larger than the Fermi level :-) )
    Ekin_max = E.max()

    # some constants (see deg2kz() for details)
    hsq_2m = 3.80998194907763662131527 # hbar^2 / 2me
    m2_hsq = 0.26246843511741342307750 # 2me / hbar^2
    sqfac  = 0.51231673320067676965494 # sqrt ( 2m / hbar^2)
    
    _d2k = lambda deg, Ekin:  sqfac * np.sqrt(Ekin_max)*np.sin(deg*np.pi/180.0)

    # axes limits of the k-space data
    ik_d_lim = _d2k (np.array([ideg_d[0], ideg_d[-1]]), Ekin_max)
    ik_t_lim = _d2k (np.array([ideg_t[0], ideg_t[-1]]), Ekin_max)
    if isinstance (idata, wave.Wave):
        odata.dim[1].lim = ik_d_lim
        odata.dim[2].lim = ik_t_lim
    

    # rectangular, evenly-spaced grid in k coordinates
    oe, okd, okt = np.broadcast_arrays (E[:,None,None],
                                        np.linspace (start = ik_d_lim[0],
                                                     stop  = ik_d_lim[1],
                                                     num   = len(ideg_d))[None,:,None],
                                        np.linspace (start = ik_t_lim[0], 
                                                     stop  = ik_t_lim[1],
                                                     num   = len(ideg_t))[None,None,:])


    print "done."

    print "Calculating reverse coordinates... ",
    log.info ("Calculating reverse coordinates")
    
    # Reverse transformations: this is where the magic happens
    # (see notes above). Everything else is just house keeping. :-)
    odeg_d = np.arcsin (okd / (sqfac*np.sqrt(oe)) )  # !!! here, this is still in rad
    odeg_t = np.arcsin (okt / (sqfac*np.sqrt(oe)*np.cos(odeg_d )) ) * 180.0/np.pi # deg
    odeg_d *= (180.0 / 3.1415926535)  # conversion rad->deg

    # Some of the coordinates above may end up as NaNs and choke the
    # interpolator. Map the positions and clean-up the data later.
    nan_map = np.isnan(odeg_d) + np.isnan(odeg_t)
    odeg_d[nan_map] = ideg_d[0] # safe polar coordinates to...
    odeg_t[nan_map] = ideg_t[0] # ...use with the interpolator.
    print "done."

    
    print "Interpolating data... ",
    log.info ("Interpolating")
    
    # map_coordinates() takes index coordinates.
    ocoord_index = np.broadcast_arrays (np.arange(idata.dim[0].size)[:,None,None],
                                        idata.dim[1].x2i(odeg_d),
                                        idata.dim[2].x2i(odeg_t))
    spni.map_coordinates (idata, ocoord_index, output=odata.view(np.ndarray),
                          order=degree, mode='constant', cval=fill)

    print "done."
    
    # Clean up points that previously had NaN coordinates.
    odata[nan_map] = fill

    return odata


def deg2ky_single (wav, **kwargs):
    '''
    Quick'n'dirty wrapper around *deg2ky()* that convers a
    single 2D wave by re-packing it into a 3D wave.
    A proper dedicated algorithm could be written for the
    2D case at normal emission (which is easier to do than the
    3D case). But the speed improvement would be marginal, at
    the expense of having to maintain two different
    code paths for the same work.

    Parameters:
       - `wav`: The 2D wave to transform.
       - `axes`: one of 'ed' or 'de', showing the order of the
                 (e)nergy and (d)etector axes.

    For valid *kwargs*, see *deg2ky()*.
    '''
    if kwargs.has_key('axes'):
        kwargs['axes'] += "t"
    return deg2ky (wave.dstack([wav, wav, wav]), **kwargs)[0]


def deg2kz (*args, **kwargs):
    '''
    Converts a 3D wave from the natural coordinates of an ARPES
    measurement (energy / degrees / Ekin-scan) into k-space
    (inverse space) coordinates.The 3D wave is assumed to be put
    together by stacking 2D waves (ARPES scans in
    energy x detector coordinates) measured at different photon
    energies at normal emission (tilt = 0) angles along
    the 3rd dimension (photon energy coordinate).

    It is useful for transforming 2D scans measured at different
    photon energies, resulting in ky*kz Fermi Surface maps. 
        
    For a similar function to be applied to  kx*ky Fermi Surface maps
    see *deg2kz()*. For a similar function to be applied to a single
    2D arpes slice, see *deg2ky_single()*.
    

    Following parametes:
    - `(unnamed)`:   The 3D data, either as a Wave or as an ndarray,
      containing the deg_tilt*deg_detector*E dependent data
      (i.e. the intensity in dependence of energy, the tilt
      angle, and the detector angle). If data is a Wave,
      axis information is extracted from the wave scaling.
      'e', 'd' and 't' parameters below override internal
      Wave data. If data is an ndarray, then 'e', 'd' and 't'
      parameters are mandatory.
      
    - `axes`:        Combination of the letters 'e', 'd' and 'x'
      describing the current axes layout
      in terms of (e)nergy, (d)etector or e(x)citation energy
      (i.e. beam energy E=hv of the excitation beam).
      Default is 'edx'.
      
    - `energy, e`:   Values to go along with the (kinetic or binding)
      energy axis.
    
    - `detector, d`: Values for the detector axis.
    
    - `exbeam, x`:   Values for the excitation energy axis.
    
    - `Phi`:         The work function defined as
      *E_hv = E_final - Phi = E_kin - |E_bind| - Phi*.
      This is usually a material specific constant, in most
      ARPES applications depending on the measurement device.
      Defaults to Phi=4.352(1) eV, which is ok for a
      Scienta R4000 analyzer.
      
    - `v0`: The inner potential of the crystal V0 in eV, usually
      somewhere between 8 and 15 eV. Defaults to 12.5 eV.
      
    - `degree`: The spline degree to use for interpolation. Default is 3.
    
    - `fill`:        Value to fill invalid data points. Defaults to min(data).


    Notes on energy dependence for *polar* -> *kz* transformations
    ==============================================================
    The transformation *polar* -> *kz* depends strongy on the excitation
    energy, as follows:
    
       . *kx* transformation (i.e. detector angle transformation)
         depend on the kinetic energy of the electrons, defined
         as *E_kin = E_beam - |E_bind| - Phi*, where
         *E_beam* is the energy of the excitation beam, *Phi* is
         the work function and *E_bind* is the binding energy
         within the solid, usually related to the Fermi level.
         (For ARPES, this is energy is negative.)
         
       . Additionally, *kz* transformations depend on the excitation
         energy directly via... (here comes V0 etc. need some more
         reading.)

    The input wave, usually, will have two energy axes: the kinetic
    energy axis and the excitation energy axis. The kinetic energy
    axis is usually wrong for all but one single slice, as all
    scans in a *kz* series are performed at different excitation
    energies, and thus their electrons have different kinetic
    energy ranges. The "kinetic energy axis" therefore is interpreted
    more like a "renormalized binding energy". *deg2kz()* makes
    use of it in the following fashion:
    
      a) to determine energy steps from one kx*kz slice to the next.
         Hereby the energy is treated as "kinetic", i.e. stronger
         bound states below E_F are expected to receive _smaller_ 
         energies. Note that this is in contrast some binding
         energy definitions of ARPES, but consistent with energy
         definitions where Inverse PES or a combination of IPES and
         ARPES is employed.

      b) to determine the offset between the excitation energy
         specified by the *z* axis and the kinetic energy axis.
         This offset is calculated once, for the first wave,
         and then retained for the rest of the slices.

    In other words, energy determination for the i-th
    ARPES slice (i.e. data set measured as a *E_kin* x *detector_angle*
    image) is performed as follows:
    
      1) The offset *eoffs* = *E_beam[0]* - *max(E_kin)* is calculated
         once, for the 0-th  ARPES slice.

      2) For all subsequent slices i, *max(E_kin)* = *E_beam[i]* - *offs*
         is used.
    

    Notes on coordinate transformation, also valid for *deg2ky()*
    =============================================================
    
          Conversion from polar coordinates to k-space is performed from the
          input data, which usually lies on a rectangular, regularly spaced
          grid, onto another rectangular, regularly spaced grid in k-space.
          Now, the shape of a rectangular data surface in polar coordinates
          is, usually, non-rectangular in k-space, and vice versa. This
          means that some kind interpolation has to take place at one point
          or another.
          
          Here, two fundamentally different strategies can be employed:
          
             (a) Forward transformation method
             ---------------------------------
             To each input data point in the original, rectangular
             polar coordinate system, the corresponding coordinate
             in k-space is calculated; this yields an "intermediate
             data set", which consists of the original data values
             arranged on an usually non-rectangular grid.
             From the intermediate non-rectangular data set, an
             the final, rectangular k-space grid is obtained by
             interpolation.
             In other words, the input data re-interpreted as being in
             in k-space coordinates, albeit on a non-rectangular
             "random" grid, on which an interpolation to a rectangular,
             regular grid is to be applied using the polar -> k-space
             transformation rules.
             
             An interpolation algorithm from a *random grid* to a
             *regular rectangular grid* is required.

    
             (b) Reverse transformation method
             ---------------------------------
             For every output point in the final, rectangular k-space
             grid the corresponding would-be coordinate in a
             non-rectangular, intermediate polar grid is calculated
             using a *reverse* (i.e. k-space -> polar) transformation.
             Then the input data is interpolated from a rectangular
             regular polar grid to the non-rectangular intermediate
             grid, the latter ultimately corresponding point-by-point
             with the final k-space system.
             In other words, the input data is regarded as regular
             polar-space data, and a non-rectantular polar-space
             (identical to a regular k-space) version of that data
             is obtined by interpolation.
             
             An interpolation algorithm from a *regular rectangular grid*
             to a *random grid* is required.

          Downside of method (a) is the fact that interpolation methods
          from random grids tend to be numerically less efficient and
          more difficult to implement, while method (b) has the downside
          of requiring a closed expression for the reverse transformation.
          Also, sometimes, the reverse coordinate transformation may be
          computationally more expensive than the forward transformation.
          
          However, since in our case the reverse transformation is easy
          to calculate (although computationally slightly more expensive
          than the forward transformation), we're using method (b)
          for **deg2ky()** and **deg2kz()** to gain advantage of the
          numerical stability that comes with it.
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

    if args[0].ndim != 3:
        raise ValueError ("Input has to be a 3D array of values. "
                          "Use numpy.newaxis for proper casting of 2D arrays!")

    # rotate data into 'edt' axis configuration
    axes = kwargs.setdefault('axes', 'edt')
    idata = wave.transpose(args[0], (axes.find('e'), axes.find('d'), axes.find('x')))
    
    E      = _param ('energy',   'e', idata.dim[0].range)
    ideg   = _param ('detector', 'd', idata.dim[1].range)
    iex    = _param ('exbeam',   'x', idata.dim[2].range)
    Phi    = _param ('Phi',      'Phi', 4.352)
    V0     = _param ('V0',       'v0', 12.5)
    fill   = _param ('fill',     'fill', idata.min())
    degree = _param ('degree',   'deg', 3)
    
    odata = idata.copy()

    #
    # Some useful constants:
    #
    # hbar = 1.054571726*10^-34 Js     [kg m^2 s^-1]
    #      = 6.58211928*10^-16 eVs
    # m_e  = 9.10938291*10^-31  kg     [kg]
    #
    # Unit conversion factors
    # eV_J = 1.602176565*10^-19 J/eV (eV <-> Joule conversion factor)
    # m_A  = 10^-10 m/A              (meter <-> Angstroem conversion factor)
    #
    # => hbar^2 / 2m = hbar**2 / 2*m_e        * (1/eV_J)  * (1/m_A^2)
    #                = [kg m^4 s^-2]  [kg^-2] * [eV J^-1] * [m^-2 A^2]
    #                = [kg m^2 s^-2  * m^2]   * [eV J^-1] * [m^-2 A^2]
    #                = [ J * m^2 ]            * [eV J^-1] * [m^-2 A^2]
    #                = [ eV * A^-2 ]
    #                = 3.80998194907763662131527 eV * A^-2
    #
    # Testing: in conjunction with a k calculation:
    #
    #     k = sqrt (2m*E / hbar^2  * sin(...)) 
    #       = sqrt ([eV * A^-2] * eV) = 1/A^-1    :-)
    #
    # Useful constants in the transformation formulae:
    #
    hsq_2m = 3.80998194907763662131527 # hbar^2 / 2me
    m2_hsq = 0.26246843511741342307750 # 2me / hbar^2
    sqfac  = 0.51231673320067676965494 # sqrt ( 2m / hbar^2)

    _rad   = lambda deg: deg * np.pi / 180.0
    _deg   = lambda rad: rad * 180.0 / np.pi

    # forward transformations (needed for boundary calculations)
    _dx2ky = lambda deg, hv:  np.sqrt( m2_hsq * (hv-Phi) ) * np.sin(_rad(deg))
    _dx2kz = lambda deg, hv:  np.sqrt( m2_hsq * ( 
                                                  (hv-Phi) * (1-np.sin(_rad(deg))**2) - V0 
                                                )
                                     )

    # axes limits of the k-space data
    ik_d_lim = tuple( _dx2ky (np.array([ideg[0], ideg[-1]]), max(iex)))
    ik_x_lim = (_dx2kz ( max(ideg[0], ideg[-1]), min(iex[0], iex[-1])),
                _dx2kz ( 0,                      max(iex[0], iex[-1])) )

    # rectangular, evenly-spaced grid in k coordinates;
    kaxis_d = np.linspace (start=ik_d_lim[0], stop=ik_d_lim[1], num=len(ideg))
    kaxis_x = np.linspace (start=ik_x_lim[0], stop=ik_x_lim[1], num=len(iex))

    # for some funny reason, we need to invert det/exb order here...
    okx, okd = np.meshgrid(kaxis_x, kaxis_d)

    # Reverse transformations: this is where the magic happens
    # (see notes above). Everything else is just house keeping. :-)
    oex, odeg   = np.meshgrid (iex, ideg)
    oex  = V0 + hsq_2m*(okx**2 + okd**2)
    odeg = _deg( np.arcsin (np.sign(okd)*np.sqrt (hsq_2m * okd**2 / oex)) )

    # Some of the coordinates above may end up as NaNs. Filter them out
    # before interpolation, and replace the points with 'fill' values later.
    nan_map = np.isnan(odeg) + np.isnan(oex)
    odeg[nan_map] = ideg[0]
    oex[nan_map]  = iex[0]
    
    ocoord_index = [idata.dim[1].x2i(odeg), idata.dim[2].x2i(oex)]

    #
    # Here's another version, commented out. It has the advantage of being
    # easier to read and calculates all points in one single map_coordinates
    # call. However it's slower, it unesessarily interpolates in the 0-th dim.
    #
    ##ocoord_index = np.broadcast_arrays(np.arange(idata.dim[0].size)[:,None,None],
    ##                                   idata.dim[1].x2i(odeg)[None,:,:],
    ##                                   idata.dim[2].x2i(oex)[None,:,:])
    ##odata = spni.map_coordinates (idata, ocoord_index, mode='constant', cval=fill)

    for idat, odat in zip(idata, odata):
        spni.map_coordinates (idat, ocoord_index, output=odat,
                              order=degree, mode='constant', cval=fill)

    # erase values at coordinates that were originally NaNs
    odata[:,nan_map] = fill

    if isinstance (idata, wave.Wave):
        odata = odata.view(wave.Wave)
        odata.dim[1].lim = ik_d_lim
        odata.dim[2].lim = ik_x_lim
    
    return odata


def align2d (a, b, iregion=(0, -1, 0, -1), xregion=None, 
             xshift=None, ishift=None, step=0.5, offset=0,
             stretch=None, weighted=False, fitmode='lsq'):
    '''
    Aligns two 2D waves (*a* and *b*) by shifting them with
    respect to one another systematically over a region of
    maximum *shift* units in *step* steps in either direction,
    and checking the least error squares of the intensity
    quotient field, i.e. optimizes *(a/b' - avg(a/b'))^2*
    for various versions of *b'*, where *b'* is *b* shifted
    in either direction by a certain amount.
    The shifting parameter set which delivers the
    best-matching *b'* (i.e. the one with the smallest error
    squares as specified above) is returned as the optimal
    shifting parameter.

    **Note:** This is a brute-force aligning method, it will compute
              all options and return the best there is within the
              specified grid. It works well, but is very expensive
              (scales with O(n^2)) and will be unusable for too large
              searching regions.

    Parameters:

      * `a`: (Wave-like) the wave to align to.

      * `b`: (Wave-like) the wave to align.

      * `xregion`: (4-tuple: (left, right, top, bottom)) Specifies
        the region in which to check for best
        least-squares matching (axis coordinates,
        takes precedence over *ireg* if specified).
        Any value can be specified as None, in which
        case the corresponding axis limit (offset or end),
        will be substituted.

      * `iregion`: (4-tuple) the region in which to check for best
        least-squares matching (index coordinates)

      * `xshift`:  (2-tuple) Shifting parameters (one per dimension).
        Describes the distance (specified in axis units) 
        in x and y direction to traverse looking for a best
        match. The shifting will take place from 
        *-shift* to *+shift* in either direction.
        If *shift* is None, or either of the components
        are None, it will be replaced by the respective
        width/height of the *region* parameter.
        Has precedence over *ishift*.

      * `ishift`:  (2-tuple) Same as xshift, only in index units.               

      * `step`:   (float or 2-tuple) Step to use per dimension, as a
        factor of the respective dimension granularity (i.e.
        dim-offset). Default is 0.5, which means checking
        with half-index precision.           

      * `offset`: (float or 2-tuple) Known offset of the 
        feature to be matched in *b* over *a*, specified in
        index coordiantes. This parameter can be used to reduce
        the search space in case that a rough estimate of the
        shifting between the two.

      * `stretch`: (number) If specified, the best shifting offset
        will not be extracted from the grid within *xshift* specified
        *step*. Instead, the resulting scoring map will be
        interpolated (i.e. "stretched") by a factor of `stretch`,
        and the ideal shift will be computed from the resulting
        array. With "decent" data, this should gain a fairly
        precise estimate of the shift even with a relatively coarse
        step size. Default (*stretch*=None) disables this feature.

      * `weighted`: (boolean) If True, the error square for each
        data point will be weighted by the data point intensity
        at the respective position. The idea is to have more
        intense data points carry more influence. Use with care,
        it often gives worse results if activated.

      * `fitmode`: (one of 'lsq', 'para-lsq' or 'perp-lsq')
        The mode in which the fitness of a particular shift
        combination is calculated, meaning as follows:
        
          - 'lsq': (the default) means that a least-squares
            sill be calculated over the quotient field of
            all points, minimizing the overall jitter of the
            intensity quotient
            
          - 'para-lsq': the quotient field will be integrated
            parallel to axis 0 (i.e. along axis 1)
            
          - 'perp-lsq': the quotient field will be integrated
            perpendicular to axis 0 (i.e. along axis 0)
            prior to LSQ calculations

        Default is 'lsq', which should be best for most cases.
        'para-lst' and 'perp-lsq' may give better results when
        data is to be aligned only along one direction.

    Region coordinates are all specified in the coordinate system
    of *a*. The waves need not have the same number of points, as
    long as the region of interest *xreg* or *ireg* fits well
    within both waves.
    The procedure does not check or care whether the waves have
    the same granularity (i.e. axis delta parameters). The calling
    instance needs to handle cases with different deltas on its own,
    if required.

    Returns: a 3-tuple *((s0, s1), scores)* where:
        * *(s0, s1)*:  is a 2-tuple with the shifting offsets used
          specified in index units, not including what was specified
          in *offset*. (I.e. if *offset* was (1, 1), and the
          algorithm resulted in a shfting offset of (3, 3), this part
          of the return value will read (3, 3), and *not* (4, 4)! )
        
        * *scores*:    is an numpy.ndarray of the shape [shift[0], shift[0]]
          containing the square-root of the scores normalized
          per data point (lowest score ist best match).
    '''

    # prepare region indices -- working with interger index
    # coordinates, but xreg takes precedence if specified.
    reg = list(iregion)
    if xregion is not None:
        for i in range(len(xregion)):
            idim = i/2
            if xregion[i] is not None:
                reg[i] = a.dim[idim].x2i_rnd(xregion[i])
            else:
                reg[i] = -1 if i%2 else 0

    for i in range(len(reg)):
        if reg[i] == -1 or reg[i] is None:
            reg[i] = a.dim[i/2].size if i%2 else 0

    # calculating step size. working with a 2-tuple internally.
    if not hasattr(step, "__len__"):
        step = (step, step)

    if not hasattr (offset, "__len__"):
        offset = (offset, offset)

    # calculating region size
    if xshift is not None:
        if not hasattr(xshift, "__len__"):
            xshift = (xshift, xshift)
        ishift = [ (x/d.delta if x is not None \
                    else d.size) \
                   for d, x in zip(a.dim, xshift) ]
            
    if not hasattr(ishift, "__len__"):
        ishift = [ishift, ishift]
        
    for i in range(len(ishift)):
        if ishift[i] is None:
            ishift[i] = reg[i*2+1] - reg[i*2]
    ishift = tuple(ishift)

    # slicing indexer for the region of interest
    indexer = tuple ([ slice(reg[2*i], reg[2*i+1]) \
                       for i in range(0, len(reg)/2)])

    # shifting coordinates (1D and 2D)
    _shx = np.arange(-ishift[0], ishift[0], step=step[0])
    _shy = np.arange(-ishift[1], ishift[1], step=step[1])
    if (_shx.size == 0):
        _shx = np.array([0])
    if (_shy.size == 0):
        _shy = np.array([0])
    _shx_2d, _shy_2d = np.broadcast_arrays (_shx[:,np.newaxis], _shy[np.newaxis,:])

    # score matrix
    scores = np.ndarray ([_shx.size, _shy.size], dtype=np.float64)

    # Speed improvement by stripping the info[] block? Probably not much...
    ##b = b.view(np.ndarray).view(wave.Wave)

    # the hard work... :-)
    ar = a[indexer]
    for sx, sy, i in zip(_shx_2d.flat, _shy_2d.flat, range(_shx_2d.size)):
        ##br = wave.regrid(b, {'shift': sx+offset[0]}, {'shift': sy+offset[1]},
        ##                 units='index').view(np.ndarray)[indexer]
        br = wave.regrid(b, {'shift': sx+offset[0]}, {'shift': sy+offset[1]},
                         indexer=indexer, units='index').view(np.ndarray)
        q = ar/br
        err = q-np.average(q)
        if weighted:
            err /= np.abs(ar)

        if fitmode == 'para-lsq':
            err = err.sum(1)
        elif fitmode == 'perp-lsq':
            err = err.sum(0)
            
        sq_err = err**2
        
        score = math.sqrt( sq_err.sum() / ar.size )
        scores.flat[i] = score

    # sort out the best match;
    # now, the first approach is to simply take the smallest
    # value in "scores". however, more elegant would be to
    # interpolate "scores" to a denser number of points
    # and find a minimum of those :-)

    # simple approach:
    best_i = np.argmin(scores)
    best_shift = (_shx_2d.flat[best_i], _shy_2d.flat[best_i])

    # interpolation approach:
    if stretch is not None:
        stretch_factor = stretch
        scores_w = scores.view(wave.Wave)
        scores_w.dim[0].lim = (-ishift[0], ishift[0])
        scores_w.dim[1].lim = (-ishift[1], ishift[1])
        new_shape = (scores.shape[0]*stretch_factor if scores.shape[0] > 1 else 1,
                     scores.shape[1]*stretch_factor if scores.shape[1] > 1 else 1)
        scores_stretch = wave.regrid (scores_w,
                                      {'numpts': new_shape[0]} if new_shape[0] > 1 else None,
                                      {'numpts': new_shape[1]} if new_shape[1] > 1 else None)
        best_i = np.argmin(scores_stretch.view(np.ndarray))
        ix = int(best_i % scores_stretch.shape[0])
        iy = int(best_i / scores_stretch.shape[0])
        best_shift = (scores_stretch.dim[0].i2x(ix),
                      scores_stretch.dim[1].i2x(iy))

    return best_shift, scores


def fermi_guess_efi (data, axis=0, fac=5, cnt=3):
    '''
    Finds the Fermi level in an 1D or 2D set of data
    by integrating along !axis and subsequently averaging
    over increasing distances from one end of the data,
    then from the other, until the current point is 
    larger than the average plus *fac* times the
    standard deviation for at leat *cnt* points.

    The ratio behind the algorithm is that data beyond the Fermi
    level will eventually become numerically zero and stay that
    way (with the exception of some background noise, of course).
    Therefore, the Fermi level is found by walking the data
    starting at one end (the one with less intensity) until
    the current data becomes much larger than the average
    plus a number of standard deviations.

    Works best with low-temperature data.

    Paremters:
    
      - `data`:   (array-like) the 1D or 2D data set.
      
      - `axis`:   (integer) the energy axis in a 2D data set.
        I.e. 2D data will be integrated along momentum
        (!*axis*) and treated as 1D data.
        
      - `fac`: (float) How many standard deviation must a data
        point differ from the average background in order
        to qualify as "relevant intensity" (i.e. Fermi the
        Fermi level).

      - `cnt`: (integer) For how many points to differ. Default is 3.

    Returns the guessed position of the Fermi level.
    '''
    
    if data.ndim == 2:
        ax = int(not axis)
        ndata = data.mean(axis=ax).view(np.ndarray)
    else:
        ndata = data.view(np.ndarray)

    data_max  = ndata.max()
    data_min  = ndata.min()
    data_span = data_max - data_min

    # which end do we need to look at?
    if ndata[0] < ndata[-1]:
        revert = False
    else:
        revert = True
        ndata = ndata[::-1]
        
    found = 0

    # to avoid initial triggering (i.e. at the very first points)
    # by a too small a standard deviation because of too few points,
    # we start by summing up over a certain number of points. *cnt*
    # sounds like a good estimate here...
    for i in range(cnt,ndata.size):
        sub = ndata[0:i]
        val = ndata[i]
        avg = np.mean(sub)
        std = np.std(sub)
 
        # Find the thermally populated end of the data, i.e. where
        # values start significantly increasing.
        if (val) > (avg + std*fac):
            found += 1
            ##print "found", i, found, cnt, val, avg, std
        else:
            found = 0

        # If values have been juming standard deviation for a number
        # of consecutive points, then this must be the Fermi level :-)
        if found >= cnt:
            return i if not revert else ndata.size-i
            

    return None


def fermi_guess_ef (*args, **kwargs):
    '''
    Wrapper for *fermi_guess_efi()* which returns axis units
    instead of index.
    '''
    if len(args) <= 1:
        axis = 0
    elif kwargs.has_key('axis'):
        axis = kwargs['axis']
    else:
        axis = args[1]
    return args[0].dim[axis].i2x(fermi_guess_efi(*args, **kwargs))


def align_stack_ax (dlist, axis=0,
                    xcheck=None, xsearch=None, xcenter=None,
                    icheck=None, isearch=None, icenter=None,
                    step=0.25, debug=False, maxerr=5,
                    fitmode='lsq', stretch=10):
    '''
    Calls **align2d()** on a stack of waves to align them to
    along a specified axis. Useful to align a bunch of 2D ARPES
    waves to their Fermi level by assuming that they change
    only slightly from scan to scan.
    
    The recipe is the following:
    
      1. Take the first pair of waves

      2. Get the approximative positions of the Fermi levels

      3. Shift the 2nd wave to have its Fermi level close
         to the same index as the 1st

      4. Use **align2d()** on the pair

      5. Realign the 2nd wave according to the result of step 4

      6. Proceed to the next wave.

    
    Parameters:
     - `dlist`: (sequence of length N) List of waves to align.

     - `axis`: (integer) the axis along which to search in every
       element (wave) within dlist. Defaults to 0.

     - `xcheck`: (float) Width of the area in which
       the function should verify if data is matching, in
       axis coordinates. Default is 10.

     - `xsearch`: (float) Width of the area in which
       the algorithm should search for matching data, in
       axis coordinates (number of data points).
       Default is 10.

     - `xcenter`: (number or 1D array-like, length same as *dlist*)
       The position at which the identification feature is to be
       expected each wave, in axis units. (In the case of a Fermi
       level align, for example, it could be the approximate position
       of the Fermi level).
       Specifying a single value will use that value for all input
       waves.
       However, specifying a sequence has the advantage that
       known differences in alignment can used to reduce the *xsearch*
       area. A multi-pass use of *align_stack()*, with a coarse step
       over a large area at first, and with subsequently smaller steps
       over smaller areas later, could thus significantly reduce
       computing time ;-)

     - `icenter`: Same as *icenter*, only in axis coordinates.
        If specified, *xcenter* has precedence.

     - `icheck`: Same as *icheck*, only in axis coordinates.
        If specified, *xcheck* has precedence.

     - `isearch`: Same as *isearch*, only in axis coordinates.
        If specified, *xsearch* has precedence.

     - `step`: (floating point) Granularity with which
       to perform the search, as a factor of dimension
       granularity (i.e. as a fraction of the index, or the
       dimension delta, which has the same meaning).

     - `debug`: (boolean) If True, the shift amounts (in index
       coordinates) and the scorings will be returned, too.

     - `maxerr`: (integer) Maximum number of iterations to try and
       work around MemoryError exceptions of scipy :-(
       For some reason, scipy.ndimage.interpolation seems to raise
       MemoryErrors at random points within the computation. Since
       the probability becomes rather high if one tries to align
       a large number of waves in *dlist*, this option will keep
       retrying a step if it fails with MemoryError until it
       succeeds, up to a maximum number of counts. Default is 5.

     - `fitmode`: *fitmode* to pass to align2d(). Default is
       'para-lsq' (different from the align2d() default).

     - `stretch`: Passed to *align2d()*.

    Returns: a tuple *(shifts, scores)*, where:
      - *shifts*: is an array of length *len(dlist)*
        containing the shift of every image to its previous image,
        without accounting for differences specified in *icenter*.
        By definition, *shifts[0]* is 0, since the first image
        does not have any shift.
        If *icenter* was specified as a sequence, then the difference
        between images *i* and *j* should be calculated as:
        
          *(icenter[j]-icenter[i])+(shifts[j]-shifts[i])*

      - *scores*: an array of length *len(dlist)* with the scoring
        values for all shifts. The first element has scoring 0, by
        definition.
    '''

    # calculate known offsets in 'icenter'
    if xcenter is not None:
        if not hasattr(xcenter, "__len__"):
            xcenter = [xcenter] * len(dlist)            
        icenter = np.array([w.dim[axis].x2i(x) for w, x in zip(dlist, xcenter)])
    else:
        if icenter is None:
            icenter = [dlist[0].shape[axis] / 2] * len(dlist)
        elif not hasattr(icenter, "__len__"):
            icenter = [icenter] * len(dlist)

    if xsearch is not None:
        isearch = xsearch / dlist[0].dim[axis].delta
    elif isearch is None:
        isearch = dlist[0].shape[0]/2

    if xcheck is not None:
        icheck = xcheck / dlist[0].dim[axis].delta
    elif icheck is None:
        icheck = isearch

    # List of shifts for every image to the previous one,
    # without accounting for the the offets specified
    # in 'icenter'. By definition, the first image has no
    # shift.
    shifts = [0.0]
    scores = [0.0]
    
    for d0, d1, pos0, pos1 in zip(dlist, dlist[1:], icenter, icenter[1:]):

        if (axis != 0):
            d0 = d0.swapaxes(0,axis)
            d1 = d1.swapaxes(0,axis)
        
        success = False
        trycnt = 0
        while not success and trycnt < maxerr:
            try:
                offs = pos1 - pos0
                shift, score = align2d(d0, d1,
                                       iregion=(pos0-icheck, pos0+icheck, None, None),
                                       ishift=(isearch, 0),
                                       step=(step, None),
                                       offset=(pos1-pos0, 0),
                                       stretch=stretch, fitmode=fitmode)

                score_min = score.min()

                msg =  "%s: shift %f (input offset %f, score=%f, isearch=%f, icheck=%f..%f)"  \
                        % (d1.infs("name"), shift[0], offs, score_min, 
                           isearch, pos0-icheck, pos0+icheck)

                if debug:
                    print msg
                log.info (msg)
                
                #new_data1 = wave.regrid (d1, {'offset': shift[0]+offs}, units='index')
                success = True
                
            except MemoryError:
                msg = "MemoryError, retry %d/%d..." % (trycnt, maxtries)
                log.error ("MemoryError, retry %d/%d..." % (trycnt, maxtries))
                print msg
                success = False
                trycnt += 1

        #print "shapes:", data1.shape, new_data1.shape
        #new_list.append (new_data1)
        
        shifts.append (shift[0])
        scores.append (score_min)

        
    return np.array(shifts), np.array(scores)
    
        
def align_to_ef (dlist, axis=0, step=0.25, search=0.005, check=None,
                 passes=1, guess_opts={}, regrid_opts={}, symrun=False,
                 fitmode='lsq', stretch=10):
    '''
    Aligns a stack of ARPES images to their Fermi levels.
    Uses align_stack(), align2d(), fermi_guess_ef() and
    paul.base.wave.regrid() to do this.
    Basic idea is as follows:
    
      1. Caculate a rough position of the Fermi level
         using fermi_guess_ef().
         
      2. Calculate a fine-grained alignment offset using
         align_stack() (which, itself, uses align2d()).

      3. Remove a systematic offset that will have been
         induced in step 2. (Remember that after step 1,
         usually any shiting intrinsic to the data should
         have been removed, leaving only systematic algorithmic
         drifts.)

      4. Subsequently add all offsets calculated in 1 and 2
         to align the complete stack of data to the same Ef
         level.

    Parameters:

      - `dlist`: (sequence) The waves to align.
      
      - `axis`': (integer) The energy axis in each wave.
      
      - `step`: (number) Fraction of the dimension granularity with
        which to search using align_stack()
        
      - `search`: (number) Search region around Ef in eV, defaults to 5 meV.
      
      - `check`: (number) Checking region around Ef in eV, defaults to *search*x4.
        See align_stack() for a meaning of *search* and *check*

      - `passes`: (number) Number of passes. The 1-pass version represents
        the algorithm as above. On each additional pass, steps (2) and (3)
        are repeated, using the offsets of the previous pass
        as an estimate starting position. On each supplementary pass,
        step size and search distance is being halvened. Sometimes this
        increases accuracy. Use with care.
        
      - `guess_opts`: (dictionary) Optional parameters to pass to fermi_guess_efi().
      
      - `regrid_opts`: (dictionary) Optional parameers to pass to regrid(),
        when doing the final shift of the wave.

      - `symrun`: (boolean) If True, each step (2) will be performed
        twice: once with the regular order of waves, once with the order
        reversed. Depending on the data, this may improve quality and/or
        counteract the average algorithmical drift.

      - `fitmode`: (string) *fitmode* parameter for *align2d()*, default
        is 'lsq'.

    Returns: tuple *(aligned, offsets)* with:
    
      - *aligned*: sequence of length *len(dlist)* containing the
        shifted waves
        
      - *offsets*: the offsets by which the waves were shifted
    '''
    
    # step 1: calculate rough offsets
    guess_opts.update({'axis': axis})
    efi = [fermi_guess_efi(w, **guess_opts) for w in dlist]

    if check is None:
        check = search * 4.0
    
    # 
    fac = 1.0
    while passes > 0:
        log.info ("Calculating offsets for image sequence (%d pass%s to go). This may take a looong while..." 
                  % (passes, "es" if passes > 1 else ""))

        # step 2: calculate fine-grained offsets
        shifts, scores = align_stack_ax (dlist, axis=0, icenter=efi, xsearch=search*fac,
                                         xcheck=check, step=step*fac, fitmode=fitmode,
                                         stretch=stretch)
        if symrun:
            log.info ("Calculating offsets for reverse image sequence...")
            shifts2, scores2 = align_stack_ax (dlist[::-1], axis=0, icenter=efi[::-1],
                                               xsearch=search*fac,
                                               xcheck=check, step=step*fac,
                                               fitmode=fitmode, stretch=stretch)
            shifts[1:] -= (shifts2[1:])[::-1]
            shifts /= 2

        # step 3: remove drift
        drift = shifts[1:].mean()
        shifts[1:] -= drift

        log.info ("Average drift: %f per scan." % drift)

        # step 4: caculate cumulative offset to first wave
        offs_sum = 0.0
        offsets = [0.0]
        for e0, e1, s0, s1 in zip(efi[0:-1], efi[1:], shifts[0:-1], shifts[1:]):
            offs_prev = (e1-e0) + (s1-s0)  # local offset, to previous wave
            offs_sum += offs_prev          # global offset, to the first wave
            offsets.append (offs_sum)

        efi += offsets
        passes -= 1
        fac *= 0.5



    # (step 5: calculate output waves :-) )
    log.info ("Regridding waves...")
    out = []
    for iw, offs in zip(dlist, offsets):
        log.info ("%s: offset %f" % (iw.infs("name"), offs))
        regrid_opts.update ({'units': 'index'})
        ow = wave.regrid (iw.swapaxes(0, axis), {'shift': offs}, **regrid_opts)
        out.append (ow)

    return out, np.array(offsets)


def align_kz_uru2si2 (data, ss_index='auto', ef_index='auto',
                      deg_search=5.0, e_search=0.005,
                      reverse_energy='auto', e_axis=0):
    '''
    Convenience function to align a stack of URu2Si2 scans in
    taken with a kz succession in both angle and energy direction.
    This is a highly specialized function, possibly not useful
    for another situation.

    The recipe is the following:

      1. Perform an intensity normalization based on nose above Ef
      2. Align make a rough Ef alignment
      3. Make an angle alignment (symmetric run)
      4. Make a fine-grained Ef alignment (symmetric run)

    Returns as a tuple: (aligned stack, energy shifts, momentum shifts).
    Shifts are cumulative (i.e. absolute), relative to the first wave.
    '''

    deg_axis = int (not e_axis)

    if ef_index == 'auto':
        ef_index = fermi_guess_efi (data[0], axis=e_axis)

    if reverse_energy == 'auto':
        foo = data[0].swapaxes(e_axis, 0)
        reverse_energy == data[0][0] > data[0][-1]

    norm_pos = (0, ef_index/2) if not reverse_energy \
      else (ef_index + (data[0].shape[e_axis] - ef_index)/2, data[0].shape[e_axis]-1)


    # intensity profile normalization along the momentum axis
    # (i.e. integrate along energy axis).
    # this function will substact -1 from the overall intensity
    # to get rid of the noise above Ef. however, various steps
    # in align2d() misbehave if intensity is too close to 0,
    # so we'll add back the 1 and remove it later.
    log.info ("Step 1: Intensity normalization at indices %d...%d" % norm_pos)
    norm = [norm_by_noise(w, axis=e_axis, ipos=norm_pos)+1 for w in data]
    
    log.info ("Step 2: First Ef alignemnt (this will take a while)")
    alg1, shf1 = align_to_ef (norm, axis=e_axis, step=0.5, search=e_search)
    
     
    if ss_index == 'auto':
        ss_index = np.argmax(alg1[0]) % alg1[0].shape[deg_axis]
    log.info ("Step 3: Angle alignment using SS at %d "\
              "(this will also take a while)" % ss_index)
    ss_width = round(float(deg_search) / alg1[0].dim[deg_axis].delta)
     
    log.info ("Step 3a: (normal sequence aligment)")
    shf1_deg, sc1_deg = align_stack_ax (alg1, axis=deg_axis, icenter=ss_index,
                                        isearch=ss_width, icheck=ss_width, step=1)
    #print shf1_deg
     
    log.info ("Step 3b: (reverse sequence alignment)")
    shf2_deg, sc2_deg = align_stack_ax (alg1[::-1], axis=deg_axis, icenter=ss_index,
                                        isearch=ss_width, icheck=ss_width, step=1)
    #print shf2_deg
     
    _rel_shf = (shf1_deg[1:] - shf2_deg[1:][::-1]) / 2.0
    shf_deg = np.cumsum(np.insert(_rel_shf, obj=0, values=np.array([0]), axis=0))
    
    #print shf_deg, shf_deg.shape

    #shf_deg = np.zeros (len(data))  # DEBUG
    #alg1 = data                     # DEBUG
    
    if deg_axis == 1:
        log.info ("Step 3c: (re-gridding 2nd dimension)")
        deg = [wave.regrid(w, None, {'shift': s}, units='index')
               for w, s in zip(alg1, shf_deg)]
    else:
        log.info ("Step 3c: (re-gridding 1st dimension)")
        deg = [wave.regrid(w, {'shift': s}, None, units='index') 
               for w, s in zip(alg1, shf_deg)]

        
    log.info ("Step 4: Fine-grained Ef alignment (this will take even longer)")
    alg2, shf2 = align_to_ef (deg, axis=e_axis, step=0.25, symrun=True)


    # return the aligned samples -- but don't forget to remove
    # the "1" that was added in the normalization step :-)
    return [w-1 for w in alg2], shf1+shf2, shf_deg, norm, alg1, shf1_deg, shf2_deg
    
    

if __name__ == "__main__":

    log = logging.getLogger ("paul")

    fmt = logging.Formatter('%(levelname)s: %(funcName)s: %(message)s')
    ch  = logging.StreamHandler()
    ch.setFormatter(fmt)
    ch.setLevel (logging.DEBUG)
    log.addHandler (ch)
    log.setLevel (logging.DEBUG)

    #foo = e(pts=(5, 4), mrel=-2)
    #foo = e()
    #pprint (foo)
    #print foo.info['axes']
    #e(pts=5, klim=1.0)
    #e(pts=5, llim=(-1.0, 0.5))
    #e(pts=(4, 5), klim=((1.0), (1,0)))
