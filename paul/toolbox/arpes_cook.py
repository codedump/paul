#!/usr/bin/python

'''
Creates the normalized and FDD-ed versions of the data
from the raw versions in specified subdirectories.
'''

import paul.base.wave as wave
import paul.loader.igor as igor
import paul.toolbox.arpes as arpes

import numpy as np
import sys, os, argparse

def cook_scienta_get_ef (data, axis=0):
    '''
    Returns the position of the Fermi level in 'data'
    by analysing Scienta metadata.
    '''
    # this is the total energy offset of the scale:
    # i.e. usually, the scale will have an energy relative to
    # E_F = 0, but the data was measured at an E_kin,
    # i.e. relative to the vacuum energy. We need to know
    # the absolute value of the kinetic energy for the
    # polar -> kspace conversion.
    return ((data.infv('SES', 'Low Energy')  - idata.dim[axis].lim[0]) +
            (data.infv('SES', 'High Energy') - idata.dim[axis].lim[1])) / 2
    
def cook_ax_shift (data, offsets=None):
    '''
    Shifts axis offset(s) by the amount specified by the tuple 'offsets'.
    '''

    if offsets == None:
        offsets = [0] * data.ndim

    odata = data.copy()

    print "* Axis offsets:",
    for dim, offs in zip (data.dim, offsets):
        print offs,
        dim.offset += offs
    print ""
    return odata
    

def cook_norm_noise (data, axis=0, norm_range=None):
    '''
    Normalizes data along along 'not axis' by the intensity
    integrated along 'axis' in the range 'norm_range'.
    '''
    if norm_range is None:
        norm_range = (0, data.dim[axis].size)
        
    print "* Normalizing intensity along %d, range %d...%d" % (axis, norm_range[0], norm_range[1])
    return arpes.norm_by_noise (data, axis=axis, ipos=norm_range)
    

def cook_deg2ky (data, axis=0, tilt=0.0, eoffs=0.0):
    '''
    Conversion polar -> k coordinates, where 'axis' is the energy axis.
    'tilt' specifies the tilt angle at which the data was measured.
    '''
    print "* Converting to k-coordinates (total offset: %f eV)" % eoffs_total
    axtag = ('ed' if axis==0 else 'de')
    return arpes.deg2ky_single (data, axes=axtag, tilt=tilt, eoffs=0.0)

    
def cook_norm_gnd(data, axis=0, norm_range=None):
    '''
    Calculates and substracts background along 'axis' in the range 'norm_range'.
    '''
    if norm_range is None:
        norm_range = (0, data.dim[axis].size)
        
    gnd = np.average(data[norm_range[0]:norm_range[1]]) if axis \
      else np.average(data.swapaxes(0,1)[norm_range[0]:norm_range[1]])

    print "* Substracting background (%d...%d: %f)" % (norm_range[0], norm_range[1], gnd)
    return data-gnd

    
def cook_fdd(data, axis=0, kT=0, Ef=0):
    '''
    Normalizes data by the Fermi-Dirac distribution with specified
    parameters kT and E_F. 'axis' designates the energy axis.
    '''
    print "* Normalizing to Fermi-Dirac distribution (kT: %f eV)" % kT
    return arpes.norm_by_fdd (data, axis=axis, kT=kT)
    

def cook_from_beamline (data, eax=0, eoffs=0, doffs=0, tilt=None, kT=None, norm_range=None):
    '''
    Performs a complete data preparation cycle on beam-line data:
      - axis shifting
      - normalization (intensity)
      - background substraction
      - polar -> ky conversion
      - FDD normalization

    and returns a tuple with:
      - the normalized wave
      - the ky-converted wave
      - the FDD normalized wave
    '''

    
    # degree offset update (energy offset needs to be performed after the deg2ky transform)
    if doffs is not None:
        data   = cook_ax_shift (data, offsets=((0, -doffs) if eax==0 else (-doffs, 0)))
    if eoffs is not None:
        data   = cook_ax_shift (data, offsets=((-eoffs, 0) if eax==0 else (0, -eoffs)))

    if norm_range is not None:
        data   = cook_norm_noise (data, axis=eax, norm_range=norm_range)
        w_norm = data.copy()

    # polar -> ky transformation:
    if tilt is not None:
        data   = cook_deg2ky (data, axis=eax, tilt=tilt,
                              eoffs=(0.0 if eoffs is None else eoffs))
        w_ky   = data.copy()

    # fdd normalization
    if norm_range is not None:
        data   = cook_norm_gnd (data, axis=eax, norm_range=norm_range)
    if kT is not None:
        data   = cook_fdd (data, axis=0, kT=kT, Ef=0)
        w_fdd  = data.copy()

    
    return data, w_norm, w_ky, w_fdd
    
    

def _make_from_raw (idata, eoffs=0, doffs=0, kT=None, norm_range=None):
    '''
    DEPRECATED.
    '''
    
    axis = 0
    idata.dim[axis].offset -= eoffs
    phi=4.3523 # Work function

    print "* Substracting energy offset %f" % eoffs

    # this is the total energy offset of the scale:
    # i.e. usually, the scale will have an energy relative to
    # E_F = 0, but the data was measured at an E_kin,
    # i.e. relative to the vacuum energy. We need to know
    # the absolute value of the kinetic energy for the
    # polar -> kspace conversion.
    eoffs_total = ((idata.infv('SES', 'Low Energy')  - idata.dim[0].lim[0]) +
                   (idata.infv('SES', 'High Energy') - idata.dim[0].lim[1])) / 2

    print "* Normalizing intensity"
    deg     = arpes.norm_by_noise (idata, axis=axis, ipos=norm_range)
    
    deg.dim[int(not axis)].offset -= doffs

    print "* Converting to k-coordinates (total offset: %f eV)" % eoffs_total
    axtag = ('ed' if axis==0 else 'de')
    ky      = arpes.deg2ky_single (deg, axes=axtag, eoffs=eoffs_total)

    gnd = np.average(ky[norm_range[0]:norm_range[1]])
    print "* Substracting background (%d...%d: %f)" % (norm_range[0], norm_range[1], gnd)
    ky -= gnd

    print "* Normalizing to Fermi-Dirac distribution (kT: %f eV)" % kT
    oky_fdd = arpes.norm_by_fdd (ky, axis=axis, kT=kT)
    
    return deg, ky, oky_fdd
    

if __name__ == "__main__":
    
    if len (sys.argv) < 5:
        print "  Usage: %s IBW-File E_offset deg_offset kT" % sys.argv[0]
        sys.exit()

    infile = sys.argv[1]
    inw = igor.load (infile)
    
    print
    print "Input file:", infile
    
    outw, odeg, oky, oky_fdd = cook_from_beamline (inw, eax=0,
                                                   eoffs=float(sys.argv[2]),
                                                   doffs=float(sys.argv[3]),
                                                   kT=float(sys.argv[4]),
                                                   norm_range=(130, 150))
    owaves = [odeg, oky, oky_fdd]
    ofiles = [infile.replace(".ibw", "g.ibw"),
              infile.replace(".ibw", "g_ky.ibw"),
              infile.replace(".ibw", "g_ky_fdd.ibw")]

    for of, ow in zip (ofiles, owaves):
        print "* Writing %s" % of
        igor.wave_write (ow, of)
