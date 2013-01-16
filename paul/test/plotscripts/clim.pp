#!/usr/bin/python

import numpy as np
from paul.base.wave import *

import matplotlib as mpl

def populate (*args, **kwargs):
    can = kwargs['can']
    can.reset()
    w = kwargs['wav']

    if isinstance(w, list) or isinstance(w, tuple):
        wav = w[0]
        log.debug ("List")
    else:
        wav = w if isinstance(w, Wave) else w.view(Wave)
        log.debug ("Array")

    kwargs['axes'] = can.axes
    can.axes.imshow (wav, extent=wav.imlim, interpolation='none')
    decorate (*args, axes=can.axes, wav=wav)
    can.draw()
    

def decorate (*args, **kwargs):
    ax = kwargs['axes'].axes
    w = kwargs['wav']
    
    clim = (w.min()+(w.max()-w.min())*0.00, w.min()+(w.max()-w.min())*0.55)

    ax.images[0].set_clim (clim)
    ax.set_ylim (w.dim[0].min, w.dim[0].max)
    ax.axvline (0, ls=':')
    ax.axhline (0, ls=':')

    ax.xaxis.set_minor_locator (mpl.ticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator (mpl.ticker.MultipleLocator(5))

    ax.xaxis.grid (True, which='both')
    
    
