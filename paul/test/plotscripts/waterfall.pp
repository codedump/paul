#!/usr/bin/python

'''
Simple script to display a 2D wave as a waterfall diagram.
'''

import paul.toolbox.atrix as atrix
import paul.base.wave as wave
import paul.toolbox.mpltrix as mpltrix
import numpy as np
import matplotlib as mpl
import pprint as pp

import logging
log = logging.getLogger (__name__)

#reload(mpltrix)
#reload(wave)
#reload(atrix)

def populate(*args, **kwargs):
    can = kwargs['can']
    can.reset()

    comp_ax = 0
    comp_step = 1
    comp_norm = False

    offs = (0, 0.03)
    xlim = (0, 0)
    ylim = (0, 0)
    
    xlim = (11.625, 11.66)
    ylim = (-10, 12)
    ylim = (2, 38)

    if comp_ax:
        _wav = kwargs['wav'][0]
    else:
        _wav = kwargs['wav'][0].swapaxes(0,1)

    can.axes.axvline (11.6485, ls=':')

    #print "First input element:"
    #pp.pprint (_wav)
    
    # create waterwall compression
    wav = atrix.ncomp (_wav/_wav.infv('FDD', 'V_max', default=_wav.max()),
                       axis=comp_ax,
                       step=comp_step,
                       intg=comp_step,
                       norm=comp_norm) 
    
    kwargs['wav'] = [wav]

    # surpress data above Ef+5*kT if sample was FDD normalized
    if ('FDD' in wav.info):
        wav = wav(slice(wav.dim[0].offset, wav.infv('FDD', 'kT', default=0)*5), slice(None))



    #print "Compressed element:"
    #pp.pprint (wav)
    #print "Is array:", isinstance(wav, np.ndarray)
    #print "Is wave:", isinstance(wav, wave.Wave)
    #print "Limits:", wav.lim, _wav.lim

    # plot the data
    kwargs['lines'] = mpltrix.imwater (can.axes, wav,
                                       axis=comp_ax, 
                                       offs=offs, 
                                       xlim=xlim, 
                                       ylim=ylim)

    #print "axis limits:", w.dim[comp_ax].lim
    #can.axes.set_ylim (11.15, 11.35)
    #can.axes.set_ylim (w.dim[0].lim)

    decorate (*args, **kwargs)
    can.draw()
    

def decorate (*args, **kwargs):
    pass
