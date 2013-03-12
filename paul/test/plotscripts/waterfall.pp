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
import math

import logging
log = logging.getLogger (__name__)

#reload(mpltrix)
#reload(wave)
#reload(atrix)

def populate(*args, **kwargs):
    can = kwargs['can']
    can.reset()

    comp_ax = 0
    comp_step = 7
    comp_norm = False

    offs = (0, -0.02)
    xlim = (0, 0)
    ylim = (0, 0)
    
    #xlim = (11.625, 11.66)
    #xlim = (27.60, 27.66)
    #xlim = (-0.04, 0.01)
    
    #ylim = (-10, 12)
    #ylim = (2, 38)
    #ylim = (-0.5, 0.7)
    #ylim = (-0.1, 0.1)

    # swap axes to move the comp-axis to 0
    input_wave = kwargs['wav'][0]
    val_max = input_wave.infv('FDD', 'V_max', default=input_wave.argmax())
    if comp_ax:
        _wav = input_wave / val_max
    else:
        _wav = input_wave.swapaxes(0,1) / val_max

        # Surpress data above Ef+5*kT if sample was FDD normalized
    if ('FDD' in _wav.info):
        indexer = [slice(None),slice(None)]
        kt_lim   = input_wave.infv('FDD', 'kT', default=0)*5
        kt_index = math.ceil(input_wave.dim[0].x2i(kt_lim))
        indexer[int(not comp_ax)] = slice(0, kt_index)

        _wav = _wav[indexer]
    
    # waterfall compression
    wav = atrix.ncomp (_wav, step=comp_step, intg=comp_step, norm=comp_norm) 
    
    lines = mpltrix.imwater (can.axes, wav, offs=offs, xlim=xlim, ylim=ylim)

    decorate (*args, wav=[wav], axes=can.axes, lines=lines)
    can.draw()
    

def decorate (*args, **kwargs):
    ax = kwargs['axes']
    ax.axvline (19.6482, ls=':')
    #ylim = ax.get_ylim()
    #ax.set_ylim (ylim[1], ylim[0])
