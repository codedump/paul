#!/usr/bin/python

'''
Simple script to display a 2D wave as a waterfall diagram.
'''

import paul.toolbox.atrix as atrix
import paul.base.wave as wave
import paul.toolbox.mpltrix as mpltrix
import paul.toolbox.arpes as arpes
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

    input_wave = kwargs['wav'][0]

    comp_ax = 0
    comp_step = 5
    comp_norm = False
    do_fdd = False
    do_fdd = True

    offs = (0, 0.15)
    xlim = (0, 0)
    ylim = (0, 0)
    
    #xlim = (11.625, 11.66)
    #xlim = (27.60, 27.66)
    xlim = (-0.06, 0.02)
    
    ylim = (-0.3, 0.6)
    #ylim = (2, 38)
    #ylim = (-0.5, 0.7)
    #ylim = (-0.1, 0.1)

    # FDD-normalization, if required
    if do_fdd and not kwargs['wav'][0].info.has_key ('FDD'):
        input_wave = arpes.norm_by_fdd (input_wave, axis=0, Ef=0, kT=input_wave.infv('kT'))

    # calculate data limits and normalization factor
    data_lim = (input_wave.infv('FDD', 'V_min', default=np.nanmin(input_wave)),
                input_wave.infv('FDD', 'V_max', default=np.nanmax(input_wave)))
    norm = max(data_lim)    # normalize to maximum intensity (typically SS)
    norm = input_wave(0,0)  # normalize to value at E_F.
    
    # swap axes to move the comp-axis to 0
    work_wave = (1.0/norm) * (input_wave if comp_ax == 1 else input_wave.swapaxes(0,1))
    

    # Surpress data above Ef+5*kT if sample was FDD normalized
    if ('FDD' in work_wave.info):
        indexer = [slice(None),slice(None)]
        kt_lim   = input_wave.infv('FDD', 'kT', default=0)*5
        kt_index = math.ceil(input_wave.dim[0].x2i(kt_lim))
        indexer[int(not comp_ax)] = slice(0, kt_index)
        work_wave = work_wave[indexer]
    
    # waterfall compression
    displ_wave = atrix.ncomp (work_wave, step=comp_step, intg=comp_step, norm=comp_norm) 

    # displaying
    lines = mpltrix.imwater (can.axes, displ_wave, offs=offs, xlim=xlim, ylim=ylim)

    # pass on to graph decoration...
    decorate (*args, wav=[displ_wave], axes=can.axes, lines=lines)
    can.draw()
    

def decorate (*args, **kwargs):
    ax = kwargs['axes']
    ax.axvline (26.6460, ls=':')
    #ylim = ax.get_ylim()
    #ax.set_ylim (ylim[1], ylim[0])
