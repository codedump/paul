#!/usr/bin/python

'''
Simple script to display a 2D wave as a waterfall diagram.
'''

from paul.toolbox.mpltrix import imwater
from paul.toolbox.atrix import ncomp

import matplotlib as mpl

import logging
log = logging.getLogger (__name__)

def populate(*args, **kwargs):
    can = kwargs['can']
    can.reset()

    comp_ax = 0
    comp_step = 5
    comp_norm = False

    offs = (0, 0.1)
    xlim = (0, 0)
    ylim = (0, 0)

    # create waterwall compression
    wav = [ ncomp (w/w.infv('FDD', 'V_max', default=w.max()),
                   axis=comp_ax,
                   step=comp_step,
                   intg=comp_step,
                   norm=comp_norm) 
            for w in kwargs['wav']]
    
    kwargs['wav'] = wav

    # surpress data above Ef+5*kT if sample was FDD normalized
    if ('FDD' in wav[0].info):
        w = wav[0](slice(wav[0].dim[0].offset,
                         wav[0].infv('FDD', 'kT', default=0)*5),
                   slice(None))
    else:
        w = wav[0]

    # plot the data
    kwargs['lines'] = imwater (can.axes, w,
                               axis=comp_ax, 
                               offs=offs, 
                               xlim=xlim, 
                               ylim=ylim)

    decorate (*args, **kwargs)
    can.draw()
    

def decorate (*args, **kwargs):
    pass
