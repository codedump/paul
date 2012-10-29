#!/usr/bin/python

import math
import numpy as np

import logging
log = logging.getLogger(__name__)

def populate (*args, **kwargs):
    '''
    This one will just populate the graph -- put all
    beautifying elements in decorate()!
    '''
    fig = kwargs['can']
    wav = kwargs['wav']

    log.debug ("%d waves to plot" % len(wav))

    #new_wav = []
    #for i in range(0,len(wav),2):
    #    log.debug ("%d %d %d" % ( i, i+1, len(wav)))
    #    w0 = wav[i]   / wav[i].max()
    #    w1 = wav[i+1] / wav[i+1].max()
    #    wd = w1 / w0
    #    wd.info['name'] = 'LV/LH ratio'
    #    new_wav.append(w0)
    #    new_wav.append(w1)
    #    new_wav.append(wd)
    #old_wav = wav
    #wav = new_wav

    fig.clear()
    prows = math.ceil(math.sqrt( len(wav) ))
    pcols = math.ceil(len(wav) / prows)

    log.debug ("Figure: %dx%d, %d plots" % (prows, pcols, len(wav)))
    
    fig.axes = []
    for i in range(len(wav)):
        w = wav[i]
        log.debug ("Subplot %d (%dx%d): %s" % (i+1, prows, pcols, 
                                               w.info['name']))
        ax = fig.fig.add_subplot (prows, pcols, i+1)
        ax.imshow(w, aspect='auto', extent=w.imlim)
        decorate (wav=[w], axis=ax)
        if not ((i+1) > (len(wav)-pcols) != 0):
            ax.set_xlabel (r'')
        if not ((i) % pcols == 0):
            ax.set_ylabel (r'')
        fig.axes.append(ax)
    
    fig.fig.subplots_adjust(          top=0.95,
                            left=0.1,             right=0.95,
                                      bottom=0.1,
                            hspace=0.2, wspace=0.2)
    

def decorate(*args, **kwargs):
    '''
    Set decorative elements, called once for each (sub-)plot.
    '''
    ax = kwargs['axis']
    wav = kwargs['wav']
    w = wav[0]

    if w.info.has_key ('T'):
        temp = w.info['T'][0]+" K"
    else:
        temp =  'n/a'
    
    ax.set_title ("Sample %s, T = %s" % (w.info['name'], temp))
    ax.axhline(ls=':', color=(1, 1, 1, 0.8))
    #ax.set_ylim (w.dim[0].lim) # make graph correct side up
    ax.set_ylim (w.dim[0].lim[0], 0.02)

    if (w.dim[1].lim[1] > 1):
        ax.set_xlim (-10, 10)
    else:
        ax.set_xlim (-0.3, 0.3)

    ax.set_xlabel (r'$k_{\parallel} (\AA^{-1})$')
    ax.set_ylabel (r'$E$ (eV)')



#
# default wave list, in case we have nothing else to plot
#
