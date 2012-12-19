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

    fig.clear()
    pcols = math.ceil(math.sqrt( len(wav) ))
    prows = math.ceil(len(wav) / pcols)
    
    fig.axes = []
    for i in range(len(wav)):
        w = wav[i]
        log.debug ("Subplot %d (%dx%d): %s" % (i+1, prows, pcols, 
                                               w.info['name']))
        ax = fig.fig.add_subplot (prows, pcols, i+1)
        ax.imshow(w, aspect='auto', extent=w.imlim, interpolation='none')
        decorate (wav=[w], axis=ax)

        # keep xlabels only for the bottom row
        if (len(wav)-pcols) > i:
            ax.set_xlabel (r'')

        # keep ylabels only for the left column
        if not ((i % pcols) == 0):
            ax.set_ylabel (r'')

        fig.axes.append(ax)
    
    fig.fig.subplots_adjust(          top=0.94,
                            left=0.1,             right=0.95,
                                      bottom=0.12,
                            hspace=0.2, wspace=0.2)
    

def decorate(*args, **kwargs):
    '''
    Set decorative elements, called once for each (sub-)plot.
    '''
    ax = kwargs['axis']
    w  = kwargs['wav'][0]
    
    ax.set_title (r'%s   T=%.0f K' % (w.infs('name'), w.infv('T')))
    ax.axhline(ls=':', color=(1, 1, 1, 0.8))
    ax.set_ylim (w.dim[0].lim) # make graph correct side up

    ax.set_xlabel (r'$k_{\parallel} (\AA^{-1})$')
    ax.set_ylabel (r'$E$ (eV)')

    # adjust color limits, in case we have an FDD sample
    if not np.isnan(w.infv('FDD', 'V_max')):
        ax.images[0].set_clim (w.infv('FDD', 'V_min'),
                               w.infv('FDD', 'V_max'))
