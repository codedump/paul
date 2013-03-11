#!/usr/bin/python

import paul.toolbox.arpes as arpes
import pprint as pp
import numpy as np

reload (arpes)

def populate (*args, **kwargs):
    kwargs['can'].reset()
    fig = kwargs['fig']
    wav = kwargs['wav']
    ax = fig.axes[0]

    if not wav[0].info.has_key ('FDD'):
        #clim = (np.nanmin(wav[0]), np.nanmax(wav[0]))
        foo = arpes.norm_by_fdd (wav[0], axis=0, Ef=0, kT=0.001500)
        ax.plot (wav[0].dim[0].range, foo)
    else:
        foo = wav[0]
        
    clim = (foo.infv('FDD', 'V_min'),
            foo.infv('FDD', 'V_max'))
    
    ax.imshow (foo, extent=wav[0].imlim)
    ax.images[0].set_clim (clim)
    ax.set_ylim (ax.get_ylim()[::-1])
    
    decorate (wav=[foo], axes=fig.axes)
    
    kwargs['can'].draw()


def decorate (*args, **kwargs):
    '''
    Valid kwargs: can, wav, fig, axes. 
    '''
    pass
