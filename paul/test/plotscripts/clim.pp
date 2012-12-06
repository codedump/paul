#!/usr/bin/python

def decorate (*args, **kwargs):
    ax = kwargs['can'].axes
    w = kwargs['wav']
    
    clim = (w.min()+(w.max()-w.min())*0.00, w.min()+(w.max()-w.min())*0.90)

    ax.images[0].set_clim (clim)
    ax.set_ylim (w.dim[0].min, w.dim[0].max)
    pass
    
