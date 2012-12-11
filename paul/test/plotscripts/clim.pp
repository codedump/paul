#!/usr/bin/python

def populate (*args, **kwargs):
    can = kwargs['can']
    can.reset()
    wav = kwargs['wav'][0]
    kwargs['axes'] = can.axes
    can.axes.imshow (wav, extent=wav.imlim, interpolation='none')
    decorate (*args, **kwargs)
    can.draw()
    

def decorate (*args, **kwargs):
    ax = kwargs['can'].axes
    w = kwargs['wav'][0]
    
    clim = (w.min()+(w.max()-w.min())*0.00, w.min()+(w.max()-w.min())*0.90)

    ax.images[0].set_clim (clim)
    ax.set_ylim (w.dim[0].min, w.dim[0].max)
    pass
    
