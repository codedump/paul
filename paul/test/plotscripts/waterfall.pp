#!/usr/bin/python

from matplotlib.collections import LineCollection
import matplotlib as mp
import numpy as np
import pprint as pp
from paul.toolbox.slicers import WaterfallSlicer

import logging
log = logging.getLogger (__name__)

class GlobalVars:
    pass

def reload(*args, **kwargs):
    '''
    called in case of reload -- we don't want to lose our GUI settings
    '''
    vars = kwargs['vars']
    can  = kwargs['can']
    win  = kwargs['win']
    global G
    G = GlobalVars()
    G.gui = vars.gui
    G.gui.sliced.connect (waterfall)
    G.fig = can

def unload(*args, **kwargs):
    global G
    vars = kwargs['vars']
    can  = kwargs['can']
    win  = kwargs['win']
    vars.gui = G.gui
    vars.fig = G.can
    G.gui.sliced.disconnect (waterfall)

def init(*args, **kwargs):
    log.debug ("Initializing waterfall plotscript")
    vars = kwargs['vars']
    can  = kwargs['can']
    win  = kwargs['win']
    
    global G
    G = GlobalVars()
    G.gui = WaterfallSlicer(tparent=win)
    G.gui.sliced.connect (waterfall)
    G.fig = can


def populate(*args, **kwargs):
    global G
    G.fig = kwargs['can']
    wav   = kwargs['wav']
    G.gui.slice(wav[0] / wav[0].max())


def waterfall (ws):
    global G

    G.fig.reset()

    if ws.ndim is not 2:
        log.error ("Expecting 2D data set, got %dD" % (ws.ndim))
        return

    ax2    = G.gui.val_axis.value() # axis to retain
    ax1    = not ax2
    xoffs = G.gui.xoffset
    yoffs = G.gui.yoffset
    x = np.arange(start=ws.dim[ax1].offset,
                  stop=ws.dim[ax1].end,
                  step=ws.dim[ax1].delta)

    log.debug ("Axis limits: %f, %f" % (ws.dim[ax1].min, ws.dim[ax1].max))

    xmin = ws.dim[ax1].min
    xmax = ws.dim[ax1].max
    ymin = ws.min() + (yoffs*ws.shape[ax2]) * (yoffs<0)
    ymax = ws.max() + (yoffs*ws.shape[ax2]) * (yoffs>0)
    G.fig.axes.set_xlim (xmin, xmax)
    G.fig.axes.set_ylim (ymin, ymax)

    lines = LineCollection([zip(x, w) for w in ws.swapaxes(0, ax2)],
                           offsets=(xoffs, yoffs))

    #lines.set_array(x)
    G.fig.axes.add_collection(lines)
    text = "%s" % ws.info['name']
    if ws.info.has_key ("T"):
        text += ", T=%s K" % ws.info['T'][0]
    G.fig.axes.add_artist (mp.text.Text(xmin+(xmax-xmin)/12.0,
                                        ymax-(ymax-ymin)/12.0,
                                        text))
    G.fig.draw()


def decorate (can, wav):
    pass
