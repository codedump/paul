#!/usr/bin/python

from matplotlib.collections import LineCollection
import numpy as np
import pprint as pp
from paul.toolbox.slicers import WaterfallSlicer

import logging
log = logging.getLogger (__name__)

class GlobalVars:
    pass

def init(fig, view):
    log.debug ("Initializing waterfall plotscript")
    global G
    G = GlobalVars()
    G.gui = WaterfallSlicer(tparent=view)
    G.gui.sliced.connect (waterfall)
    G.fig = fig

def exit (fig, view):
    log.debug ("Exiting waterfall plotscript")
    global G
    del G.gui

def populate(fig, ws):
    global G
    G.fig = fig
    G.gui.slice(ws)

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
    x = np.arange(start=ws.axOffset(ax1),
                  stop=ws.axEndpoint(ax1),
                  step=ws.axDelta(ax1))

    log.debug ("Axis limits: %f, %f" % (ws.axMin(ax1), ws.axMax(ax1)))

    G.fig.axes.set_xlim (ws.axMin(ax1), ws.axMax(ax1))
    G.fig.axes.set_ylim (ws.min() + (yoffs*ws.shape[ax2]) * (yoffs<0),
                         ws.max() + (yoffs*ws.shape[ax2]) * (yoffs>0))

    lines = LineCollection([zip(x, w) for w in ws.swapaxes(0, ax2)],
                           offsets=(xoffs, yoffs))

    lines.set_array(x)
    G.fig.axes.add_collection(lines)
    G.fig.draw()


def decorate (fig, waves):
    pass
