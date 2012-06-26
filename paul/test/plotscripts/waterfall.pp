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

def init(fig, view):
    log.debug ("Initializing waterfall plotscript")
    global G
    G = GlobalVars()
    G.gui = WaterfallSlicer(tparent=view)
    G.gui.sliced.connect (waterfall)
    G.fig = fig
    G.fig.reset()


def exit (fig, view):
    log.debug ("Exiting waterfall plotscript")
    global G
    del G.gui

def populate(fig, ws):
    global G
    G.fig = fig
    G.gui.slice(ws / ws.max())

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

    xmin = ws.axMin(ax1)
    xmax = ws.axMax(ax1)
    ymin = ws.min() + (yoffs*ws.shape[ax2]) * (yoffs<0)
    ymax = ws.max() + (yoffs*ws.shape[ax2]) * (yoffs>0)
    G.fig.axes.set_xlim (xmin, xmax)
    G.fig.axes.set_ylim (ymin, ymax)

    lines = LineCollection([zip(x, w) for w in ws.swapaxes(0, ax2)],
                           offsets=(xoffs, yoffs))

    #lines.set_array(x)
    G.fig.axes.add_collection(lines)
    G.fig.axes.add_artist (mp.text.Text(xmin+(xmax-xmin)/12.0,
                                        ymax-(ymax-ymin)/12.0,
                                        ws.info['name']))
    G.fig.draw()
    #print ws.info


def decorate (fig, waves):
    pass
