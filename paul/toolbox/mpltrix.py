#!/usr/bin/python

from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
import numpy as np

import logging
log = logging.getLogger(__name__)

'''
Plot tricks for matplotlib plotting
'''


def gridplot (grid, loc, rowspan=1, colspan=1):
    '''
    Returns a matplotlib.gridspec.SubplotSpec for a subplot.
    The resulting object can then be added to a matplotlib.figure
    using the add_subplot() method.
    '''
    gridspec = GridSpec (grid[0], grid[1])
    subplotspec = gridspec.new_subplotspec(loc, rowspan, colspan)
    return subplotspec


def plot_water (fig_ax, wlist, offs=(0, 0), xlim=(0,0), ylim=(0,0)):
    '''
    Creates a waterfall plot on the matplotlib Axes instance
    *fig_ax* from a list of 1D waves specified by *wlist*.
    With the axis limits *xlim* and *ylim*, or auto-calculated 
    axes limits if none are specified.
    (Auto-calculated limits are quite good BTW, maybe sometimes slightly
    too large along the y-axis, depending on the data;
    mostly useful and always encompassing the whole area.)
    
    Returns a LineCollection object, which can be used to further
    manipulate the waterfall plot appearance.
    '''
    if len(wlist) == 0:
        return

    x = np.arange(start=wlist[0].axOff(0),
                  stop=wlist[0].axEnd(0),
                  step=wlist[0].axDelta(0))

    if xlim == (0, 0):
        xlim = (min([w.axMin(0) for w in wlist]) - (offs[0]*len(wlist))*(offs[0]<0),
                max([w.axMax(0) for w in wlist]) + (offs[0]*len(wlist))*(offs[0]>0))

    if ylim == (0, 0):
        ylim = (min([w.min() for w in wlist]) + (offs[1]*len(wlist)) * (offs[1]<0),
                max([w.max() for w in wlist]) + (offs[1]*len(wlist)) * (offs[1]>0))

    lines = LineCollection([zip(x, w) for w in wlist], offsets=offs)

    fig_ax.set_xlim (xlim)
    fig_ax.set_ylim (ylim)
    fig_ax.add_collection (lines)
    return lines


