#!/usr/bin/python

from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from paul.base.wave import Wave
import numpy as np
from pprint import pprint

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


def plotwater (fig_ax, wlist, axis=0, offs=(0, 0), xlim=(0,0), ylim=(0,0)):
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
    
    if isinstance(wlist, Wave):
        return imwater(fig_ax, wlist, axis, offs, xlim, ylim, autoscale=True)

    xlist = [np.arange(start=w.axOff(), stop=w.axEnd(), step=w.axDelta()) for w in wlist]
    
    if xlim == (0, 0):
        xlim = (min([w.axMin(0) for w in wlist]) - (offs[0]*len(wlist))*(offs[0]<0),
                max([w.axMax(0) for w in wlist]) + (offs[0]*len(wlist))*(offs[0]>0))
    if ylim == (0, 0):
        ylim = (min([w.min() for w in wlist]) + (offs[1]*len(wlist)) * (offs[1]<0),
                max([w.max() for w in wlist]) + (offs[1]*len(wlist)) * (offs[1]>0))

    lines = LineCollection([zip(x, w) for x, w in zip(xlist,wlist)], offsets=offs)

    if xlim is not None:
        fig_ax.set_xlim (xlim)
    if ylim is not None:
        fig_ax.set_ylim (ylim)

    fig_ax.add_collection (lines)
    return lines

plot_water = plotwater   # define an alias, for interface compatibility



def imwater (fig_ax, wlist, axis=0, offs=(0, 0), xlim=(0,0), ylim=(0,0), autoscale=True,
             ignore=[]):
    '''
    Same as plotwater(), but designed to work for 2D waves.
    if *autoscale* is True, then resulting line collection will have
    the same y-scaling as the original 2D wave.
    
    For a 2D wave, the X and Y axis limits are retained (i.e. a 2D
    image plot or a waterfall plot will have the same ranges).
    This is done by re-interpreting the *offs* parameter. The
    algorithm is roughly the following:
       1) Calculate the would-be y-axis span (which is 
           N * offs[y], where N is the number if 1D slices)
       2) Calculate the real y-axis span of the 2 wave
       3) Calculate the ratio between spans of (1) and (2)
       4) Scale *offs*[y] and the signal intensity by the ratio at (3).

    If specified, *ignore* is a list of line indices (in the final 'waterfall'
    line set) which to exclude from plotting.
    '''
    if len(wlist) == 0:
        return


    if not isinstance(wlist, Wave):
        new_wlist = np.vstack (wlist)  # this only works if waves have
                                       # the same number of points
        wlist = new_wlist
        

    if len(wlist.shape) != 2:
        err = "Only 2D waves can be plotted as waterfalls (dim = %d here)" % len(wlist.shape)
        log.error (err)
        raise ValueError (err)

    min_i  = wlist.ax(0).rx2i(wlist.axMin()) # index of k|| = min
    zero_i = wlist.ax(0).rx2i(0)             # index of k|| = 0
    max_i  = wlist.ax(0).rx2i(wlist.axMax()) # index of k|| = max
    x = np.arange(start=wlist.axOff(1),
                  stop=wlist.axEnd(1),
                  step=wlist.axDelta(1))
    ylim_wouldbe = (0, wlist.shape[0]*offs[1])
    ylim_data = (wlist.axMin(0), wlist.axMax(0))
    axzoom = (ylim_data[1]-ylim_data[0]) / (ylim_wouldbe[1]-ylim_wouldbe[0])
    offs = (offs[0], offs[1]*axzoom)
    data_offset = -(zero_i - min_i) * offs[1]
    
    wlist *= axzoom
    wlist += data_offset
    
    if xlim == (0, 0):
        xlim = (wlist.axMin(1) - (offs[0]*len(wlist))*(offs[0]<0),
                wlist.axMax(1) + (offs[0]*len(wlist))*(offs[0]>0))
    if ylim == (0, 0):
        ylim = (ylim_data[0], ylim_data[1])

    lines = LineCollection([zip(x, w) for w in wlist], offsets=offs)

    if xlim is not None:
        fig_ax.set_xlim (xlim)
    if ylim is not None:
        fig_ax.set_ylim (ylim[0]-(ylim[1]-ylim[0])*0.05,
                         ylim[1]+(ylim[1]-ylim[0])*0.05)

    fig_ax.add_collection (lines)
    return lines


