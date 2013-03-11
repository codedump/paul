#!/usr/bin/python

from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
import paul.base.wave as wave
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
    
    if isinstance(wlist, wave.Wave):
        return imwater(fig_ax, wlist, axis, offs, xlim, ylim, autoscale=True)

    xlist = [np.arange(start=wave.WCast(w).dim[0].offset,
                       stop=wave.WCast(w).dim[0].end,
                       step=wave.WCast(w).dim[0].delta) for w in wlist]
    
    if xlim == (0, 0):
        xlim = (min([wave.WCast(w).dim[0].min for w in wlist]) - (offs[0]*len(wlist))*(offs[0]<0),
                max([wave.WCast(w).dim[0].max for w in wlist]) + (offs[0]*len(wlist))*(offs[0]>0))
    if ylim == (0, 0):
        ylim = (min([np.nanmin(w) for w in wlist]) + (offs[1]*len(wlist)) * (offs[1]<0),
                max([np.nanmin(w) for w in wlist]) + (offs[1]*len(wlist)) * (offs[1]>0))

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

    if len(ignore) > 0:
        log.error ("NOT IMPLEMENTED: 'ignore' parameter!  "
                   "(Hint: You should rather preselect input wave range.)")


    if not isinstance(wlist, wave.Wave):
        new_wlist = np.vstack (wlist)  # this only works if waves have
                                       # the same number of points
        wlist = new_wlist

    #wlist.swapaxes(0, axis)

    if len(wlist.shape) != 2:
        err = "Only 2D waves can be plotted as waterfalls (dim = %d here)" % len(wlist.shape)
        log.error (err)
        raise ValueError (err)

    ## first, calculate the X axis values
    x = np.arange(start = wlist.dim[1].offset,
                  stop  = wlist.dim[1].end,
                  step  = wlist.dim[1].delta)

    
    ## Do the auto-scaling magic... :-)
    
    # some indices along the y axis...
    min_i  = wlist.ax(0).x2i_rnd(wlist.dim[0].min) # index of k|| = min
    zero_i = wlist.ax(0).x2i_rnd(0)                # index of k|| = 0
    max_i  = wlist.ax(0).x2i_rnd(wlist.dim[0].max) # index of k|| = max
    
    # the "natural" y scaling, that would be in effect if we didn't change anything
    ylim_wouldbe = (0, wlist.shape[0]*abs(offs[1]))

    # the original y-range of the data
    ylim_data    = (wlist.dim[0].min, wlist.dim[0].max)


    # zoom factor from "original" to "natural" scaling
    axzoom       = (ylim_data[1]-ylim_data[0]) / (ylim_wouldbe[1]-ylim_wouldbe[0])

    # recalculate the offs parameter to have the same effect on our
    # "rescaled" data that the original parameter would have had
    # on the original data.
    offs         = (offs[0], offs[1]*axzoom)

    # Usually, LineCollection would start displaying
    # data at y=0. Here we calculate an offset we need to
    # apply to the data in order to have the lines aligned
    # at the correct position on the y axis.
    data_yshift  = -(zero_i - min_i) * offs[1]
                                            

    # Scale the wave intensity to have it appear as big
    # (compared to the 'offs' parameter) as it would have been
    # if we didn't interfere. Then apply the y-shifting
    wlist *= axzoom
    wlist += data_yshift

    ## set the proper limis (only if user didn't specify his own).
    if xlim == (0, 0):
        xlim = (wlist.dim[1].min - (offs[0]*len(wlist))*(offs[0]<0),
                wlist.dim[1].max + (offs[0]*len(wlist))*(offs[0]>0))
    if ylim == (0, 0):
        ylim = (ylim_data[0], ylim_data[1])

    ## ...then go for the actual work.
    lines = LineCollection([zip(x, w) for w in wlist], offsets=offs)

    if xlim is not None:
        fig_ax.set_xlim (xlim)
    if ylim is not None:
        fig_ax.set_ylim (ylim[0]-(ylim[1]-ylim[0])*0.05,
                         ylim[1]+(ylim[1]-ylim[0])*0.05)

    fig_ax.add_collection (lines)
    return lines


