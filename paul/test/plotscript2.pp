#!/usr/bin/python

import logging
log = logging.getLogger (__name__)
set_win = None

from PyQt4 import QtGui

def init(canvas):
    log.debug ("INIT")
    global set_win
    set_win = QtGui.QLabel ("Look ma', I can create Qt Widgets!")
    return

def exit(canvas):
    log.debug ("EXIT")
    global set_win
    del set_win
    return

def decorate(can, wav):
    '''
    Called when a newly plotted 2D image (or 1D graph) is to be decorated.
    '''
    if not hasattr(can, 'axes'):
        return

    if not len(wav) > 0:
        return

    global set_win
    set_win.show()

    data = wav[0]

    can.axes.set_xlabel (r'$k_{||}$ (deg.)')
    can.axes.set_ylabel (r'E$_{total}$ (meV)')

    yax = can.axes.get_ylim()
    can.axes.set_ylim (min(yax), max(yax))

    log.debug ("Image values: \t\t%f, %f" % (data.min(), data.max()))
    dmin = data.min()
    dmax = data.max()

    if not len(can.axes.images) > 0:
        return

    img = can.axes.images[0]

    col_min = 0.00 * (dmax-dmin) + dmin
    col_max = 0.59 * (dmax-dmin) + dmin
    log.debug ("Color scale values: \t%f, %f" % (col_min, col_max))
    
    
    img.set_clim (col_min, col_max)
    #can.fig.colorbar(img)
    
