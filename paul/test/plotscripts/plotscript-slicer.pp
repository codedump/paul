#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

from PyQt4 import QtGui, QtCore
from paul.viewer.viewerwindow import ViewerWindow
from paul.base.wave import Wave

from paul.toolbox.slicers import SingleSlicer, WaterfallSlicer
from paul.toolbox.widgets import ValueWidget


class UiElements(QtCore.QObject):
    pass

class PlotElements:
    pass

GUI = UiElements()
PLOT = PlotElements()


@QtCore.pyqtSlot('double')
def setColor(i):
    ''' this one is called when the spin box color value changed '''
    global GUI, PLOT

    if not hasattr(PLOT, 'waves') or len(PLOT.waves) == 0:
        return

    data = PLOT.waves
    if not len(PLOT.canvas.axes.images) > 0:
        return
    img = PLOT.canvas.axes.images[0]

    log.debug ("values: %f, %f" % (GUI.col_min.value(), GUI.col_max.value()))

    dmin = data.min()
    dmax = data.max()

    #log.debug ("min/max: %f, %f" % (dmin, dmax))

    col_min = GUI.col_min.value() * (dmax-dmin) + dmin
    col_max = GUI.col_max.value() * (dmax-dmin) + dmin

    #log.debug ("absolute colors: %f, %f" % (col_min, col_max))

    img.set_clim (col_min, col_max)
    PLOT.canvas.draw()


@QtCore.pyqtSlot()
def addSlice():
    global GUI, PLOT
    GUI.slices.append(SingleSlicer(axis=0, viewer=ViewerWindow()))
    GUI.slices[-1].viewer.resize (400, 350)
    GUI.slices[-1].slice(wave=PLOT.waves[0])
    GUI.slices[-1].viewer.show()
     
@QtCore.pyqtSlot()
def addWaterfall():
    global GUI, PLOT
    GUI.slices.append(WaterfallSlicer(axis=0, viewer=ViewerWindow()))
    GUI.slices[-1].viewer.resize (500, 350)
    GUI.slices[-1].slice(wave=PLOT.waves[0])
    GUI.slices[-1].viewer.show()


def init(*args, **kwargs):
    '''
    Initialize the plotscript: this is a more complicated one. We are
    creating a fully function Qt GUI for this one, which will get
    displayed as a floating window, next to the Viewer :-)
    (In principle, it should be even possible to create this window
    as a Dock to the Viewer...)
    '''
    log.debug ("INIT")
    can = kwargs['can']
    win = kwargs['win']
    vars = kwargs['vars']
    global GUI, PLOT

    # some variables
    GUI.slices = []
    GUI.waves = []
    GUI.mainwin = win
    PLOT.canvas = can
    PLOT.waves = []

    # main window layout, resulting in:
    #   . ctrlwin: window containing parameter control widgets
    #   . plotwin: plotting window (ViewerWindow instance)

    # some setup variables (graph coloring etc)
    GUI.col_min = ValueWidget (None, "Color min: ", QtGui.QDoubleSpinBox(), -99, 99, 0.01, 0.00, setColor)
    GUI.col_max = ValueWidget (None, " max: ", QtGui.QDoubleSpinBox(), -99, 99, 0.01, 1.00, setColor)
    GUI.toolbar = QtGui.QToolBar()
    GUI.mainwin.addToolBarBreak()
    GUI.mainwin.addToolBar(GUI.toolbar)
    for w in GUI.col_min, GUI.col_max:
        GUI.toolbar.addWidget (w)
    GUI.toolbar.addAction ("+S", addSlice)
    GUI.toolbar.addAction ("+W", addWaterfall)


def exit(*args, **kwargs):
    log.debug ("EXIT")
    global GUI
    GUI.mainwin.removeToolBar (GUI.toolbar)
    for s in GUI.slices:
        log.debug ("Removing slicer %s" % s)
        #GUI.mainwin.removeToolBar (s.tools) # sometimes the __del__ function
        #                                    # is not called -- need to remove
        #                                    # the toolbar by hand.
        del s
    del GUI.slices
    del GUI
    return


def decorate(*args, **kwargs):
    '''
    Called when a newly plotted 2D image (or 1D graph) is to be decorated.
    '''
    can = kwargs['can']
    wav = kwargs['wav']
    global GUI, PLOT
    if not hasattr(can, 'axes'):
        log.debug ("No axes.")
        return
    if not len(wav) > 0:
        log.debug ("No waves.")
        return

    PLOT.waves = wav

    can.axes.set_xlabel (r'$k_{||}$ ($^\circ$)')
    can.axes.set_ylabel (r'E$_{total}$ (meV)')

    yax = can.axes.get_ylim()
    can.axes.set_ylim (min(yax), max(yax))
    #can.fig.colorbar(can.axes.images[0])

    for s in GUI.slices:
        if (s.viewer.isVisible()):
            log.debug ("Calling slice window %s (for axis %d)" % (s, s.slice_axis))
            s.slice(wave=wav[0])
        else:
            log.debug ("Removing hidden slice window %s (for axis %d)" % (s, s.slice_axis))
            GUI.slices.remove(s)


    
