#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

from PyQt4 import QtGui, QtCore
from paul.viewer.viewerwindow import ViewerWindow

class UiElements:
    pass

class PlotElements:
    pass

GUI = UiElements()
PLOT = PlotElements()

@QtCore.pyqtSlot('double')
def colValueChanged(i):
    ''' this one is called when the spin box color value changed '''
    global GUI, PLOT

    if not hasattr(PLOT, 'waves') or len(PLOT.waves) == 0:
        return

    data = PLOT.waves[0]
    if not len(PLOT.canvas.axes.images) > 0:
        return
    img = PLOT.canvas.axes.images[0]

    dmin = data.min()
    dmax = data.max()

    col_min = GUI.col_min.value() * (dmax-dmin) + dmin
    col_max = GUI.col_max.value() * (dmax-dmin) + dmin
    img.set_clim (col_min, col_max)
    PLOT.canvas.draw()

@QtCore.pyqtSlot('int')
def slicingChanged(i):
    ''' Called when the slicing parameters changes.'''
    global GUI, PLOT
    if not hasattr(PLOT, 'waves') or len(PLOT.waves) == 0:
        return
    data = PLOT.waves[0]
    x_from = GUI.slice_from.value()
    x_to = x_from + GUI.slice_delta.value()
    PLOT.cur_slice = data[x_from]
    GUI.slicewin.plotWaves([PLOT.cur_slice])
    #GUI.slicewin.refresh()
    

def init(canvas):
    '''
    Initialize the plotscript: this is a more complicated one. We are
    creating a fully function Qt GUI for this one, which will get
    displayed as a floating window, next to the Viewer :-)
    (In principle, it should be even possible to create this window
    as a Dock to the Viewer...)
    '''
    log.debug ("INIT")
    global GUI, PLOT

    # main window layout, resulting in:
    #   . ctrlwin: window containing parameter control widgets
    #   . plotwin: plotting window (ViewerWindow instance)
    GUI.mainwin = QtGui.QWidget()
    GUI.hbox = QtGui.QHBoxLayout(GUI.mainwin)
    GUI.splitter = QtGui.QSplitter(GUI.mainwin)
    GUI.hbox.setContentsMargins (0, 0, 0, 0)
    GUI.hbox.addWidget (GUI.splitter)    
    GUI.ctrlwin = QtGui.QWidget(GUI.splitter)
    GUI.ctrlbox = QtGui.QGridLayout(GUI.ctrlwin)
    GUI.ctrlwin.setLayout (GUI.ctrlbox)
    GUI.mainwin.setWindowTitle ("2D Slicer")

    # some setup variables (graph coloring etc)
    GUI.col_min = QtGui.QDoubleSpinBox (GUI.ctrlwin)
    GUI.col_max = QtGui.QDoubleSpinBox (GUI.ctrlwin)
    GUI.col_max.setValue (1)
    row = 0
    for w in GUI.col_min, GUI.col_max:
        w.valueChanged.connect(colValueChanged)
        w.setSingleStep (0.01)
        w.setRange (-99, 99)
        GUI.ctrlbox.addWidget (w, row, 1)
        row += 1

    # slice display: we'll create a ViewerWindow() instance here
    # -- after all, that's what it is made for :-)
    # however, if 'Automagic' plot script will be loaded, then
    # there exists the risk that we will end up in a recursion
    # (e.g. if this plotscript is in the Automagic path).
    # We need to break recurtion by starting up with 'none' plotscript.
    GUI.slicewin = ViewerWindow(plotscript='(none)')
    GUI.slicewin.setParent (GUI.splitter)

    GUI.slice_from  = QtGui.QSpinBox (GUI.ctrlwin)
    GUI.slice_delta = QtGui.QSpinBox (GUI.ctrlwin)
    for w in GUI.slice_from, GUI.slice_delta:
        w.valueChanged.connect (slicingChanged)

    # layout stuff
    GUI.ctrlbox.addWidget (QtGui.QLabel ("Color min: ", GUI.ctrlwin), 0, 0)
    GUI.ctrlbox.addWidget (GUI.col_min, 0, 1)
    GUI.ctrlbox.addWidget (QtGui.QLabel ("Color max: ", GUI.ctrlwin), 1, 0)
    GUI.ctrlbox.addWidget (GUI.col_max, 1, 1)
    GUI.ctrlbox.addWidget (QtGui.QLabel ("Slice from: ", GUI.ctrlwin), 2, 0)
    GUI.ctrlbox.addWidget (GUI.slice_from, 2, 1)
    GUI.ctrlbox.addWidget (QtGui.QLabel ("Slice delta: ", GUI.ctrlwin), 3, 0)
    GUI.ctrlbox.addWidget (GUI.slice_delta, 3, 1)

    GUI.mainwin.show()

    PLOT.canvas = canvas
    PLOT.waves = []


def exit(canvas):
    log.debug ("EXIT")
    global GUI
    del GUI
    return


def decorate(can, wav):
    '''
    Called when a newly plotted 2D image (or 1D graph) is to be decorated.
    '''
    global GUI, PLOT
    if not hasattr(can, 'axes'):
        log.debug ("No axes.")
        return
    if not len(wav) > 0:
        log.debug ("No waves.")
        return
    if not GUI.mainwin.isVisible():
        GUI.mainwin.show()

    PLOT.waves = wav

    for w in GUI.slice_from, GUI.slice_delta:
        w.setRange (0, wav[0].shape[0])
    slicingChanged(-1)

    can.axes.set_xlabel (r'$k_{||}$ (deg.)')
    can.axes.set_ylabel (r'E$_{total}$ (meV)')

    yax = can.axes.get_ylim()
    can.axes.set_ylim (min(yax), max(yax))
    #can.fig.colorbar(can.axes.images[0])

    log.debug ("Displayed.")
    
