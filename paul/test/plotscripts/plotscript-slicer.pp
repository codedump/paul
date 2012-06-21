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

# simple container to hold a label and a SpinBox
class ValueWidget(QtGui.QWidget):
    def __init__ (self, parent=None, label='value', spin=None, vmin=0, vmax=99, vinc=1, vdef=0, slot=None):
        QtGui.QWidget.__init__ (self, parent)
        self.label = QtGui.QLabel (label, self)
        if spin is None:
            spin = QtGui.QSpinBox(self)
        self.spin = spin
        self.spin.setRange (vmin, vmax)
        self.spin.setSingleStep (vinc)
        self.spin.setValue (vdef)
        if slot is not None:
            self.spin.valueChanged.connect (slot)
        self.box = QtGui.QHBoxLayout(self)
        self.box.addWidget (self.label)
        self.box.addWidget (self.spin)
    
    def value(self):
        return self.spin.value()


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

    log.debug ("values: %f, %f" % (GUI.col_min.value(), GUI.col_max.value()))

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
    delta = x_to - x_from
    if delta <= 0:
        delta = 1
    PLOT.cur_slice = data[x_from:x_to,0:].sum(0) / (x_to-x_from)
    GUI.slicewin.plotWaves([PLOT.cur_slice])
    #GUI.slicewin.refresh()
    

def init(canvas, mainwin):
    '''
    Initialize the plotscript: this is a more complicated one. We are
    creating a fully function Qt GUI for this one, which will get
    displayed as a floating window, next to the Viewer :-)
    (In principle, it should be even possible to create this window
    as a Dock to the Viewer...)
    '''
    log.debug ("INIT")
    global GUI, PLOT

    GUI.viewer = mainwin

    # main window layout, resulting in:
    #   . ctrlwin: window containing parameter control widgets
    #   . plotwin: plotting window (ViewerWindow instance)

    # some setup variables (graph coloring etc)
    GUI.col_min = ValueWidget (None, "Color min: ", QtGui.QDoubleSpinBox(), -99, 99, 0.01, 0.00, colValueChanged)
    GUI.col_max = ValueWidget (None, "Color max: ", QtGui.QDoubleSpinBox(), -99, 99, 0.01, 1.00, colValueChanged)
    GUI.toolbar = QtGui.QToolBar()
    mainwin.addToolBar(QtCore.Qt.RightToolBarArea, GUI.toolbar)

    GUI.slice_from  = ValueWidget (None, "X-from: ", QtGui.QSpinBox(), 0, 99, 1, 0, slicingChanged)
    GUI.slice_delta = ValueWidget (None, "X-to: ", QtGui.QSpinBox(), 0, 99, 1, 0, slicingChanged)
    for w in GUI.col_min, GUI.col_max, GUI.slice_from, GUI.slice_delta:
        GUI.toolbar.addWidget (w)

    GUI.slicewin = ViewerWindow()
    GUI.slicewin.show()

    PLOT.canvas = canvas
    PLOT.waves = []


def exit(canvas, mainwin):
    log.debug ("EXIT")
    GUI.viewer.removeToolBar (GUI.toolbar)
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

    PLOT.waves = wav

    for w in GUI.slice_from, GUI.slice_delta:
        w.spin.setRange (0, wav[0].shape[0])
    slicingChanged(-1)

    can.axes.set_xlabel (r'$k_{||}$ ($^\circ$)')
    can.axes.set_ylabel (r'E$_{total}$ (meV)')

    yax = can.axes.get_ylim()
    can.axes.set_ylim (min(yax), max(yax))
    #can.fig.colorbar(can.axes.images[0])

    log.debug ("Displayed.")
    
