#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

from PyQt4 import QtGui, QtCore

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
    log.debug ("Color scale values: \t%f, %f (%f, %f)" %
               (col_min, col_max, GUI.col_min.value(), GUI.col_max.value()))
    img.set_clim (col_min, col_max)
    PLOT.canvas.draw()
    

def init(canvas):
    log.debug ("INIT")
    global GUI
    GUI.win = QtGui.QWidget()
    GUI.win.setWindowTitle ("Volatile Plot Settings")
    GUI.layout = QtGui.QGridLayout(GUI.win)
    GUI.win.setLayout (GUI.layout)

    GUI.col_min = QtGui.QDoubleSpinBox (GUI.win)
    GUI.col_max = QtGui.QDoubleSpinBox (GUI.win)
    GUI.col_max.setValue (1)
    row = 0
    for w in GUI.col_min, GUI.col_max:
        w.valueChanged.connect(colValueChanged)
        w.setSingleStep (0.01)
        GUI.layout.addWidget (w, row, 1)
        row += 1


    GUI.layout.addWidget (QtGui.QLabel ("Color min:", GUI.win), 0, 0)
    GUI.layout.addWidget (GUI.col_min, 0, 1)
    GUI.layout.addWidget (QtGui.QLabel ("Color max:", GUI.win), 1, 0)
    GUI.layout.addWidget (GUI.col_max, 1, 1)

    return

def exit(canvas):
    log.debug ("EXIT")
    global GUI
    del GUI
    return

def decorate(can, wav):
    '''
    Called when a newly plotted 2D image (or 1D graph) is to be decorated.
    '''
    if not hasattr(can, 'axes'):
        return

    if not len(wav) > 0:
        return

    global GUI, PLOT
    GUI.win.show()
    PLOT.canvas = can
    PLOT.waves = wav

    can.axes.set_xlabel (r'$k_{||}$ (deg.)')
    can.axes.set_ylabel (r'E$_{total}$ (meV)')

    yax = can.axes.get_ylim()
    can.axes.set_ylim (min(yax), max(yax))

    #can.fig.colorbar(img)
    
