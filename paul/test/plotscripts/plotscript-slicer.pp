#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

from PyQt4 import QtGui, QtCore
from paul.viewer.viewerwindow import ViewerWindow
from paul.base.wave import Wave


class UiElements(QtCore.QObject):
    pass

class PlotElements:
    pass

GUI = UiElements()
PLOT = PlotElements()

class ValueWidget(QtGui.QWidget):
    '''
    Simple container to hold a label and a SpinBox.
    Used for convenience.
    '''
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
        self.box.setContentsMargins (0, 0, 0, 0)
        self.box.setSpacing (0)
        self.box.addWidget (self.label)
        self.box.addWidget (self.spin)
    
    def value(self):
        return self.spin.value()


class WaveSlice (ViewerWindow):
    '''
    Represents an (N-1)-dim slice of an N-dimensional wave.
    This class holds both the GUI elements (most prominent
    being the ViewerWindow used to display the slice, and
    a tool bar used to control the slicing parameters) aswell
    as the functionality (i.e. _creating_ the slice in the
    first place).
    '''

    slicing = QtCore.pyqtSignal ([int, int, int])

    def __init__(self, parent=None, axis=0, master=None, tp=None):
        ViewerWindow.__init__(self, parent)
        
        if parent is not None:
            self.setParent (parent)

        self.master_wave = master  # reference to the wave we will be slicing
        self.slice_wave  = None
        self.slice_axis  = axis
        
        # initialize GUI
        self.tools = QtGui.QToolBar()
        self.val_axis =  ValueWidget (None, "Axis: ",  QtGui.QSpinBox(), 0, 99, 1, 0, self.slice)
        self.val_from  = ValueWidget (None, "From: ",  QtGui.QSpinBox(), 0, 99, 1, 0, self.slice)
        self.val_delta = ValueWidget (None, "Delta: ", QtGui.QSpinBox(), 1, 99, 1, 0, self.slice)
        self.slice_label = QtGui.QLabel ("<info>")

        # attach toolbar to parent, by default (or to us, if no parent)
        self.tools_parent = self
        if tp is not None and hasattr(tp, 'addToolBar'):
            self.tools_parent = tp
        self.tools_parent.addToolBarBreak()
        self.tools_parent.addToolBar(self.tools)
        for w in [ self.val_axis, self.val_from, self.val_delta, self.slice_label ]:
            self.tools.addWidget(w)


    def __del__(self):
        # need to remove the toolbar by hand (because it might be attributed to
        # the parent window...)
        log.debug ("Deleting slice for axis %d" % self.slice_axis)
        self.tools_parent.removeToolBar (self.tools)



    @QtCore.pyqtSlot('int')
    def slice (self, vfrom=0, vdel=0, wave=None):
        '''
        Called when the slicing parameters or the master wave change.
        '''
        #log.debug ("Slicing along axis %d" % self.slice_axis)

        if wave is not None:
            # set a new master wave and adjust spinbox limits
            self.master_wave = wave
        if self.master_wave is None:
            return
        if self.slice_axis != self.val_axis.value():
            self.slice_axis = self.val_axis.value()

        # need to re-set the range (mostly because slice_axis or master_wave changed
        self.val_from.spin.setRange (0, self.master_wave.shape[self.slice_axis])
        self.val_delta.spin.setRange (1, self.master_wave.shape[self.slice_axis])
        self.val_axis.spin.setRange (0, self.master_wave.ndim-1)

        xfrom = self.val_from.value()
        delta = self.val_delta.value()
        xto = xfrom + delta
        log.debug ("Slicing axis %d, from %d to %d" % (self.slice_axis, xfrom, xto))
        
        # we want to slice along an arbitrary axis, here's the strategy:
        #   . swap the axis 0 with slice_axis
        #   . slice along axis 0
        #   . swap back (slice_axis with 0)
        #   . sum along slice_axis (we want to integrate out the
        #     dimension we slice along...)
        self.slice_wave = self.master_wave.swapaxes(0,self.slice_axis)[xfrom:xto].swapaxes(self.slice_axis,0).sum(self.slice_axis)
        self.slice_label.setText(" ax%d: %5.2f : %5.2f"
                                 % (self.slice_axis,
                                    self.slice_wave.i2f(self.slice_axis, xfrom),
                                    self.slice_wave.i2f(self.slice_axis, xto)))
        self.plotWaves (self.slice_wave)


@QtCore.pyqtSlot('double')
def setColor(i):
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


@QtCore.pyqtSlot()
def addSlice():
    global GUI, PLOT
    GUI.slices.append(WaveSlice(axis=0))
    GUI.slices[-1].resize (400, 350)
    GUI.slices[-1].slice(wave=PLOT.waves)
    GUI.slices[-1].show()
     

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

    # some variables
    GUI.slices = []
    GUI.waves = []
    GUI.mainwin = mainwin
    PLOT.canvas = canvas
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
    GUI.toolbar.addAction ("+Slice", addSlice)


def exit(canvas, mainwin):
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

    can.axes.set_xlabel (r'$k_{||}$ ($^\circ$)')
    can.axes.set_ylabel (r'E$_{total}$ (meV)')

    yax = can.axes.get_ylim()
    can.axes.set_ylim (min(yax), max(yax))
    #can.fig.colorbar(can.axes.images[0])

    for s in GUI.slices:
        if (s.isVisible()):
            log.debug ("Calling slice window %s (for axis %d)" % (s, s.slice_axis))
            s.slice (wave=wav)
        else:
            log.debug ("Removing hidden slice window %s (for axis %d)" % (s, s.slice_axis))
            GUI.slices.remove(s)


    
