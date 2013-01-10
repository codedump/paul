#!/usr/bin/python

from PyQt4 import QtGui, QtCore

from paul.viewer.viewerwindow import ViewerWindow
from paul.toolbox.widgets import ValueWidget
from paul.base.wave import Wave
from paul.toolbox.atrix import ncomp

import logging
log = logging.getLogger (__name__)


'''
Some base classes for wave slicing in plotscripts.
'''

class BaseSlicer(QtCore.QObject):
    '''
    Base class for a slicer. SingleSlice and MultiSlice will be derived
    from this class.
    This class implements basic UI functionality, leaving for
    the sub-classes only the UI-initialization and the
    actual slicing to be done.
    For now, it GUI and functionality are merged (maybe we should
    split them?)
    Typically, slicing will be controlled by a number of values
    specified in spin boxes, to be created by the respective
    sub-classes. In order to have a snappy UI even if the
    slicing/plotting job is more work intensive, we don't trigger
    the slicing on every spin box update. Instead, on spin box updates,
    we trigger a single-shot timer with a short timeout value (~0.1 sec).
    On timeout, we do the slicing.
    This way, on single spin-box changes, the response is quick.
    On multiple changes (like auto-repeating keystrokes), only
    the last change will affect the slicing, saving lots of time.
    For this to work, it is crucial that spin boxes don't bind to
    slice(), but to triggerSlicing() instead
    '''

    # Signal emitted when a slicing has occured.
    # Typically, it is used by external viewers to display
    # the slice_wave.

    sliced = QtCore.pyqtSignal(list)

    def __init__(self, parent=None, master=None, tparent=None, viewer=None):
        QtCore.QObject.__init__(self, parent)
        '''
        Parameters:
           parent:  The plotting widget parent (where the plotting widget will be shown)
           master:  Master wave (can be changed later)
           tparent: Toolbar parent, i.e. to which QMainWindow to attach
                    the toolbar to.
           viewer:  The viewer in which to show the slices.
                    If none is specified, only a signal will be emitted.
        '''
        if viewer is not None:
            self.viewer = viewer
            self.sliced.connect (self.viewer.plotWaves)

        if parent is not None:
            if hasattr(self, 'viewer'):
                self.viewer.setParent (parent)

        self.master_wave = master  # reference to the wave we will be slicing
        self.slice_wave  = None    # the result (i.e. wave to be plotted)
        self.tools       = QtGui.QToolBar("Slicing specs")   # subclasses should attach all UI elements here

        # the slicing timer: will be activated by triggerSlicing(),
        # will call slice() on timeout.
        self.timer = QtCore.QTimer()
        self.timer.setSingleShot (True)
        self.timer.timeout.connect(self.slice)
        
        # attach toolbar to parent, by default (or to us, if no parent)
        if tparent is not None:
            self.tools_parent = tparent
            log.debug ("Toolbar parent: explicit, %s" % tparent)
        elif hasattr(self, 'viewer'):
            self.tools_parent = self.viewer
            log.debug ("Toolbar parent: implicit, viewer %s" % self.viewer)
        else:
            self.tools_parent = QtGui.QMainWindow()
            log.debug ("Toolbar parent: implicit, new window")
        self.tools_parent.addToolBarBreak()
        self.tools_parent.addToolBar(self.tools)


    def __del__(self):
        '''
        Need to remove the toolbar by hand because it might be attributed to
        a different parent window.
        '''
        log.debug ("Deleting slice for axis %d" % self.slice_axis)
        self.tools_parent.removeToolBar (self.tools)


    @QtCore.pyqtSlot()
    @QtCore.pyqtSlot(int)
    def triggerSlicing(self, val=0):
        '''
        The slicing timer. The timeout() signal is connected to self.slice().
        '''
        self.timer.start (100)  # wait 100 ms for repeated key input

        
    def _prepareMaster (self, wave):
        '''
        To be called internally, from self.slice(), when the master
        wave changes. Typically it would be used to do setup work 
        associated with the new wave (i.e. set new spinbox values...)
        '''
        self.master_wave = wave        
    Base_prepareMaster = _prepareMaster


    @QtCore.pyqtSlot()      
    def slice(self, wave=None):
        '''
        The function that does the actual slicing. Subclasses
        need to implement this.
        '''
        pass


class AxisSlicer(BaseSlicer):
    '''
    Base class for an axis-based slicer. This one typically
    slices perpendicular to a specified wave axis.
    (Opposed to a free slicer, where the the slicing takes
    place not perpendicular to a particular axis, but along
    a freely specified line.)
    '''

    def __init__ (self, parent, axis=0, master=None, tparent=None, viewer=None):
        '''
        Parameters are those of BaseSlicer and the following:
          axis:  Axis perpendicular to which to slice
        '''
        BaseSlicer.__init__(self, parent, master, tparent, viewer)
        self.slice_axis = axis  # axis perpendicular to which we will slice
        self.val_axis =  ValueWidget (None, "Axis: ",  QtGui.QSpinBox(),
                                      0, 99, 1, axis, self._prepareAxis)
        self.tools.addWidget(self.val_axis)


    @QtCore.pyqtSlot (int)
    def _prepareAxis (self, axis):
        ''' Make sure val_axis and slice_axis are in sync. '''
        self.slice_axis = self.val_axis.value()
        self._prepareMaster(self.master_wave)  # adjust spinbox limits
        self.slice(self.master_wave)
    Axis_prepareAxis = _prepareAxis


    def _prepareMaster(self, wave):
        self.Base_prepareMaster(wave)
        if self.master_wave is None:
            return
        self.val_axis.spin.setRange (0, self.master_wave.ndim-1)
    Axis_prepareMaster = _prepareMaster



class SingleSlicer (AxisSlicer):
    '''
    Represents an (N-1)-dim slice of an N-dimensional wave.
    '''

    # need to reimplement some methods, but we wish to retain
    # their functionality:

    def __init__(self, parent=None, axis=0, master=None, tparent=None, viewer=None):
        AxisSlicer.__init__(self, parent, axis, master, tparent, viewer)

        # initialize GUI
        self.val_from  = ValueWidget (None, "From: ",  QtGui.QSpinBox(), 0, 99, 1, 0,    self.triggerSlicing)
        self.val_delta = ValueWidget (None, "Delta: ", QtGui.QSpinBox(), 1, 99, 1, 1,    self.triggerSlicing)
        self.slice_label = QtGui.QLabel ("<info>")

        # attach toolbar to parent, by default (or to us, if no parent)
        for w in [ self.val_from, self.val_delta, self.slice_label ]:
            self.tools.addWidget(w)


    def _prepareMaster(self, wave):
        self.Axis_prepareMaster (wave)
        if self.master_wave is None:
            return
        self.val_from.spin.setRange (0, self.master_wave.shape[self.slice_axis])
        self.val_delta.spin.setRange (1, self.master_wave.shape[self.slice_axis])


    @QtCore.pyqtSlot()
    def slice (self, wave=None):
        # everything we do depends strongly on the selected slice_axis,
        # so if that one changes, we can aswell trigger a _prepareMaster()
        # session.
        if wave is not None:
            self._prepareMaster(wave)
        if self.master_wave is None:
            return

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
        log.debug ("%d axes: %s" % (self.slice_wave.ndim, [str(a) for a in self.slice_wave.dim]))
        self.slice_label.setText(" ax%d: %5.2f : %5.2f"
                                 % (self.slice_axis,
                                    self.master_wave.dim[self.slice_axis].i2x (xfrom),
                                    self.master_wave.dim[self.slice_axis].i2x (xto)))
        self.sliced.emit ([self.slice_wave])


class CompressingSlicer (AxisSlicer):
    '''
    This class generates a compressed version of the specified wave.
    It is useful for creating waterfall diagrams along specified axes.
      axis:  The slice axis
      step:  Distance between display slices
      norm:  Normalize: if specified, the wave will be normalized
             to the maximum intensity point prior to intensity
             manipulation.
      intg:  Specifies whether slices between step N and N+step
             should be added up.
    '''
    

    def __init__ (self, parent=None, axis=0, master=None, tparent=None, viewer=None):
        AxisSlicer.__init__(self, parent, axis, master, tparent, viewer)

        self.old_intg = 1 # need these to link intg/step 
        self.old_step = 1

        self.xoffset = 0  # property that will be needed at plot time
        self.yoffset = 0 

        # Initialize GUI.
        self.val_step = ValueWidget (None, "Step: ", QtGui.QSpinBox(), 1, 99, 1, 1,
                                     self.triggerSlicing)
        self.chk_norm = QtGui.QCheckBox ("Norm")
        self.chk_intg = QtGui.QCheckBox ("Intg")
        

        self.chk_norm.stateChanged.connect (self.triggerSlicing)
        self.chk_intg.stateChanged.connect (self.triggerSlicing)

        # attach toolbar to parent, by default (or to us, if no parent)
        for w in [ self.val_step, self.chk_norm, self.chk_intg ]:
            self.tools.addWidget(w)


    def _prepareMaster(self, wave):
        self.Axis_prepareMaster (wave)
        if self.master_wave is None:
            return
        self.val_step.spin.setRange (1, self.master_wave.shape[self.slice_axis])
            

    @QtCore.pyqtSlot()
    def slice (self, wave=None):
        '''
        Creates a slice-wise compressed version of the input wave.
        (I.e. one with a reduced size in the specified dimension,
        that has beend periodically sliced and integrated.)
        '''
        if wave is not None:
            self._prepareMaster(wave)
        if self.master_wave is None:
            return

        ax = self.slice_axis
        self.old_step = step = self.val_step.value()

        if self.chk_intg.isChecked():
            intg = step
        else:
            intg = 1
        
        self.slice_wave = ncomp(self.master_wave, ax, step, intg, self.chk_norm.isChecked())
        self.sliced.emit ([self.slice_wave])
