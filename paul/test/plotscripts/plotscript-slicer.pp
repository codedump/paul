#!/usr/bin/python

import logging
log = logging.getLogger (__name__)

from PyQt4 import QtGui, QtCore
from paul.viewer.viewerwindow import ViewerWindow
from paul.base.wave import Wave

from paul.toolbox.slicers import SingleSlicer, CompressingSlicer
from paul.toolbox.widgets import ValueWidget
import matplotlib as mp

import math


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

    data = PLOT.waves[0]
    if not len(PLOT.canvas.axes.images) > 0:
        return
    img = PLOT.canvas.axes.images[0]
    
    dmin = data.infv('FDD', 'V_min', default=data.min())
    dmax = data.infv('FDD', 'V_max', default=data.max())

    col_min = GUI.col_min.value() * (dmax-dmin) + dmin
    col_max = GUI.col_max.value() * (dmax-dmin) + dmin

    img.set_clim (col_min, col_max)
    if i >= 0:
        PLOT.canvas.draw()


@QtCore.pyqtSlot()
def addSlice():
    global GUI, PLOT
    GUI.slices.append(SingleSlicer(axis=0, viewer=ViewerWindow()))
    GUI.slices[-1].viewer.resize (400, 350)
    GUI.slices[-1].slice(wave=PLOT.waves[0])
    GUI.slices[-1].viewer.show()
    GUI.slices[-1].sliced.connect(updateIndicators)

    # Each Slice object will carry along its own marker rectangle.
    # The marker rectangle will be drawn on the axis of the original
    # (unsliced) wave upon decorate(). It will be added to the 
    # axes in the decorate() function, because it is independent
    # of the currently selected un-sliced wave (i.e. upon selection
    # of a new wave, the slicers and thus the marker rects need
    # to survive).
    # Upon updateIndicators(), the marker rects of the respective
    # rectangles will receive their proper coordinates.
    GUI.slices[-1].marker_rect = mp.patches.Rectangle ((0, 0), 0, 0, hatch='/',
                                                       edgecolor=(0.5, 0.7, 1.0, 0.7),
                                                       lw=0.5, fill=False)

     
@QtCore.pyqtSlot()
def addWaterfall():
    global GUI, PLOT
    GUI.slices.append(CompressingSlicer(axis=0, viewer=ViewerWindow()))
    GUI.slices[-1].viewer.resize (500, 350)
    GUI.slices[-1].slice(wave=PLOT.waves[0])
    GUI.slices[-1].viewer.show()

def rmSlice(s):
    global GUI
    if hasattr(s, "marker_rect"):
        i = PLOT.canvas.axes.patches.index (s.marker_rect)
        if i > 0:
            del PLOT.canvas.axes.patches[i]
    GUI.slices.remove(s) 
    

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
    GUI.toolbar = QtGui.QToolBar("Slicing and colors")
    GUI.mainwin.addToolBarBreak()
    GUI.mainwin.addToolBar(GUI.toolbar)
    for w in GUI.col_min, GUI.col_max:
        GUI.toolbar.addWidget (w)
    GUI.toolbar.addAction ("+S", addSlice)
    GUI.toolbar.addAction ("+C", addWaterfall)


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


@QtCore.pyqtSlot()
def updateIndicators():
    global PLOT, GUI
    wav = PLOT.waves[0]

    for s in GUI.slices:
        if hasattr(s, "marker_rect"):

            # add marker rect to current axes, if not present
            if not (s.marker_rect in PLOT.canvas.axes.patches):
                s.marker_rect.set_transform (PLOT.canvas.axes.transAxes)
                s.marker_rect.set_axes (PLOT.canvas.axes)
                PLOT.canvas.axes.add_patch (s.marker_rect)
            
            sax   = s.slice_axis
            sfrom = float(s.val_from.value())  / wav.dim[sax].size
            ssize = float(s.val_delta.value()) / wav.dim[sax].size

            sx = sfrom*(sax==1)
            sy = sfrom*(sax==0)
            sw = ssize*(sax==1) + (sax!=1)
            sh = ssize*(sax==0) + (sax!=0)
            
            ## rectangle orientation depends on which axis we're slicing...
            #if sax == 1:
            #    sx = sfrom
            #    sy = 0
            #    sw = ssize
            #    sh = 1
            #else:
            #    sx = 0
            #    sy = sfrom
            #    sw = 1
            #    sh = ssize
            
            log.debug ("Slice marker rect: %f %f %f %f" % (sx, sy, sw, sh))
            s.marker_rect.set_xy ((sx, sy))
            s.marker_rect.set_width (sw)
            s.marker_rect.set_height (sh)
        else:
            print "...no marker rect"
    PLOT.canvas.draw()
    
    
def decorate(*args, **kwargs):
    '''
    Called when a newly plotted 2D image (or 1D graph) is to be decorated.
    '''
    can = kwargs['can']
    wav = kwargs['wav']
    global GUI, PLOT
    
    PLOT.waves = wav

    can.axes.set_xlabel (r'$k_{||}$ ($^\circ$)')
    can.axes.set_ylabel (r'E$_{total}$ (meV)')
    can.axes.set_ylim (PLOT.waves[0].dim[0].lim)
    setColor(-1)
    
    log.debug ("Slices: %s" % str(GUI.slices))
    for s in GUI.slices:
        if (s.viewer.isVisible()):
            # calculate slice, if slicer window is visible
            log.debug ("Calling slice window %s (for axis %d)" % (s, s.slice_axis))
            s.slice(wave=PLOT.waves[0])
                
        else:
            log.debug ("Removing hidden slice window %s (for axis %d)" % (s, s.slice_axis))
            rmSlice (s)

            #print "Patches:", PLOT.canvas.axes.patches
        


    
