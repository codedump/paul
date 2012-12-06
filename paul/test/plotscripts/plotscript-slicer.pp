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

    if not len(PLOT.canvas.axes.images) > 0:
        return
    img  = PLOT.canvas.axes.images[0]
    data = PLOT.cut_wav \
      if (hasattr(PLOT, 'cut_wav') and PLOT.cut_wav is not None) \
      else PLOT.waves[0]
    
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
    GUI.cut_dim = ValueWidget (None, " 3D visual plane: ", QtGui.QSpinBox(), 0, 2, 1, 1, set3DDisplay)
    GUI.toolbar = QtGui.QToolBar("Slicing and colors")
    GUI.mainwin.addToolBarBreak()
    GUI.mainwin.addToolBar(GUI.toolbar)
    GUI.tool_act = [ GUI.toolbar.addWidget (w)
                         for w in GUI.col_min, GUI.col_max, GUI.cut_dim ]
    GUI.cut_dim_act = GUI.tool_act[2]
    
    GUI.toolbar.addAction ("+S", addSlice)
    GUI.toolbar.addAction ("+C", addWaterfall)


def exit(*args, **kwargs):
    log.debug ("EXIT")
    global GUI
    GUI.mainwin.removeToolBar (GUI.toolbar)
    for s in GUI.slices:
        log.debug ("Removing slicer %s" % s)
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

            if PLOT.waves[0].ndim == 2:
                sax = s.slice_axis
                wav = PLOT.waves[0]
                
            elif PLOT.waves[0].ndim == 3:
                wav = PLOT.cut_wav
                other_dims = range(PLOT.waves[0].ndim)
                other_dims.remove(PLOT.cut_dim)
                sax = 0 if (s.slice_axis == min(other_dims)) else 1
                #if (PLOT.cut_dim == 0):
                #    sax = 0 if s.slice_axis == 1 else 1
                #elif (PLOT.cut_dim == 1):
                #    sax = 0 if s.slice_axis == 0 else 1
                #elif (PLOT.cut_dim == 2):
                #    sax = s.slice_axis
    
            else:
                return # nothing to do, 1D wave
            
            sfrom = float(s.val_from.value())  / wav.dim[sax].size
            ssize = float(s.val_delta.value()) / wav.dim[sax].size

            sx = sfrom*(sax==1)
            sy = sfrom*(sax==0)
            sw = ssize*(sax==1) + (sax!=1) if (PLOT.cut_dim != s.slice_axis) else 0
            sh = ssize*(sax==0) + (sax!=0) if (PLOT.cut_dim != s.slice_axis) else 0
                        
            log.debug ("Slice marker rect: %f %f %f %f" % (sx, sy, sw, sh))
            s.marker_rect.set_xy ((sx, sy))
            s.marker_rect.set_width (sw)
            s.marker_rect.set_height (sh)
        else:
            print "...no marker rect"
    PLOT.canvas.draw()


@QtCore.pyqtSlot('int')
def set3DDisplay (dim):
    '''
    Takes the currently displaying waves (PLOT.waves[0]), assuming
    that it's a 3D wave, and generates a cut perpendicular to *dim*,
    in order to generate a displaying version suitable for holding
    slice markers.
    '''
    global PLOT, GUI

    log.debug ("Generating 3D display proxy plane")

    wav = PLOT.waves[0]
    ax = PLOT.canvas.axes
    
    PLOT.cut_dim = dim  # dimension perpendicular to screen
    PLOT.cut_fac = 0.5  # where to cut, relative on cut_dim
    PLOT.cut_wav = None # the representative wave

    dlim = wav.dim[PLOT.cut_dim].lim
    
    if dlim[0] < 0 and dlim[1] > 0:
        PLOT.cut_fac = wav.dim[PLOT.cut_dim].x2i(0) / wav.dim[PLOT.cut_dim].size
        log.debug ("Calculated proxy plane @ %d%% along axis %d" 
                   % (PLOT.cut_fac, PLOT.cut_dim))

    index = [slice(None)] * wav.ndim
    index[PLOT.cut_dim] = wav.dim[PLOT.cut_dim].size*PLOT.cut_fac

    PLOT.cut_wav = wav[index]
        
    ax.imshow (PLOT.cut_wav, extent=PLOT.cut_wav.imlim)
    ax.text (0, 1.05, "3D wave: representative cut at %d%% axis %d"
             % (PLOT.cut_fac*100, PLOT.cut_dim), transform=ax.transAxes)

    decorate()
    updateIndicators()
    

def populate (*args, **kwargs):

    global GUI, PLOT
    
    PLOT.waves = kwargs['wav']
    
    w0 = kwargs['wav'][0]
    wav = w0 if isinstance(w0, Wave) else w0.view(Wave)
    kwargs['can'].reset()
    ax = kwargs['can'].axes

    if wav.ndim == 1:
        log.debug ("1D data (line)")
        GUI.cut_dim_act.setEnabled(False)
        ax.plot (wav.dim[0].lim, wav)
        
    elif wav.ndim == 2:
        log.debug ("2D data (image)")
        GUI.cut_dim_act.setEnabled(False)
        ax.imshow (wav, extent=wav.imlim)
        
    elif wav.ndim == 3:
        log.debug ("3D data (volume)")
        set3DDisplay (2)
        GUI.cut_dim_act.setEnabled(True)
        kwargs['cut'] = PLOT.cut_wav
        
    decorate (*args, **kwargs)

    
def decorate(*args, **kwargs):
    '''
    Called when a newly plotted 2D image (or 1D graph) is to be decorated.
    '''
    global GUI, PLOT

    ax = PLOT.canvas.axes
    wav = PLOT.waves[0]

    if wav.ndim == 2:
        ax.set_xlabel (r'$k_{||}$ ($^\circ$)')
        ax.set_ylabel (r'E$_{total}$ (meV)')
        ax.set_ylim (wav.dim[0].lim)
        setColor(-1)
    elif PLOT.waves[0].ndim == 3:
        ax.set_xlabel (r'')
        ax.set_ylabel (r'')
        ax.set_ylim (PLOT.cut_wav.dim[0].lim)
        setColor(-1)
    
    log.debug ("Slices: %s" % str(GUI.slices))
    for s in GUI.slices:
        if (s.viewer.isVisible()):
            # calculate slice, if slicer window is visible
            log.debug ("Calling slice window %s (for axis %d)" % (s, s.slice_axis))
            s.slice(wave=wav)
                
        else:
            log.debug ("Removing hidden slice window %s (for axis %d)" % (s, s.slice_axis))
            rmSlice (s)
        


    
