
import logging
log = logging.getLogger (__name__)

import sys, os, random
from PyQt4 import QtCore, QtGui

import matplotlib
matplotlib.use('Qt4Agg')
import pylab

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from paul.loader import igor
from paul.base import wave


class MatplotlibWidget(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, name=None, width=5, height=4, dpi=100,
                 bgcolor=None):
        QtGui.QWidget.__init__(self, parent)

	self.parent = parent

        self.fig = Figure(figsize=(width, height), dpi=dpi, 
                          facecolor=bgcolor, edgecolor=bgcolor)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        self.setSizePolicy(QtGui.QSizePolicy.Expanding,
                           QtGui.QSizePolicy.Expanding)
        self.updateGeometry()


    def sizeHint(self):
        w = self.fig.get_figwidth()*self.fig.get_dpi()
        h = self.fig.get_figheight()*self.fig.get_dpi()
        return QtCore.QSize(w, h)

    @QtCore.pyqtSlot()
    def clear(self):
        ''' Clears the figure. '''
        self.fig.clear()


    @QtCore.pyqtSlot()
    def reset(self):
        '''
        Clears the figure and prepares a new plot.
        The difference between this and clear() is that
        the latter only clears the figure, while this
        also prepares the canvas for the plot commands.
        '''
        self.clear()
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(False)  # We want the axes cleared every time plot() is called


    @QtCore.pyqtSlot (wave.Wave)
    def plot1d (self, data, redraw=True):
        '''
        Called to plot the specified wave. It will be passed
        to matplotlib's plot() as it is. This typically means
        that if it's a higher-D matrix, it will be plotted as a
        series of 1D graphs.
        '''
        self.axes.plot(data)
        if redraw == True:
                self.draw()


    @QtCore.pyqtSlot (wave.Wave)
    def plot2d (self, data, redraw=True):
        '''
        Called to plot the specified 2D wave. Uses matplotlib's
        imshow() to show the specified image.
        '''
        if hasattr(data, 'imgLim'):
            self.axes.imshow(data, aspect='auto', extent=data.imgLim())
        else:
            self.axes.imshow(data, aspect='auto')

        if redraw == True:
                self.draw()


    @QtCore.pyqtSlot(wave.Wave)
    def plot(self, data, redraw=True):
        '''
        Convenience wrapper for plot1d() or plot2d().
        Assuming that 'data' is one single Wave (or ndarray object),
        it calls plot1d() or plot2d(), depending on the dimensionality
        of the data.
        '''
        if not hasattr(data, 'ndim'):
            log.error ("Don't know how to plot data: %s" % data)
            return

        if data.ndim == 1:
            self.plot1d(data, redraw)
        elif data.ndim == 2:
            self.plot2d(data, redraw)
        else:
            log.error ("Don't know what to do with %d-dimensional data."
                       % (data.ndim))
        
