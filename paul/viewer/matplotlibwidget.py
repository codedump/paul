
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

        #self.tools = NavigationToolbar (self, self.parent)
        #self.tools.show()

    def sizeHint(self):
        w = self.fig.get_figwidth()
        h = self.fig.get_figheight()
        return QtCore.QSize(w, h)

    def minimumSizeHint(self):
        return QtCore.QSize(100, 100)

    def compute_initial_figure(self):
       	pass

    # Called to plot file specified by the given full path
    @QtCore.pyqtSlot ('QString')
    def plotWave (self, data):

        self.fig.clear()
        self.axes = self.fig.add_subplot(111)
        self.axes.hold(False)  # We want the axes cleared every time
                               # plot() is called

        if (data.ndim == 1):
            log.debug ("1D plot")
            self.axes.plot(data)

        elif (data.ndim == 2):
            log.debug ("2D imshow")

            # use the correct axis scaling, if this happens to be a Wave() object.
            if hasattr(data, 'imgLim'):
                self.axes.imshow(data, aspect='auto', extent=data.imgLim())
            else:
                self.axes.imshow(data, aspect='auto')
        else:
            log.error ("Not implemented for dim > 2")
        self.draw()
